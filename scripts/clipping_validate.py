#!/usr/bin/env python3
"""Clipping per assembler: input checks, anvi'o candidate normalization and HiFi read validation.

Subcommands:
  check              Check the clipping inputs before any analysis
  normalize          Normalize anvi'o clipping candidates to exact positions
  validate           Validate candidates with clipped and spanning HiFi reads
"""
from __future__ import annotations

import argparse
import bisect
import csv
import math
import os
import re
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Set, Tuple


FIELDS = (
    "Assembler",
    "Assembly_FASTA",
    "HiFi_all_BAM",
    "HiFi_primary_BAM",
    "HiFi_primary_supplementary_BAM",
    "Output_directory",
)


def check_file(path, label):
    if not path or not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"{label} missing or empty: {path or '<unset>'}")


def check_inputs(assembly, mag_dir, manifest):
    check_file(assembly, "All-assemblies FASTA (In)")
    if assembly.suffix.lower() == ".gz":
        raise ValueError("All-assemblies FASTA must be uncompressed for samtools faidx")
    if not mag_dir or not mag_dir.is_dir():
        raise ValueError(
            f"MAG_DIR missing: {mag_dir or '<unset>'}. "
            "Use MAGverify clipping -g /path/to/MAGs"
        )
    mags = [p for p in mag_dir.iterdir() if p.is_file() and p.suffix.lower() in (".fa", ".fna", ".fasta")]
    if not mags:
        raise ValueError(f"No MAG FASTA files (.fa/.fna/.fasta) in MAG_DIR: {mag_dir}")
    check_file(manifest, "ASSEMBLER_MAPPING_MANIFEST (run 01_mapping.sh first)")
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(FIELDS) - set(reader.fieldnames or []):
            raise ValueError(
                "Invalid mapping manifest header. Expected tab-delimited columns: "
                + ", ".join(FIELDS)
            )
        seen = set()
        n = 0
        for line_no, row in enumerate(reader, start=2):
            n += 1
            assembler = (row.get("Assembler") or "").strip()
            if not assembler or not all(c.isalnum() or c in "._-" for c in assembler):
                raise ValueError(f"Invalid assembler name at manifest line {line_no}: {assembler!r}")
            if assembler in seen:
                raise ValueError(f"Duplicate assembler at manifest line {line_no}: {assembler}")
            seen.add(assembler)
            if None in row or any(not (row.get(f) or "").strip() for f in FIELDS):
                raise ValueError(f"Missing manifest field at line {line_no}. Re-run 01_mapping.sh")
            source = Path(row["Assembly_FASTA"])
            bam = Path(row["HiFi_all_BAM"])
            check_file(source, f"Assembly FASTA for {assembler}")
            check_file(bam, f"HiFi_all_BAM for {assembler} (run 01_mapping.sh first)")
            check_file(bam.parent / ".cache" / (bam.name + ".inputs.sha256"), f"HiFi_all_BAM cache for {assembler} (run 01_mapping.sh first)")
    if n == 0:
        raise ValueError("Mapping manifest has no assembler rows. Run 01_mapping.sh first")
    return len(mags), n


def _check_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--mag-dir", default="")
    parser.add_argument("--manifest", default="")
    args = parser.parse_args()
    try:
        mag_count, assembler_count = check_inputs(
            Path(args.assembly), Path(args.mag_dir) if args.mag_dir else None,
            Path(args.manifest) if args.manifest else None,
        )
    except ValueError as error:
        parser.exit(1, f"ERROR: {error}\n")
    print(f"Clipping inputs found: {mag_count} MAG FASTAs; {assembler_count} assembler(s).")
    print("File presence checked. BAM integrity and contig identities are verified during analysis.")


def cmd_check(argv):
    """Check the clipping inputs before any analysis"""
    sys.argv = [sys.argv[0] + " check", *argv]
    return _check_main()


OUTPUT_FIELDS = (
    "contig",
    "length",
    "anvio_pos_0based",
    "pos_1based",
    "relative_pos",
    "coverage",
    "clipping_count",
    "clipping_ratio",
    "default_candidate",
    "Clipping_candidate",
    "candidate_clipping",
)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anvio-clipping", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--default-output", required=True)
    parser.add_argument("--stringent-output", required=True)
    parser.add_argument("--stringent-min-depth", type=int, default=10)
    parser.add_argument("--stringent-ratio", type=float, default=1.0)
    return parser.parse_args()


def atomic_write(path: str, rows: Iterable[Mapping[str, object]]) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        delete=False,
        dir=str(destination.parent),
        prefix=destination.name + ".",
        suffix=".tmp",
    )
    try:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.close()
        os.replace(handle.name, destination)
    finally:
        try:
            handle.close()
        except Exception:
            pass
        if os.path.exists(handle.name):
            os.unlink(handle.name)


def load_fai(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                fail(f"malformed FAI row {line_number}")
            try:
                length = int(fields[1])
            except ValueError as error:
                fail(f"invalid FAI length at row {line_number}: {error}")
            if not fields[0] or length < 1 or fields[0] in lengths:
                fail(f"invalid or duplicate FAI sequence: {fields[0]}")
            lengths[fields[0]] = length
    if not lengths:
        fail("FAI contains no sequences")
    return lengths


def resolve_columns(fields: Sequence[str]) -> Dict[str, str]:
    normalized = {field.strip().lower(): field for field in fields}
    aliases = {
        "contig": ("contig",),
        "length": ("length",),
        "pos": ("pos", "position"),
        "relative_pos": ("relative_pos", "relative_position"),
        "cov": ("cov", "coverage"),
        "clipping": ("clipping", "clipping_coverage", "clipping_count"),
        "clipping_ratio": ("clipping_ratio",),
    }
    result: Dict[str, str] = {}
    for label, choices in aliases.items():
        for choice in choices:
            if choice in normalized:
                result[label] = normalized[choice]
                break
        if label not in result:
            fail(f"anvi'o clipping table lacks column for {label}")
    return result


def _normalize_main() -> None:
    args = parse_args()
    if args.stringent_min_depth < 1:
        fail("stringent minimum depth must be at least 1")
    if not 0 < args.stringent_ratio <= 1:
        fail("stringent clipping ratio must be in (0,1]")
    for path, label in ((args.anvio_clipping, "anvi'o clipping table"), (args.fai, "FAI")):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"{label} is missing or empty: {path}")

    lengths = load_fai(args.fai)
    rows: List[Dict[str, object]] = []
    seen = set()
    with open(args.anvio_clipping, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = resolve_columns(reader.fieldnames or [])
        for line_number, raw in enumerate(reader, start=2):
            contig = (raw.get(columns["contig"]) or "").strip()
            try:
                length = int(raw[columns["length"]])
                pos0 = int(raw[columns["pos"]])
                relative_pos = float(raw[columns["relative_pos"]])
                coverage = int(raw[columns["cov"]])
                clipping = int(raw[columns["clipping"]])
                ratio = float(raw[columns["clipping_ratio"]])
            except (TypeError, ValueError) as error:
                fail(f"invalid numeric value in anvi'o row {line_number}: {error}")
            if contig not in lengths or lengths[contig] != length:
                fail(f"anvi'o contig/length is stale at row {line_number}: {contig}")
            if not 0 <= pos0 < length or (contig, pos0) in seen:
                fail(f"invalid or duplicate anvi'o position: {contig}:{pos0}")
            seen.add((contig, pos0))
            if coverage < 1 or clipping < 1 or clipping > coverage:
                fail(f"invalid clipping counts at {contig}:{pos0}")
            expected_ratio = clipping / coverage
            if not math.isclose(ratio, expected_ratio, rel_tol=0, abs_tol=5.0e-7):
                fail(f"inconsistent clipping ratio at {contig}:{pos0}")
            expected_relative = pos0 / length
            if not math.isclose(relative_pos, expected_relative, rel_tol=0, abs_tol=5.0e-7):
                fail(f"inconsistent relative position at {contig}:{pos0}")
            stringent = coverage >= args.stringent_min_depth and math.isclose(
                ratio, args.stringent_ratio, rel_tol=0, abs_tol=5.0e-7
            )
            rows.append(
                {
                    "contig": contig,
                    "length": length,
                    "anvio_pos_0based": pos0,
                    "pos_1based": pos0 + 1,
                    "relative_pos": f"{relative_pos:.12g}",
                    "coverage": coverage,
                    "clipping_count": clipping,
                    "clipping_ratio": f"{ratio:.6f}",
                    "default_candidate": "Y",
                    "Clipping_candidate": "Y" if stringent else "N",
                    "candidate_clipping": "Y",
                }
            )

    rows.sort(key=lambda row: (str(row["contig"]), int(row["anvio_pos_0based"])))
    stringent_rows = [row for row in rows if row["Clipping_candidate"] == "Y"]
    atomic_write(args.default_output, rows)
    atomic_write(args.stringent_output, stringent_rows)
    print(f"anvi'o default candidates: {len(rows)}")
    print(f"Stringent candidates: {len(stringent_rows)}")


def cmd_normalize(argv):
    """Normalize anvi'o clipping candidates to exact positions"""
    sys.argv = [sys.argv[0] + " normalize", *argv]
    return _normalize_main()


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

MATCH_OPS = {"M", "=", "X"}

REF_OPS = {"M", "D", "N", "=", "X"}

CLIP_OPS = {"S", "H"}

FINAL_MIN_CLIPPED_QNAMES = 2


def _validate_parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--high-confidence-min-coverage", type=int, default=10)
    parser.add_argument("--high-confidence-ratio", type=float, default=1.0)
    return parser.parse_args()


def cigar_ops(cigar: str) -> List[Tuple[int, str]]:
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not operations or "".join(f"{length}{op}" for length, op in operations) != cigar:
        fail(f"invalid CIGAR: {cigar}")
    return operations


def terminal_clip(operations: Sequence[Tuple[int, str]], left: bool) -> int:
    values = operations if left else reversed(operations)
    total = 0
    for length, op in values:
        if op in CLIP_OPS:
            total += length
        else:
            break
    return total


def matched_runs(ref_start0: int, operations: Sequence[Tuple[int, str]]) -> List[Tuple[int, int]]:
    """Return intervals containing actual M/= /X aligned reference bases."""
    runs: List[Tuple[int, int]] = []
    ref_pos = ref_start0
    run_start = None
    run_end = None
    for length, op in operations:
        if op in MATCH_OPS:
            start, end = ref_pos, ref_pos + length
            if run_start is None:
                run_start, run_end = start, end
            elif run_end == start:
                run_end = end
            else:
                runs.append((run_start, run_end))
                run_start, run_end = start, end
            ref_pos = end
        elif op in {"I", "P"}:
            continue
        elif op in {"D", "N"}:
            if run_start is not None and run_end is not None:
                runs.append((run_start, run_end))
                run_start = run_end = None
            ref_pos += length
        elif op in {"S", "H"}:
            if run_start is not None and run_end is not None:
                runs.append((run_start, run_end))
                run_start = run_end = None
        else:
            fail(f"unsupported CIGAR operation: {op}")
    if run_start is not None and run_end is not None:
        runs.append((run_start, run_end))
    return runs


def candidate_idn_event_lengths(
    ref_start0: int, operations: Sequence[Tuple[int, str]], pos0: int
) -> Dict[str, List[int]]:
    """Return candidate-local I/D/N event lengths for audit only."""
    events: Dict[str, List[int]] = defaultdict(list)
    ref_pos = ref_start0
    for length, op in operations:
        if op in MATCH_OPS:
            ref_pos += length
        elif op == "I":
            if ref_pos in {pos0, pos0 + 1}:
                events["I"].append(length)
        elif op in {"D", "N"}:
            if ref_pos <= pos0 < ref_pos + length:
                events[op].append(length)
            ref_pos += length
        elif op in {"S", "H", "P"}:
            continue
        else:
            fail(f"unsupported CIGAR operation: {op}")
    return events


def matched_bp_each_side(
    runs: Sequence[Tuple[int, int]], pos0: int
) -> Tuple[int, int]:
    """Count actual M/= /X aligned reference bases on each side of a candidate."""
    left = 0
    right = 0
    for start, end in runs:
        if start < pos0:
            left += max(0, min(end, pos0) - start)
        if end > pos0 + 1:
            right += max(0, end - max(start, pos0 + 1))
    return left, right


def update_qname_max(
    store: Dict[int, Dict[str, int]], row_index: int, qname: str, value: int
) -> None:
    """Retain the largest audit value observed for one QNAME at one site."""
    previous = store[row_index].get(qname)
    if previous is None or value > previous:
        store[row_index][qname] = value


def length_stats(values: Iterable[int]) -> Tuple[str, str, str]:
    """Return min/median/max strings, or NA when no values are present."""
    ordered = sorted(values)
    if not ordered:
        return "NA", "NA", "NA"
    med = statistics.median(ordered)
    med_text = str(int(med)) if float(med).is_integer() else f"{med:.1f}"
    return str(ordered[0]), med_text, str(ordered[-1])


def stream_sam(command: Sequence[str]) -> Iterator[List[str]]:
    process = subprocess.Popen(
        list(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="utf-8", errors="replace",
    )
    assert process.stdout is not None and process.stderr is not None
    for raw in process.stdout:
        fields = raw.rstrip("\n\r").split("\t")
        if len(fields) >= 11:
            yield fields
    process.stdout.close()
    stderr = process.stderr.read()
    code = process.wait()
    if code != 0:
        fail(f"command failed ({code}): {' '.join(command)}\n{stderr}")


def _validate_atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False,
        dir=str(destination.parent), prefix=destination.name + ".", suffix=".tmp",
    )
    try:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.close()
        os.replace(handle.name, destination)
    finally:
        try:
            handle.close()
        except Exception:
            pass
        if os.path.exists(handle.name):
            os.unlink(handle.name)


def candidate_indices_in_run(
    positions: Sequence[int], indices: Sequence[int], start: int, end: int
) -> Iterator[int]:
    lo = bisect.bisect_left(positions, start)
    hi = bisect.bisect_left(positions, end)
    for offset in range(lo, hi):
        yield indices[offset]


def _validate_main() -> None:
    args = _validate_parse_args()
    if args.threads < 1:
        fail("invalid thread count")
    if args.high_confidence_min_coverage < 1:
        fail("high-confidence minimum coverage must be at least 1")
    if not 0 < args.high_confidence_ratio <= 1:
        fail("high-confidence clipping ratio must be in (0,1]")
    for path, label in ((args.candidates, "candidate table"), (args.bam, "BAM")):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"{label} is missing or empty: {path}")

    with open(args.candidates, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        base_fields = list(reader.fieldnames or [])
        required = {
            "contig", "length", "anvio_pos_0based", "pos_1based", "coverage",
            "clipping_count", "clipping_ratio", "default_candidate",
            "Clipping_candidate",
        }
        missing = required - set(base_fields)
        if missing:
            fail(f"candidate table lacks columns: {sorted(missing)}")
        rows = [dict(row) for row in reader]

    by_contig: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    lengths: Dict[str, int] = {}
    seen = set()
    for index, row in enumerate(rows):
        contig = row["contig"].strip()
        try:
            length = int(row["length"])
            pos0 = int(row["anvio_pos_0based"])
        except ValueError as error:
            fail(f"invalid candidate coordinate: {error}")
        if not contig or length < 1 or not 0 <= pos0 < length or (contig, pos0) in seen:
            fail(f"invalid or duplicate candidate: {contig}:{pos0}")
        seen.add((contig, pos0))
        if contig in lengths and lengths[contig] != length:
            fail(f"inconsistent candidate contig length: {contig}")
        lengths[contig] = length
        by_contig[contig].append((pos0, index))

    for contig in by_contig:
        by_contig[contig].sort()
    positions_by_contig = {
        contig: [value[0] for value in values] for contig, values in by_contig.items()
    }
    indices_by_contig = {
        contig: [value[1] for value in values] for contig, values in by_contig.items()
    }
    index_by_position = {
        contig: {pos: index for pos, index in values} for contig, values in by_contig.items()
    }

    overlapping_qnames: Dict[int, Set[str]] = defaultdict(set)
    primary_clipped: Dict[int, Set[str]] = defaultdict(set)
    secondary_clipped: Dict[int, Set[str]] = defaultdict(set)
    supplementary_clipped: Dict[int, Set[str]] = defaultdict(set)
    spanning: Dict[int, Set[str]] = defaultdict(set)
    spanning_with_i: Dict[int, Set[str]] = defaultdict(set)
    spanning_with_d: Dict[int, Set[str]] = defaultdict(set)
    spanning_with_n: Dict[int, Set[str]] = defaultdict(set)
    spanning_min_flank_by_qname: Dict[int, Dict[str, int]] = defaultdict(dict)
    spanning_i_length_by_qname: Dict[int, Dict[str, int]] = defaultdict(dict)
    spanning_d_length_by_qname: Dict[int, Dict[str, int]] = defaultdict(dict)
    spanning_n_length_by_qname: Dict[int, Dict[str, int]] = defaultdict(dict)
    recomputed_depth: Dict[int, int] = defaultdict(int)
    recomputed_clipping: Dict[int, int] = defaultdict(int)
    primary_records: Dict[int, int] = defaultdict(int)
    secondary_records: Dict[int, int] = defaultdict(int)
    supplementary_records: Dict[int, int] = defaultdict(int)

    if rows:
        command = ["samtools", "view", "-@", str(args.threads), "-F", "4", args.bam]
        for fields in stream_sam(command):
            qname, contig, cigar = fields[0], fields[2], fields[5]
            if contig not in by_contig or cigar == "*":
                continue
            try:
                flag = int(fields[1])
                ref_start0 = int(fields[3]) - 1
            except ValueError as error:
                fail(f"invalid SAM numeric field: {error}")
            if flag & 0x4:
                continue
            if ref_start0 < 0:
                continue
            if flag & 0x100:
                alignment_type = "secondary"
            elif flag & 0x800:
                alignment_type = "supplementary"
            else:
                alignment_type = "primary"
            operations = cigar_ops(cigar)
            ref_end0 = ref_start0 + sum(length for length, op in operations if op in REF_OPS)
            if ref_end0 <= ref_start0 or ref_end0 > lengths[contig]:
                fail(f"invalid alignment span: {contig}:{ref_start0}-{ref_end0}")

            positions = positions_by_contig[contig]
            indices = indices_by_contig[contig]
            runs = matched_runs(ref_start0, operations)
            touched_indices: Set[int] = set()
            for run_start, run_end in runs:
                for row_index in candidate_indices_in_run(positions, indices, run_start, run_end):
                    overlapping_qnames[row_index].add(qname)
                    touched_indices.add(row_index)
            broad_start = ref_start0 + 1
            broad_end = ref_end0 - 1
            if broad_start < broad_end:
                for row_index in candidate_indices_in_run(positions, indices, broad_start, broad_end):
                    pos0 = int(rows[row_index]["anvio_pos_0based"])
                    has_left = any(run_start < pos0 for run_start, _run_end in runs)
                    has_right = any(run_end > pos0 + 1 for _run_start, run_end in runs)
                    if has_left and has_right:
                        spanning[row_index].add(qname)
                        left_bp, right_bp = matched_bp_each_side(runs, pos0)
                        update_qname_max(
                            spanning_min_flank_by_qname, row_index, qname,
                            min(left_bp, right_bp),
                        )
                        events = candidate_idn_event_lengths(ref_start0, operations, pos0)
                        if events.get("I"):
                            spanning_with_i[row_index].add(qname)
                            update_qname_max(
                                spanning_i_length_by_qname, row_index, qname, max(events["I"])
                            )
                        if events.get("D"):
                            spanning_with_d[row_index].add(qname)
                            update_qname_max(
                                spanning_d_length_by_qname, row_index, qname, max(events["D"])
                            )
                        if events.get("N"):
                            spanning_with_n[row_index].add(qname)
                            update_qname_max(
                                spanning_n_length_by_qname, row_index, qname, max(events["N"])
                            )

            for row_index in touched_indices:
                recomputed_depth[row_index] += 1
                if alignment_type == "primary":
                    primary_records[row_index] += 1
                elif alignment_type == "secondary":
                    secondary_records[row_index] += 1
                else:
                    supplementary_records[row_index] += 1

            clipped_positions: Set[int] = set()
            if terminal_clip(operations, left=True) > 0 and ref_start0 != 0:
                clipped_positions.add(ref_start0)
            if terminal_clip(operations, left=False) > 0 and ref_end0 != lengths[contig]:
                clipped_positions.add(ref_end0 - 1)
            position_to_index = index_by_position[contig]
            for pos0 in clipped_positions:
                row_index = position_to_index.get(pos0)
                if row_index is None:
                    continue
                recomputed_clipping[row_index] += 1
                if alignment_type == "primary":
                    primary_clipped[row_index].add(qname)
                elif alignment_type == "secondary":
                    secondary_clipped[row_index].add(qname)
                else:
                    supplementary_clipped[row_index].add(qname)

    appended = [
        "n_primary_alignment_records_covering_site",
        "n_secondary_alignment_records_covering_site",
        "n_supplementary_alignment_records_covering_site",
        "n_unique_all_overlapping_qnames",
        "n_unique_primary_clipped_qnames",
        "n_unique_secondary_clipped_qnames",
        "n_unique_supplementary_clipped_qnames",
        "n_unique_all_clipped_qnames",
        "qname_clipping_ratio",
        "n_unique_exact_spanning_qnames",
        "exact_spanning_min_flank_bp_min",
        "exact_spanning_min_flank_bp_median",
        "exact_spanning_min_flank_bp_max",
        "n_unique_spanning_with_I_qnames",
        "spanning_I_length_bp_min",
        "spanning_I_length_bp_median",
        "spanning_I_length_bp_max",
        "n_unique_spanning_with_D_qnames",
        "spanning_D_length_bp_min",
        "spanning_D_length_bp_median",
        "spanning_D_length_bp_max",
        "n_unique_spanning_with_N_qnames",
        "spanning_N_length_bp_min",
        "spanning_N_length_bp_median",
        "spanning_N_length_bp_max",
        "n_unique_spanning_with_any_IDN_qnames",
        "n_unique_clip_and_span_qnames",
        "valid_spanning_reads",
        "anvio_count_check",
        "Mixed_clipping",
        "exact_position_validation",
        "clipping_confidence",
        "support_class",
        "evidence_class",
        "validated_clipping",
        "Clipping_analyzable",
        "Clipping_analysis_reason",
        "Clipping_ge1",
        "Clipping",
        "High_confidence_clipping",
    ]
    output_fields = list(base_fields)
    for field in appended:
        if field not in output_fields:
            output_fields.append(field)

    output_rows: List[Dict[str, str]] = []
    for index, row in enumerate(rows):
        primary = primary_clipped.get(index, set())
        secondary = secondary_clipped.get(index, set())
        supplementary = supplementary_clipped.get(index, set())
        clipped = primary | secondary | supplementary
        spanning_reads = spanning.get(index, set())
        spanning_i_reads = spanning_with_i.get(index, set())
        spanning_d_reads = spanning_with_d.get(index, set())
        spanning_n_reads = spanning_with_n.get(index, set())
        spanning_any_idn_reads = spanning_i_reads | spanning_d_reads | spanning_n_reads
        span_flank_min, span_flank_median, span_flank_max = length_stats(
            spanning_min_flank_by_qname.get(index, {}).values()
        )
        i_len_min, i_len_median, i_len_max = length_stats(
            spanning_i_length_by_qname.get(index, {}).values()
        )
        d_len_min, d_len_median, d_len_max = length_stats(
            spanning_d_length_by_qname.get(index, {}).values()
        )
        n_len_min, n_len_median, n_len_max = length_stats(
            spanning_n_length_by_qname.get(index, {}).values()
        )
        clip_and_span = clipped & spanning_reads
        qname_depth = len(overlapping_qnames.get(index, set()))
        qname_ratio = len(clipped) / qname_depth if qname_depth else 0.0
        raw_depth = int(row["coverage"])
        raw_clipping = int(row["clipping_count"])
        raw_ratio = float(row["clipping_ratio"])
        stringent = row["Clipping_candidate"].strip().upper() == "Y"
        high_confidence = raw_depth >= args.high_confidence_min_coverage and math.isclose(
            raw_ratio, args.high_confidence_ratio, rel_tol=0, abs_tol=5.0e-7
        )
        if stringent != high_confidence:
            fail(
                "Clipping_candidate is inconsistent with the configured high-confidence "
                f"thresholds at {row['contig']}:{row['pos_1based']}"
            )

        reference_boundary = int(row["anvio_pos_0based"]) in {0, int(row["length"]) - 1}
        mixed = bool(clipped and spanning_reads)
        permissive_ge1 = False
        if reference_boundary:
            validated = False
            exact_validation = "BOUNDARY_UNRESOLVED"
            confidence = "NOT_VALIDATED"
            evidence = "REFERENCE_BOUNDARY_UNRESOLVED"
            support_class = "BOUNDARY_UNRESOLVED"
        else:
            permissive_ge1 = len(clipped) >= 1 and not spanning_reads
            validated = len(clipped) >= FINAL_MIN_CLIPPED_QNAMES and not spanning_reads
            if validated:
                exact_validation = "PASS"
            elif permissive_ge1:
                exact_validation = "INSUFFICIENT_CLIPPED_QNAMES"
            else:
                exact_validation = "FAIL"
            if validated and high_confidence:
                confidence = "HIGH_CONFIDENCE"
                evidence = "VALIDATED_HIGH_CONFIDENCE"
                support_class = "VALIDATED_HIGH_CONFIDENCE"
            elif validated:
                confidence = "DEFAULT_ONLY"
                evidence = "VALIDATED_DEFAULT_ONLY"
                support_class = "VALIDATED_DEFAULT_ONLY"
            elif permissive_ge1:
                confidence = "NOT_VALIDATED"
                evidence = "INSUFFICIENT_CLIPPED_QNAMES"
                support_class = "BELOW_MIN_CLIPPED_READS"
            elif spanning_reads:
                confidence = "NOT_VALIDATED"
                evidence = "EXACT_SITE_SPANNED"
                support_class = "SPANNED"
            else:
                confidence = "NOT_VALIDATED"
                evidence = "NO_CLIPPED_QNAME_AT_EXACT_SITE"
                support_class = "UNRESOLVED"

        if reference_boundary:
            count_check = "BOUNDARY_NOT_COMPARABLE"
        else:
            count_check = (
                "PASS"
                if recomputed_depth.get(index, 0) == raw_depth
                and recomputed_clipping.get(index, 0) == raw_clipping
                else "MISMATCH"
            )
        row.update(
            {
                "n_primary_alignment_records_covering_site": str(primary_records.get(index, 0)),
                "n_secondary_alignment_records_covering_site": str(secondary_records.get(index, 0)),
                "n_supplementary_alignment_records_covering_site": str(supplementary_records.get(index, 0)),
                "n_unique_all_overlapping_qnames": str(qname_depth),
                "n_unique_primary_clipped_qnames": str(len(primary)),
                "n_unique_secondary_clipped_qnames": str(len(secondary)),
                "n_unique_supplementary_clipped_qnames": str(len(supplementary)),
                "n_unique_all_clipped_qnames": str(len(clipped)),
                "qname_clipping_ratio": f"{qname_ratio:.6f}",
                "n_unique_exact_spanning_qnames": str(len(spanning_reads)),
                "exact_spanning_min_flank_bp_min": span_flank_min,
                "exact_spanning_min_flank_bp_median": span_flank_median,
                "exact_spanning_min_flank_bp_max": span_flank_max,
                "n_unique_spanning_with_I_qnames": str(len(spanning_i_reads)),
                "spanning_I_length_bp_min": i_len_min,
                "spanning_I_length_bp_median": i_len_median,
                "spanning_I_length_bp_max": i_len_max,
                "n_unique_spanning_with_D_qnames": str(len(spanning_d_reads)),
                "spanning_D_length_bp_min": d_len_min,
                "spanning_D_length_bp_median": d_len_median,
                "spanning_D_length_bp_max": d_len_max,
                "n_unique_spanning_with_N_qnames": str(len(spanning_n_reads)),
                "spanning_N_length_bp_min": n_len_min,
                "spanning_N_length_bp_median": n_len_median,
                "spanning_N_length_bp_max": n_len_max,
                "n_unique_spanning_with_any_IDN_qnames": str(len(spanning_any_idn_reads)),
                "n_unique_clip_and_span_qnames": str(len(clip_and_span)),
                "valid_spanning_reads": str(len(spanning_reads)),
                "anvio_count_check": count_check,
                "Mixed_clipping": "Y" if mixed else "N",
                "exact_position_validation": exact_validation,
                "clipping_confidence": confidence,
                "support_class": support_class,
                "evidence_class": evidence,
                "validated_clipping": "Y" if validated else "N",
                "Clipping_analyzable": "N" if reference_boundary else "Y",
                "Clipping_analysis_reason": (
                    "REFERENCE_BOUNDARY_TWO_SIDED_TEST_IMPOSSIBLE"
                    if reference_boundary else "EXACT_POSITION_EVALUATED"
                ),
                "Clipping_ge1": "Y" if permissive_ge1 else "N",
                "Clipping": "NA" if reference_boundary else ("Y" if validated else "N"),
                "High_confidence_clipping": "NA" if reference_boundary else ("Y" if validated and high_confidence else "N"),
            }
        )
        output_rows.append(row)

    _validate_atomic_write(args.output, output_fields, output_rows)
    print(f"Read-level candidates annotated: {len(output_rows)}")
    print(f"Clipping_ge1=Y audit sites: {sum(row['Clipping_ge1'] == 'Y' for row in output_rows)}")
    print(f"Exact-position validation PASS (>=2 clipped reads): {sum(row['exact_position_validation'] == 'PASS' for row in output_rows)}")
    print(f"Default-only validated sites: {sum(row['clipping_confidence'] == 'DEFAULT_ONLY' for row in output_rows)}")
    print(f"Clipping=Y sites (>=2 clipped reads): {sum(row['Clipping'] == 'Y' for row in output_rows)}")
    print(f"High_confidence_clipping=Y sites: {sum(row['High_confidence_clipping'] == 'Y' for row in output_rows)}")
    print(f"Candidates with spanning counter-evidence: {sum(int(row['n_unique_exact_spanning_qnames']) > 0 for row in output_rows)}")
    print(f"Spanning candidates with candidate I audit: {sum(int(row['n_unique_spanning_with_I_qnames']) > 0 for row in output_rows)}")
    print(f"Spanning candidates with candidate D audit: {sum(int(row['n_unique_spanning_with_D_qnames']) > 0 for row in output_rows)}")
    print(f"Spanning candidates with candidate N audit: {sum(int(row['n_unique_spanning_with_N_qnames']) > 0 for row in output_rows)}")
    mismatches = sum(row["anvio_count_check"] == "MISMATCH" for row in output_rows)
    print(f"anvi'o count-check mismatches: {mismatches}")


def cmd_validate(argv):
    """Validate candidates with clipped and spanning HiFi reads"""
    sys.argv = [sys.argv[0] + " validate", *argv]
    return _validate_main()


COMMANDS = {
    'check': cmd_check,
    'normalize': cmd_normalize,
    'validate': cmd_validate,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
