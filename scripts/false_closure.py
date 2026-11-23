#!/usr/bin/env python3
"""False closure: contig circularity, breakpoint detection and contig/MAG summaries.

Subcommands:
  resolve            Resolve which contigs each assembler reported as circular
  detect             Detect false-closure breakpoints from validated clipping and split HiFi alignments
  mag                Write the MAG table and its audit
  contig             Write the contig table
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from typing import (
    Dict,
    Iterable,
    Iterator,
    List,
    Mapping,
    NoReturn,
    Optional,
    Sequence,
    Set,
    Tuple,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from circularity_finalize import resolve_call_table  # noqa: E402


RIBOSOMAL_PATTERN = re.compile(r"^ribosom", re.IGNORECASE)

SCG_SOURCES = ("Bacteria_71", "Archaea_76")


def _source_and_family(row: Dict[str, str]) -> Tuple[str, str]:
    marker = (row.get("hmm_name") or row.get("marker") or "").strip()
    source = (row.get("HMM_collection") or row.get("source") or row.get("hmm_source") or "").strip()
    if "|" in marker:
        prefix, rest = marker.split("|", 1)
        if prefix in SCG_SOURCES:
            source, marker = prefix, rest
    return source, marker


def load_ribosomal_families(path: str) -> Dict[str, Tuple[str, ...]]:
    """Return contig -> sorted distinct ribosomal-protein HMM families."""
    per_source: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if "contig" not in fields or not ({"hmm_name", "marker"} & fields):
            raise SystemExit(f"ERROR: SCG table needs contig and hmm_name columns: {path}")
        for row in reader:
            contig = (row.get("contig") or "").strip()
            source, family = _source_and_family(row)
            if contig and family and RIBOSOMAL_PATTERN.match(family):
                per_source[contig][source or "NA"].add(family)

    families: Dict[str, Tuple[str, ...]] = {}
    for contig, by_source in per_source.items():
        best = max(by_source.items(), key=lambda item: (len(item[1]), item[0]))[1]
        families[contig] = tuple(sorted(best))
    return families


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class Alignment:
    qname: str
    flag: int
    rname: str
    mapq: int
    cigar: str
    rstart: int
    rend: int
    left_clip: int
    right_clip: int
    matched_ref_bases: int

    @property
    def secondary(self) -> bool:
        return bool(self.flag & 256)

    @property
    def supplementary(self) -> bool:
        return bool(self.flag & 2048)


def fail(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


def atomic_tsv(path: str, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for line_no, row in enumerate(rows, start=2):
                if len(row) != len(header):
                    fail(
                        f"TSV schema mismatch at {path}:{line_no}: "
                        f"{len(row)} != {len(header)}"
                    )
                writer.writerow(row)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def stream_sam(command: Sequence[str]) -> Iterator[str]:
    process = subprocess.Popen(
        list(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="utf-8", errors="replace",
    )
    assert process.stdout is not None and process.stderr is not None
    for line in process.stdout:
        yield line
    process.stdout.close()
    stderr = process.stderr.read()
    rc = process.wait()
    if rc != 0:
        fail(f"command failed ({rc}): {' '.join(command)}\n{stderr}")


def cigar_operations(cigar: str) -> List[Tuple[int, str]]:
    ops = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if not ops or "".join(f"{n}{op}" for n, op in ops) != cigar:
        fail(f"unsupported or malformed CIGAR: {cigar}")
    return ops


def parse_alignment(line: str) -> Optional[Alignment]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11 or fields[2] == "*" or fields[5] == "*":
        return None
    flag = int(fields[1])
    ops = cigar_operations(fields[5])
    ref_span = sum(n for n, op in ops if op in "MDN=X")
    matched = sum(n for n, op in ops if op in "M=X")
    if ref_span <= 0 or matched <= 0:
        return None
    left_clip = 0
    for n, op in ops:
        if op not in "SH":
            break
        left_clip += n
    right_clip = 0
    for n, op in reversed(ops):
        if op not in "SH":
            break
        right_clip += n
    rstart = int(fields[3]) - 1
    rend = rstart + ref_span
    return Alignment(
        fields[0], flag, fields[2], int(fields[4]), fields[5],
        rstart, rend, left_clip, right_clip, matched,
    )


def representative(alignments: Sequence[Alignment]) -> Optional[Alignment]:
    if not alignments:
        return None
    return sorted(
        alignments,
        key=lambda a: (
            a.secondary,
            -a.matched_ref_bases,
            a.supplementary,
            a.rstart,
            a.rend,
            a.flag,
            a.cigar,
        ),
    )[0]


def load_alignments_for_qnames(
    bam: str, qnames: Set[str]
) -> Dict[str, List[Alignment]]:
    result: Dict[str, List[Alignment]] = defaultdict(list)
    if not qnames:
        return result
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".qnames") as qf:
        qpath = qf.name
        for qname in sorted(qnames):
            qf.write(qname + "\n")
    try:
        for line in stream_sam(
            ["samtools", "view", "-F", "4", "-q", "0", "-N", qpath, bam]
        ):
            aln = parse_alignment(line)
            if aln is not None:
                result[aln.qname].append(aln)
    finally:
        os.unlink(qpath)
    return result


def summarize_external_connection(
    source: Dict[str, Alignment],
    all_by_read: Dict[str, List[Alignment]],
    source_contig: str,
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Return primary/supplementary QNAMEs supporting cross-contig split mapping."""
    supported: Set[str] = set()
    targets: Dict[str, Set[str]] = defaultdict(set)
    for qname, source_aln in source.items():
        if source_aln.secondary:
            continue
        for other in all_by_read.get(qname, []):
            if other.rname == source_contig or other.secondary:
                continue
            if not (source_aln.supplementary or other.supplementary):
                continue
            supported.add(qname)
            targets[other.rname].add(qname)
    return supported, targets


def target_count_text(targets: Dict[str, Set[str]]) -> str:
    if not targets:
        return "NA"
    return ";".join(
        f"{target}:{len(reads)}"
        for target, reads in sorted(
            targets.items(), key=lambda item: (-len(item[1]), item[0])
        )
    )


def dominant_target(
    targets: Dict[str, Set[str]], denominator: int
) -> Tuple[str, int, float]:
    """Audit-only summary; dominance is never a binary evidence gate."""
    if not targets:
        return "NA", 0, 0.0
    ranked = sorted(
        ((target, len(reads)) for target, reads in targets.items()),
        key=lambda item: (-item[1], item[0]),
    )
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        n = ranked[0][1]
        return "AMBIGUOUS", n, n / denominator if denominator else 0.0
    target, n = ranked[0]
    return target, n, n / denominator if denominator else 0.0


def load_contig_to_mag(path: str) -> Dict[str, str]:
    with open(path, encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        lower = {field.lower(): field for field in fields}
        mag_col = lower.get("mag_id") or lower.get("id")
        contig_col = lower.get("contig") or lower.get("original_contig")
        if not mag_col or not contig_col:
            fail(f"MAG contig map lacks MAG_ID/Contig columns: {path}")
        return {
            row[contig_col].strip(): row[mag_col].strip()
            for row in reader
            if row.get(contig_col, "").strip()
        }


VERSION = "resolve_false_closure_contig_circularity"

HIFIASM_META_RE = re.compile(r"^s[0-9]+\.[uc]tg[0-9]{6}([lc])$")


def norm_assembler(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.casefold())


def normalize_call(value: str) -> str | None:
    v = (value or "").strip().casefold()
    if v in {"y", "yes", "true", "1", "+", "c", "circular"}:
        return "Y"
    if v in {"n", "no", "false", "0", "-", "l", "linear"}:
        return "N"
    return None


def load_mag_map(path: str) -> Tuple[List[Tuple[str, str]], Dict[str, int]]:
    rows: List[Tuple[str, str]] = []
    counts: Dict[str, int] = defaultdict(int)
    seen: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        lower = {x.casefold(): x for x in fields}
        mag_col = lower.get("mag_id") or lower.get("id")
        contig_col = lower.get("contig") or lower.get("original_contig")
        if not mag_col or not contig_col:
            fail("MAG contig map requires MAG_ID/ID and Contig columns")
        for line_no, row in enumerate(reader, start=2):
            mag = (row.get(mag_col) or "").strip()
            contig = (row.get(contig_col) or "").strip()
            if not mag or not contig:
                fail(f"blank MAG/contig at {path}:{line_no}")
            prev = seen.get(contig)
            if prev is not None and prev != mag:
                fail(f"contig assigned to multiple MAGs: {contig}")
            if prev is None:
                seen[contig] = mag
                rows.append((mag, contig))
                counts[mag] += 1
    if not rows:
        fail(f"MAG contig map contains no rows: {path}")
    return rows, counts


def fasta_headers(path: str) -> Iterable[Tuple[str, str]]:
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\r\n")
            name = header.split(None, 1)[0]
            if name:
                yield name, header


def parse_header_circularity(path: str, assembler: str) -> Dict[str, Tuple[str, str]]:
    out: Dict[str, Tuple[str, str]] = {}
    key = norm_assembler(assembler)
    for contig, header in fasta_headers(path):
        match = re.search(r"(?:^|\s)circular\s*=\s*([^\s]+)", header, flags=re.I)
        if match:
            call = normalize_call(match.group(1))
            if call is not None:
                out[contig] = (call, "SOURCE_FASTA_HEADER_CIRCULAR_FIELD")
                continue
        if "hifiasm" in key:
            match = HIFIASM_META_RE.fullmatch(contig)
            if match:
                out[contig] = (
                    "Y" if match.group(1) == "c" else "N",
                    "HIFIASM_META_CONTIG_SUFFIX_C_L",
                )
    return out


def parse_flye_info(path: str) -> Dict[str, Tuple[str, str]]:
    out: Dict[str, Tuple[str, str]] = {}
    with open(path, encoding="utf-8-sig", errors="replace") as handle:
        header = None
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            tokens = line.lstrip("#").split()
            if header is None:
                lower = [re.sub(r"[^a-z0-9]+", "", x.casefold()) for x in tokens]
                if any(x in {"seqname", "sequence", "contig", "contigid"} for x in lower) and any(
                    x in {"circ", "circular", "iscircular"} for x in lower
                ):
                    header = lower
                    continue
                if not line.startswith("#") and len(tokens) >= 4:
                    call = normalize_call(tokens[3])
                    if call is not None:
                        out[tokens[0]] = (call, f"FLYE_ASSEMBLY_INFO:{os.path.abspath(path)}")
                continue
            if line.startswith("#"):
                continue
            if len(tokens) < len(header):
                continue
            try:
                name_idx = next(i for i, x in enumerate(header) if x in {"seqname", "sequence", "contig", "contigid"})
                circ_idx = next(i for i, x in enumerate(header) if x in {"circ", "circular", "iscircular"})
            except StopIteration:
                continue
            call = normalize_call(tokens[circ_idx])
            if call is not None:
                out[tokens[name_idx]] = (call, f"FLYE_ASSEMBLY_INFO:{os.path.abspath(path)}")
    return out


def flye_candidate_paths(fasta: str, explicit: str) -> List[str]:
    paths: List[str] = []
    if explicit:
        paths.append(explicit)
    directory = os.path.dirname(os.path.abspath(fasta))
    stem = os.path.splitext(os.path.basename(fasta))[0]
    for name in (
        "assembly_info.txt",
        f"{stem}.assembly_info.txt",
        f"{stem}_assembly_info.txt",
        "metaflye.assembly_info.txt",
        "metaFlye.assembly_info.txt",
    ):
        paths.append(os.path.join(directory, name))
    seen = set()
    return [p for p in paths if not (p in seen or seen.add(p))]


def load_source_calls(mapping_manifest: str, flye_info: str) -> Tuple[Dict[str, str], Dict[str, Tuple[str, str]], Dict[str, str]]:
    source_of: Dict[str, str] = {}
    calls: Dict[str, Tuple[str, str]] = {}
    fasta_of: Dict[str, str] = {}
    with open(mapping_manifest, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Assembler", "Assembly_FASTA"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"mapping manifest lacks columns: {sorted(missing)}")
        for row in reader:
            assembler = (row.get("Assembler") or "").strip()
            fasta = (row.get("Assembly_FASTA") or "").strip()
            if not assembler or not fasta:
                continue
            if not os.path.isfile(fasta):
                fail(f"assembly FASTA missing for {assembler}: {fasta}")
            fasta_of[assembler] = fasta
            header_calls = parse_header_circularity(fasta, assembler)
            names = []
            for contig, _header in fasta_headers(fasta):
                names.append(contig)
                prev = source_of.get(contig)
                if prev is not None and prev != assembler:
                    fail(f"contig occurs in multiple assemblers: {contig} ({prev}, {assembler})")
                source_of[contig] = assembler
            for contig, value in header_calls.items():
                calls[contig] = value

            if "flye" in norm_assembler(assembler):
                name_set = set(names)
                candidates = []
                for candidate in flye_candidate_paths(fasta, flye_info):
                    if os.path.isfile(candidate) and os.path.getsize(candidate) > 0:
                        parsed = parse_flye_info(candidate)
                        overlap = len(name_set.intersection(parsed))
                        if overlap:
                            candidates.append((overlap, candidate, parsed))
                if candidates:
                    candidates.sort(key=lambda x: (-x[0], x[1]))
                    best_overlap = candidates[0][0]
                    best = [x for x in candidates if x[0] == best_overlap]
                    if len(best) > 1 and len({os.path.realpath(x[1]) for x in best}) > 1:
                        fail(
                            "multiple Flye assembly_info tables tie for contig overlap: "
                            + ", ".join(x[1] for x in best)
                        )
                    for contig, value in best[0][2].items():
                        if contig in name_set:
                            calls[contig] = value
    if not source_of:
        fail("mapping manifest/assembler FASTAs contain no contigs")
    return source_of, calls, fasta_of


def _resolve_atomic_tsv(path: str, header: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def _resolve_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--mapping-manifest", required=True)
    parser.add_argument("--mag-circularity", required=True)
    parser.add_argument("--flye-assembly-info", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mapping, counts = load_mag_map(args.mag_contig_map)
    contig_to_mag = {contig: mag for mag, contig in mapping}
    mag_ids = set(counts)
    selected, mag_calls = resolve_call_table(args.mag_circularity, mag_ids)
    source_of, source_calls, _fasta_of = load_source_calls(
        args.mapping_manifest, args.flye_assembly_info
    )

    missing_mag_contigs = sorted(set(contig_to_mag) - set(source_of))
    if missing_mag_contigs:
        fail(f"MAG contigs absent from all assemblies: {missing_mag_contigs[:10]}")

    rows: List[List[object]] = []
    single_conflicts: List[str] = []
    for contig in sorted(source_of):
        assembler = source_of[contig]
        mag = contig_to_mag.get(contig, "NA")
        mag_count = counts.get(mag, 0) if mag != "NA" else 0
        mag_call = mag_calls.get(mag, "NA") if mag != "NA" else "NA"
        direct = source_calls.get(contig)
        if direct is not None:
            call, method = direct
        elif mag != "NA" and mag_count == 1 and mag_call in {"Y", "N"}:
            call, method = mag_call, "SINGLE_CONTIG_MAG_LEVEL_FALLBACK"
        else:
            call, method = "NA", "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED"
        conflict = "N"
        if mag != "NA" and mag_count == 1 and call in {"Y", "N"} and mag_call in {"Y", "N"} and call != mag_call:
            conflict = "Y"
            single_conflicts.append(f"{mag}:{contig}:{call}!={mag_call}")
        rows.append([
            mag, contig, mag_count, "Y" if mag != "NA" else "N", assembler, call, method,
            mag_call, conflict, selected if mag != "NA" else "NA",
        ])


    _resolve_atomic_tsv(
        args.output,
        [
            "MAG_ID", "Contig", "MAG_contig_count", "IN_FINAL_MAG_Y_N", "Source_assembler",
            "SOURCE_CONTIG_ASSEMBLER_CIRCULARITY", "Circularity_resolution_method",
            "MAG_level_assembler_circularity", "Single_contig_source_vs_MAG_conflict_Y_N",
            "MAG_level_circularity_source_table",
        ],
        rows,
    )
    if single_conflicts:
        fail(
            "contig circularity conflicts with MAG-level circularity for single-contig MAGs: "
            + ", ".join(single_conflicts[:10])
        )
    resolved = sum(1 for row in rows if row[5] in {"Y", "N"})
    circular = sum(1 for row in rows if row[5] == "Y")
    print(f"{VERSION}: contigs={len(rows)} resolved={resolved} circular={circular} unresolved={len(rows)-resolved}")
    print(f"Output: {args.output}")


def cmd_resolve(argv):
    """Resolve which contigs each assembler reported as circular"""
    sys.argv = [sys.argv[0] + " resolve", *argv]
    return _resolve_main()


_detect_VERSION = "detect_false_closure"


@dataclass(frozen=True)
class MagInfo:
    mag_id: str
    contigs: Tuple[str, ...]
    genome_length: int

    @property
    def contig_count(self) -> int:
        return len(self.contigs)

    @property
    def contig(self) -> str:
        return self.contigs[0] if len(self.contigs) == 1 else "NA"

    @property
    def eligibility(self) -> str:
        return "Y" if len(self.contigs) == 1 else "N"


@dataclass(frozen=True)
class SourceMapping:
    assembler: str
    assembly_fasta: str
    hifi_all_bam: str
    hifi_primary_supplementary_bam: str
    lengths: Dict[str, int]


@dataclass(frozen=True)
class Breakpoint:
    mag_id: str
    source_assembler: str
    contig: str
    contig_length: int
    position_1based: int
    region_type: str
    distance_to_nearest_end_bp: int
    coverage: int
    clipping_count: int
    clipping_fraction: float
    n_all_clipped_qnames_from_validation: int
    n_exact_spanning_qnames_from_validation: int


@dataclass(frozen=True)
class LargeInsertionSignal:
    mag_id: str
    source_assembler: str
    contig: str
    position_1based: int
    region_type: str
    n_spanning_insertion_qnames: int
    insertion_length_bp_min: int
    insertion_length_bp_median: str
    insertion_length_bp_max: int


def yes(value: str) -> bool:
    return (value or "").strip().upper() in {"Y", "YES", "TRUE", "1"}


def load_resolved_contig_circularity(path: str) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Load assembly-wide contig circularity calls."""
    calls: Dict[str, str] = {}
    methods: Dict[str, str] = {}
    mag_calls: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        required = {
            "MAG_ID", "Contig", "SOURCE_CONTIG_ASSEMBLER_CIRCULARITY",
            "Circularity_resolution_method", "MAG_level_assembler_circularity",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"resolved contig circularity table lacks columns {sorted(missing)}: {path}")
        for line_no, row in enumerate(reader, start=2):
            mag = (row.get("MAG_ID") or "NA").strip() or "NA"
            contig = (row.get("Contig") or "").strip()
            call = (row.get("SOURCE_CONTIG_ASSEMBLER_CIRCULARITY") or "NA").strip().upper()
            method = (row.get("Circularity_resolution_method") or "").strip()
            mag_call = (row.get("MAG_level_assembler_circularity") or "NA").strip().upper()
            if not contig or contig in calls:
                fail(f"blank/duplicate contig in resolved circularity table at line {line_no}: {contig}")
            if call not in {"Y", "N", "NA"} or mag_call not in {"Y", "N", "NA"}:
                fail(f"invalid circularity call at {path}:{line_no}: {call}/{mag_call}")
            calls[contig] = call
            methods[contig] = method or "NA"
            if mag != "NA":
                previous = mag_calls.get(mag)
                if previous is not None and previous != mag_call:
                    fail(f"MAG-level circularity differs among contig rows for {mag}")
                mag_calls[mag] = mag_call
    if not calls:
        fail(f"resolved circularity table contains no contigs: {path}")
    return calls, methods, mag_calls


def load_fai_lengths(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                fail(f"malformed FAI row {line_number}: {path}")
            name = fields[0]
            try:
                length = int(fields[1])
            except ValueError as exc:
                fail(f"invalid FAI length at {path}:{line_number}: {exc}")
            if not name or length < 1 or name in lengths:
                fail(f"blank/duplicate contig in FAI: {name}")
            lengths[name] = length
    if not lengths:
        fail(f"FAI contains no sequences: {path}")
    return lengths


def load_mag_info(map_path: str, assemblies_fai: str) -> List[MagInfo]:
    lengths = load_fai_lengths(assemblies_fai)
    grouped: Dict[str, List[str]] = defaultdict(list)
    seen_contigs: Dict[str, str] = {}
    with open(map_path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        lower = {name.lower(): name for name in fields}
        mag_col = lower.get("mag_id") or lower.get("id")
        contig_col = lower.get("contig") or lower.get("original_contig")
        if not mag_col or not contig_col:
            fail("MAG contig map requires MAG_ID/ID and Contig columns")
        for line_number, row in enumerate(reader, start=2):
            mag_id = (row.get(mag_col) or "").strip()
            contig = (row.get(contig_col) or "").strip()
            if not mag_id or not contig:
                fail(f"blank MAG/contig in MAG contig map line {line_number}")
            if contig not in lengths:
                fail(f"MAG contig absent from all-assemblies FAI: {contig}")
            previous = seen_contigs.get(contig)
            if previous is not None and previous != mag_id:
                fail(f"contig assigned to multiple MAGs: {contig}")
            if previous is None:
                grouped[mag_id].append(contig)
                seen_contigs[contig] = mag_id
    if not grouped:
        fail("MAG contig map contains no MAGs")
    return [
        MagInfo(
            mag_id=mag_id,
            contigs=tuple(contigs),
            genome_length=sum(lengths[contig] for contig in contigs),
        )
        for mag_id, contigs in sorted(grouped.items())
    ]


def load_mapping_manifest(path: str) -> Dict[str, SourceMapping]:
    mappings: Dict[str, SourceMapping] = {}
    seen_contigs: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Assembler", "Assembly_FASTA", "HiFi_all_BAM",
            "HiFi_primary_supplementary_BAM",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"mapping manifest lacks columns {sorted(missing)}: {path}")
        for line_number, row in enumerate(reader, start=2):
            assembler = (row.get("Assembler") or "").strip()
            assembly = (row.get("Assembly_FASTA") or "").strip()
            bam = (row.get("HiFi_all_BAM") or "").strip()
            ps_bam = (row.get("HiFi_primary_supplementary_BAM") or "").strip()
            if not assembler or not assembly or not bam or not ps_bam or assembler in mappings:
                fail(f"invalid/duplicate mapping manifest row {line_number}")
            fai = assembly + ".fai"
            for required_path, label in (
                (assembly, "assembly"),
                (fai, "FAI"),
                (bam, "HiFi all-mapped BAM"),
                (ps_bam, "HiFi primary+supplementary BAM"),
            ):
                if not os.path.isfile(required_path) or os.path.getsize(required_path) == 0:
                    fail(f"{label} missing for {assembler}: {required_path}")
            lengths = load_fai_lengths(fai)
            for contig in lengths:
                previous = seen_contigs.get(contig)
                if previous is not None:
                    fail(
                        "contig ID occurs in multiple assemblers: "
                        f"{contig} ({previous}, {assembler})"
                    )
                seen_contigs[contig] = assembler
            mappings[assembler] = SourceMapping(assembler, assembly, bam, ps_bam, lengths)
    if not mappings:
        fail("mapping manifest contains no assembler rows")
    return mappings


def load_breakpoints(
    path: str,
    mappings: Dict[str, SourceMapping],
    contig_to_mag: Dict[str, str],
) -> Tuple[List[Breakpoint], Set[str]]:
    """Load assembly-wide final HiFi Clipping=Y sites."""
    source_of = {contig: assembler for assembler, src in mappings.items() for contig in src.lengths}
    out: List[Breakpoint] = []
    unresolved_contigs: Set[str] = set()
    seen: Set[Tuple[str, int]] = set()
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        fields = set(reader.fieldnames or [])
        required = {
            "Assembler", "contig", "length", "pos_1based",
            "region_type", "distance_to_nearest_end_bp", "reference_boundary",
            "coverage", "clipping_count", "clipping_ratio",
            "Clipping_analyzable", "Clipping",
            "validated_clipping", "n_unique_all_clipped_qnames",
        }
        missing = required - fields
        if missing:
            fail(f"all-region clipping table lacks columns {sorted(missing)}: {path}")
        spanning_column = (
            "n_unique_exact_spanning_qnames"
            if "n_unique_exact_spanning_qnames" in fields
            else ("n_unique_clean_spanning_qnames" if "n_unique_clean_spanning_qnames" in fields else None)
        )
        if spanning_column is None:
            fail(
                "all-region clipping table lacks a spanning-QNAME column; "
                f"expected n_unique_exact_spanning_qnames: {path}"
            )
        for line_number, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            assembler = (row.get("Assembler") or "").strip()
            if contig not in source_of:
                fail(f"clipping contig absent from assemblies at line {line_number}: {contig}")
            if source_of[contig] != assembler:
                fail(f"assembler mismatch for clipping contig {contig}: {assembler} vs {source_of[contig]}")
            clipping_analyzable = (row.get("Clipping_analyzable") or "").strip().upper()
            clipping_call = (row.get("Clipping") or "NA").strip().upper()
            if clipping_analyzable not in {"Y", "N"}:
                fail(f"invalid Clipping_analyzable at line {line_number}: {clipping_analyzable}")
            if clipping_analyzable == "N":
                if clipping_call != "NA":
                    fail(f"Clipping_analyzable=N requires Clipping=NA at line {line_number}")
                unresolved_contigs.add(contig)
                continue
            if clipping_call not in {"Y", "N"}:
                fail(f"Clipping_analyzable=Y requires Clipping=Y/N at line {line_number}")
            if clipping_call != "Y":
                continue
            if not yes(row.get("validated_clipping", "")):
                fail(
                    "Clipping=Y is inconsistent with validated_clipping=N at "
                    f"{contig}:{row.get('pos_1based', '?')}"
                )
            if yes(row.get("reference_boundary", "")):
                fail(
                    "validated_clipping=Y is inconsistent with reference_boundary=Y at "
                    f"{contig}:{row.get('pos_1based', '?')}"
                )
            source = mappings[assembler]
            region = (row.get("region_type") or "").strip().upper()
            if region not in {"INTERNAL", "TERMINAL"}:
                fail(f"invalid region_type at line {line_number}: {region}")
            try:
                length = int(row["length"]); position = int(row["pos_1based"]); distance = int(row["distance_to_nearest_end_bp"])
                coverage = int(row["coverage"]); clipping_count = int(row["clipping_count"]); clipping_fraction = float(row["clipping_ratio"])
                n_clipped = int(row["n_unique_all_clipped_qnames"]); n_spanning = int(row[spanning_column])
            except (TypeError, ValueError) as exc:
                fail(f"invalid clipping row at line {line_number}: {exc}")
            if length != source.lengths[contig]:
                fail(f"clipping length disagrees with the assembler FAI for {contig}")
            if not 1 <= position <= length or distance < 0:
                fail(f"invalid clipping coordinate for {contig}: {position}")
            if n_clipped < 2 or n_spanning != 0:
                fail(
                    "Clipping=Y must mean >=2 clipped QNAMEs and zero exact "
                    f"spanning QNAMEs: {contig}:{position}"
                )
            key = (contig, position)
            if key in seen:
                fail(f"duplicate validated clipping breakpoint: {contig}:{position}")
            seen.add(key)
            out.append(Breakpoint(
                mag_id=contig_to_mag.get(contig, "NA"), source_assembler=assembler, contig=contig,
                contig_length=length, position_1based=position, region_type=region,
                distance_to_nearest_end_bp=distance, coverage=coverage,
                clipping_count=clipping_count, clipping_fraction=clipping_fraction,
                n_all_clipped_qnames_from_validation=n_clipped,
                n_exact_spanning_qnames_from_validation=n_spanning,
            ))
    return (
        sorted(out, key=lambda b: (b.source_assembler, b.contig, b.position_1based)),
        unresolved_contigs,
    )


def load_large_insertion_signals(
    path: str,
    contig_to_mag: Dict[str, str],
    min_insertion_bp: int,
    min_qnames: int,
) -> List[LargeInsertionSignal]:
    """Load the assembly-wide large read-side insertion sensitivity screen."""
    out: List[LargeInsertionSignal] = []
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        fields = set(reader.fieldnames or [])
        required = {
            "Assembler", "contig", "pos_1based", "region_type",
            "Clipping_analyzable", "n_unique_spanning_with_I_qnames",
            "spanning_I_length_bp_min", "spanning_I_length_bp_median",
            "spanning_I_length_bp_max",
        }
        missing = required - fields
        if missing:
            fail(
                "all-region clipping table lacks large-insertion audit columns "
                f"{sorted(missing)}: {path}"
            )
        for line_number, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            if (row.get("Clipping_analyzable") or "").strip().upper() != "Y":
                continue
            try:
                n_qnames = int(row.get("n_unique_spanning_with_I_qnames") or 0)
            except ValueError as exc:
                fail(f"invalid insertion-QNAME count at line {line_number}: {exc}")
            if n_qnames < min_qnames:
                continue
            min_text = (row.get("spanning_I_length_bp_min") or "NA").strip()
            median_text = (row.get("spanning_I_length_bp_median") or "NA").strip()
            max_text = (row.get("spanning_I_length_bp_max") or "NA").strip()
            if min_text.upper() == "NA" or max_text.upper() == "NA":
                fail(f"insertion-bearing spanning QNAMEs lack insertion-length audit at line {line_number}")
            try:
                min_len = int(float(min_text)); max_len = int(float(max_text)); position = int(row["pos_1based"])
            except (TypeError, ValueError) as exc:
                fail(f"invalid large-insertion audit row at line {line_number}: {exc}")
            if min_len < min_insertion_bp:
                continue
            region = (row.get("region_type") or "").strip().upper()
            if region not in {"INTERNAL", "TERMINAL"}:
                fail(f"invalid region_type at line {line_number}: {region}")
            out.append(LargeInsertionSignal(
                mag_id=contig_to_mag.get(contig, "NA"),
                source_assembler=(row.get("Assembler") or "").strip(),
                contig=contig, position_1based=position, region_type=region,
                n_spanning_insertion_qnames=n_qnames, insertion_length_bp_min=min_len,
                insertion_length_bp_median=median_text, insertion_length_bp_max=max_len,
            ))
    return sorted(out, key=lambda x: (x.source_assembler, x.contig, x.position_1based))


def load_breakpoint_sources(
    bam: str, breakpoints: Sequence[Breakpoint]
) -> Dict[Tuple[str, int], Dict[str, Dict[str, Alignment]]]:
    """Recover exact-position clipped QNAMEs from one source-assembler BAM."""
    if not breakpoints:
        return {}
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as bed:
        bed_path = bed.name
        for bp in breakpoints:
            pos0 = bp.position_1based - 1
            bed.write(f"{bp.contig}\t{pos0}\t{pos0 + 1}\n")

    by_contig: Dict[str, Dict[str, List[Alignment]]] = defaultdict(
        lambda: defaultdict(list)
    )
    try:
        for line in stream_sam(
            ["samtools", "view", "-F", "4", "-q", "0", "-L", bed_path, bam]
        ):
            aln = parse_alignment(line)
            if aln is not None:
                by_contig[aln.rname][aln.qname].append(aln)
    finally:
        os.unlink(bed_path)

    result: Dict[Tuple[str, int], Dict[str, Dict[str, Alignment]]] = {}
    for bp in breakpoints:
        pos0 = bp.position_1based - 1
        sides: Dict[str, Dict[str, Alignment]] = {"left": {}, "right": {}}
        for qname, alignments in by_contig.get(bp.contig, {}).items():
            left_candidates = [
                aln
                for aln in alignments
                if aln.right_clip > 0 and (aln.rend - 1) == pos0
            ]
            right_candidates = [
                aln
                for aln in alignments
                if aln.left_clip > 0 and aln.rstart == pos0
            ]
            rep = representative(left_candidates)
            if rep is not None:
                sides["left"][qname] = rep
            rep = representative(right_candidates)
            if rep is not None:
                sides["right"][qname] = rep
        recovered = set(sides["left"]) | set(sides["right"])
        if len(recovered) != bp.n_all_clipped_qnames_from_validation:
            fail(
                "exact clipped-QNAME recovery disagrees with 03_clipping.sh validation at "
                f"{bp.contig}:{bp.position_1based}; "
                f"validation={bp.n_all_clipped_qnames_from_validation}, recovered={len(recovered)}"
            )
        result[(bp.contig, bp.position_1based)] = sides
    return result


def false_closure_connection_rule(
    left_connected_qnames: int,
    right_connected_qnames: int,
    min_connected_qnames_per_side: int,
) -> bool:
    return (
        left_connected_qnames >= min_connected_qnames_per_side
        and right_connected_qnames >= min_connected_qnames_per_side
    )


def region_call(records: Iterable[Tuple[Breakpoint, bool, Dict[str, Set[str]], Dict[str, Set[str]]]], region: str) -> bool:
    return any(record[1] for record in records if record[0].region_type == region)


def run(args: argparse.Namespace) -> None:
    mag_info = load_mag_info(args.mag_contig_map, args.assemblies_fai)
    mag_by_id = {row.mag_id: row for row in mag_info}
    combined_lengths = load_fai_lengths(args.assemblies_fai)
    ribosomal_families = load_ribosomal_families(args.scg_tsv)
    contig_to_mag = load_contig_to_mag(args.mag_contig_map)
    mappings = load_mapping_manifest(args.mapping_manifest)
    contig_circularity, contig_circularity_method, assembler_calls = load_resolved_contig_circularity(
        args.contig_circularity
    )
    source_contigs = {contig for source in mappings.values() for contig in source.lengths}
    if set(contig_circularity) != source_contigs:
        missing = sorted(source_contigs - set(contig_circularity))[:10]
        extra = sorted(set(contig_circularity) - source_contigs)[:10]
        fail(f"resolved circularity/assembly contig mismatch: missing={missing}, extra={extra}")
    breakpoints, unresolved_clipping_contigs = load_breakpoints(
        args.clipping_all_final_tsv, mappings, contig_to_mag
    )
    large_insertion_signals = load_large_insertion_signals(
        args.clipping_all_final_tsv, contig_to_mag,
        args.large_insertion_min_bp, args.large_insertion_min_qnames,
    )
    large_insertion_by_mag: Dict[str, List[LargeInsertionSignal]] = defaultdict(list)
    for signal in large_insertion_signals:
        large_insertion_by_mag[signal.mag_id].append(signal)

    if args.large_insertion_details:
        sensitivity_rows: List[List[object]] = []
        for signal in large_insertion_signals:
            mag = mag_by_id.get(signal.mag_id)
            mag_call = assembler_calls.get(signal.mag_id, "NA") if signal.mag_id != "NA" else "NA"
            contig_call = contig_circularity.get(signal.contig, "NA")
            sensitivity_analyzable = contig_call == "Y"
            sensitivity_rows.append([
                signal.mag_id, signal.source_assembler, signal.contig,
                signal.position_1based, signal.region_type,
                signal.n_spanning_insertion_qnames,
                signal.insertion_length_bp_min, signal.insertion_length_bp_median,
                signal.insertion_length_bp_max, mag.contig_count if mag is not None else 0, contig_call, mag_call,
                contig_circularity_method.get(signal.contig, "NA"),
                "Y" if sensitivity_analyzable else "N",
                "Y" if sensitivity_analyzable else "NA",
                ("SOURCE_CIRCULAR_CONTIG_LARGE_INSERTION_SENSITIVITY"
                 if sensitivity_analyzable
                 else ("SOURCE_CONTIG_CIRCULARITY_UNRESOLVED" if contig_call == "NA"
                       else "SOURCE_CONTIG_NOT_REPORTED_CIRCULAR")),
            ])
        atomic_tsv(
            args.large_insertion_details,
            [
                "ID", "Source_assembler", "Contig", "Breakpoint_1based",
                "Region_type", "Spanning_insertion_QNAMEs",
                "Insertion_length_bp_min", "Insertion_length_bp_median",
                "Insertion_length_bp_max", "Contig_count",
                "SOURCE_CONTIG_ASSEMBLER_CIRCULARITY", "MAG_level_assembler_circularity",
                "Circularity_resolution_method",
                "LARGE_INSERTION_SENSITIVITY_ANALYZABLE_Y_N",
                "LARGE_INSERTION_SENSITIVITY_SUPPORT_Y_N_NA", "Sensitivity_reason",
            ],
            sensitivity_rows,
        )

    sources: Dict[Tuple[str, int], Dict[str, Dict[str, Alignment]]] = {}
    primary_supplementary_by_source_read: Dict[str, Dict[str, List[Alignment]]] = {}
    by_assembler: Dict[str, List[Breakpoint]] = defaultdict(list)
    for bp in breakpoints:
        by_assembler[bp.source_assembler].append(bp)

    for assembler, assembler_breakpoints in sorted(by_assembler.items()):
        bam = mappings[assembler].hifi_all_bam
        source_rows = load_breakpoint_sources(bam, assembler_breakpoints)
        sources.update(source_rows)
        candidate_qnames: Set[str] = set()
        for sides in source_rows.values():
            candidate_qnames.update(sides.get("left", {}))
            candidate_qnames.update(sides.get("right", {}))
        primary_supplementary_by_source_read[assembler] = load_alignments_for_qnames(
            mappings[assembler].hifi_primary_supplementary_bam, candidate_qnames
        )


    if args.candidate_read_output:
        candidate_rows: List[List[object]] = []
        for bp in breakpoints:
            sides = sources.get((bp.contig, bp.position_1based), {})
            for side in ("left", "right"):
                for qname in sorted(sides.get(side, {})):
                    candidate_rows.append(
                        [
                            bp.mag_id,
                            bp.source_assembler,
                            bp.contig,
                            bp.position_1based,
                            bp.region_type,
                            side.upper(),
                            qname,
                        ]
                    )
        atomic_tsv(
            args.candidate_read_output,
            [
                "ID", "Source_assembler", "Contig", "Breakpoint_1based",
                "Region_type", "Breakpoint_side", "Read_ID",
            ],
            candidate_rows,
        )

    detail_header = [
        "ID", "Source_assembler", "Contig", "Breakpoint_1based",
        "Region_type", "Distance_to_nearest_end_bp",
        "Validated_clipping", "Coverage", "Clipping_count", "Clipping_fraction",
        "Validated_all_clipped_QNAMEs", "Validated_clean_spanning_QNAMEs",
        "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_CLIPPED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_CROSS_CONTIG_CONNECTED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_TARGET_CONTIGS",
        "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_CLIPPED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_CROSS_CONTIG_CONNECTED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_TARGET_CONTIGS",
        "HIFI_PRIMARY_SUPPLEMENTARY_BILATERAL_CROSS_CONTIG_CONNECTION_SUPPORTED_Y_N",
        "Left_clipped_QNAMEs", "Left_cross_contig_connected_QNAMEs",
        "Left_target_contigs", "Right_clipped_QNAMEs",
        "Right_cross_contig_connected_QNAMEs", "Right_target_contigs",
        "Bilateral_cross_contig_connection_evidence", "False_closure_evidence",
        "Assembler_reported_circularity", "MAG_level_assembler_reported_circularity",
        "Circularity_resolution_method", "False_closure_analyzable",
        "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_Y_N",
        "False_closure", "False_closure_reason",
    ]

    detail_rows: List[List[object]] = []
    records: List[Tuple] = []
    for bp in breakpoints:
        sides = sources.get(
            (bp.contig, bp.position_1based), {"left": {}, "right": {}}
        )
        left = sides.get("left", {})
        right = sides.get("right", {})
        ps_left = {q: a for q, a in left.items() if not a.secondary}
        ps_right = {q: a for q, a in right.items() if not a.secondary}
        primary_supplementary_by_read = primary_supplementary_by_source_read.get(
            bp.source_assembler, {}
        )

        left_connected, left_targets = summarize_external_connection(
            ps_left, primary_supplementary_by_read, bp.contig
        )
        right_connected, right_targets = summarize_external_connection(
            ps_right, primary_supplementary_by_read, bp.contig
        )
        evidence = false_closure_connection_rule(
            len(left_connected),
            len(right_connected),
            args.min_linked_reads_per_side,
        )

        mag_call = assembler_calls.get(bp.mag_id, "NA") if bp.mag_id != "NA" else "NA"
        contig_call = contig_circularity.get(bp.contig, "NA")
        circularity_method = contig_circularity_method.get(bp.contig, "NA")
        false_closure_analyzable = (
            contig_call == "Y"
            and bp.contig not in unresolved_clipping_contigs
        )
        false_closure = (
            "Y" if false_closure_analyzable and evidence
            else ("N" if false_closure_analyzable else "NA")
        )

        if contig_call == "NA":
            reason = "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED"
        elif contig_call != "Y":
            reason = "SOURCE_CONTIG_NOT_REPORTED_CIRCULAR"
        elif bp.contig in unresolved_clipping_contigs:
            reason = "UNRESOLVED_HIFI_CLIPPING_CANDIDATE_ON_CIRCULAR_CONTIG"
        elif len(left_connected) < args.min_linked_reads_per_side:
            reason = "INSUFFICIENT_LEFT_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_CONNECTED_READS"
        elif len(right_connected) < args.min_linked_reads_per_side:
            reason = "INSUFFICIENT_RIGHT_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_CONNECTED_READS"
        else:
            reason = "VALIDATED_CLIPPING_WITH_BILATERAL_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_SPLIT_CONNECTIONS"

        detail_rows.append([
            bp.mag_id,
            bp.source_assembler,
            bp.contig,
            bp.position_1based,
            bp.region_type,
            bp.distance_to_nearest_end_bp,
            "Y",
            bp.coverage,
            bp.clipping_count,
            f"{bp.clipping_fraction:.6f}",
            bp.n_all_clipped_qnames_from_validation,
            bp.n_exact_spanning_qnames_from_validation,
            len(ps_left), len(left_connected), target_count_text(left_targets),
            len(ps_right), len(right_connected), target_count_text(right_targets),
            "Y" if evidence else "N",
            len(ps_left), len(left_connected), target_count_text(left_targets),
            len(ps_right), len(right_connected), target_count_text(right_targets),
            "Y" if evidence else "N", "Y" if evidence else "N",
            contig_call, mag_call, circularity_method,
            "Y" if false_closure_analyzable else "N",
            false_closure, false_closure, reason,
        ])
        records.append((bp, evidence, left_targets, right_targets))

    atomic_tsv(args.details, detail_header, detail_rows)

    by_mag: Dict[str, List[Tuple]] = defaultdict(list)
    for record in records:
        by_mag[record[0].mag_id].append(record)

    summary_header = [
        "ID", "Contig_count", "Contig", "Genome_length_bp",
        "Assembler_reported_circularity",
        "MAG_level_input_assembler_circularity",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "False_closure_analyzable",
        "Paper_ribosomal_protein_family_count",
        "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N",
        "Paper_small_genome_ribosomal_criterion_Y_N",
        "Paper_small_genome_proxy", "Paper_small_genome_proxy_contigs",
        "Premature_circularization_candidate",
        "Premature_circularization_screen",
        "Large_insertion_sensitivity_analyzable",
        "Large_insertion_sensitivity_support",
        "Large_insertion_sensitivity_site_count",
        "Large_insertion_sensitivity_sites",
        "Large_insertion_sensitivity_min_bp",
        "Large_insertion_sensitivity_min_QNAMEs",
        "All_MAG_validated_clipping_breakpoint_count",
        "Validated_clipping_breakpoint_count",
        "Validated_internal_clipping_breakpoint_count",
        "Validated_terminal_clipping_breakpoint_count",
        "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_BREAKPOINT_COUNT",
        "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_Y_N",
        "False_closure_breakpoint_count",
        "Internal_false_closure_breakpoint_count",
        "Terminal_false_closure_breakpoint_count",
        "False_closure_breakpoints", "False_closure_regions",
        "Candidate_additional_target_contigs", "Candidate_additional_target_MAGs",
        "Internal_false_closure_evidence", "Terminal_false_closure_evidence",
        "False_closure_evidence", "False_closure",
        "Min_connected_QNAMEs_per_side", "False_closure_reason",
    ]
    summary_rows: List[List[object]] = []
    for mag in mag_info:
        recs = by_mag.get(mag.mag_id, [])
        circular_contigs = [c for c in mag.contigs if contig_circularity.get(c, "NA") == "Y"]
        unresolved_circularity_contigs = [c for c in mag.contigs if contig_circularity.get(c, "NA") == "NA"]
        resolved_circularity_count = mag.contig_count - len(unresolved_circularity_contigs)
        circularity_complete = not unresolved_circularity_contigs
        circular_recs = [r for r in recs if contig_circularity.get(r[0].contig, "NA") == "Y"]
        positives = [r for r in circular_recs if r[1]]
        internal = [r for r in circular_recs if r[0].region_type == "INTERNAL"]
        terminal = [r for r in circular_recs if r[0].region_type == "TERMINAL"]
        internal_pos = [r for r in positives if r[0].region_type == "INTERNAL"]
        terminal_pos = [r for r in positives if r[0].region_type == "TERMINAL"]

        circular_clipping_unresolved = any(c in unresolved_clipping_contigs for c in circular_contigs)
        definitive_negative_possible = (
            bool(circular_contigs)
            and circularity_complete
            and not circular_clipping_unresolved
        )
        evidence = bool(positives)
        false_closure = "Y" if evidence else ("N" if definitive_negative_possible else "NA")
        analyzable = false_closure in {"Y", "N"}

        positions = (
            ";".join(
                f"{r[0].contig}|{r[0].region_type}:{r[0].position_1based}"
                for r in positives
            ) or "NA"
        )
        regions = ";".join(sorted({r[0].region_type for r in positives})) or "NA"
        targets: Set[str] = set()
        for r in positives:
            targets.update(r[2].keys())
            targets.update(r[3].keys())
        target_contigs = ";".join(sorted(targets)) if targets else "NA"
        target_mags = (
            ";".join(sorted({contig_to_mag[t] for t in targets if t in contig_to_mag})) or "NA"
        )

        large_sites = [
            site for site in large_insertion_by_mag.get(mag.mag_id, [])
            if contig_circularity.get(site.contig, "NA") == "Y"
        ]
        if large_sites:
            large_insertion_support = "Y"
        elif definitive_negative_possible:
            large_insertion_support = "N"
        else:
            large_insertion_support = "NA"
        large_insertion_analyzable = large_insertion_support in {"Y", "N"}
        large_insertion_sites_text = (
            ";".join(
                f"{site.contig}|{site.region_type}:{site.position_1based}"
                f"[n={site.n_spanning_insertion_qnames},minI={site.insertion_length_bp_min}]"
                for site in large_sites
            ) if large_sites else "NA"
        )


        circular_family_map = {c: ribosomal_families.get(c, ()) for c in circular_contigs}
        paper_family_count = max((len(v) for v in circular_family_map.values()), default=0)
        max_family_contigs = [c for c, v in circular_family_map.items() if len(v) == paper_family_count and paper_family_count > 0]
        paper_families = tuple(sorted({f for c in max_family_contigs for f in circular_family_map[c]}))
        paper_length_ok = any(combined_lengths[c] < args.premature_max_length_bp for c in circular_contigs)
        paper_ribosomal_ok = any(
            len(circular_family_map[c]) >= args.premature_min_ribosomal_families for c in circular_contigs
        )
        paper_proxy_contigs = [
            c for c in circular_contigs
            if combined_lengths[c] < args.premature_max_length_bp
            and len(circular_family_map[c]) >= args.premature_min_ribosomal_families
        ]
        paper_proxy = bool(paper_proxy_contigs)
        source_aggregate_call = (
            "Y" if circular_contigs else ("N" if circularity_complete else "NA")
        )
        mag_level_input_call = assembler_calls.get(mag.mag_id, "NA")
        premature_candidate = mag.contig_count == 1 and source_aggregate_call == "Y" and evidence
        premature_screen = premature_candidate or paper_proxy

        if evidence:
            reason = "SOURCE_CIRCULAR_CONTIG_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_STRUCTURE_SUPPORTED"
        elif not circular_contigs and unresolved_circularity_contigs:
            reason = "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED"
        elif not circular_contigs:
            reason = "NOT_APPLICABLE_NO_SOURCE_CIRCULAR_CONTIG"
        elif not circularity_complete:
            reason = "SOURCE_CONTIG_CIRCULARITY_INCOMPLETE_FOR_DEFINITIVE_NEGATIVE"
        elif circular_clipping_unresolved:
            reason = "UNRESOLVED_HIFI_CLIPPING_CANDIDATE_ON_SOURCE_CIRCULAR_CONTIG"
        elif not circular_recs:
            reason = "ANALYZABLE_SOURCE_CIRCULAR_CONTIG_WITH_NO_HIFI_CLIPPING_Y_BREAKPOINT"
        else:
            reason = "ANALYZABLE_CIRCULAR_CONTIG_CLIPPING_LACKED_BILATERAL_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_SPLIT_CONNECTIONS"

        summary_rows.append([
            mag.mag_id, mag.contig_count, mag.contig, mag.genome_length,
            source_aggregate_call, mag_level_input_call,
            len(circular_contigs), ";".join(circular_contigs) if circular_contigs else "NA",
            resolved_circularity_count, len(unresolved_circularity_contigs),
            "Y" if circularity_complete else "N",
            "Y" if analyzable else "N",
            paper_family_count, ";".join(paper_families) if paper_families else "NA",
            "Y" if paper_length_ok else "N", "Y" if paper_ribosomal_ok else "N",
            "Y" if paper_proxy else "N",
            ";".join(paper_proxy_contigs) if paper_proxy_contigs else "NA",
            "Y" if premature_candidate else "N", "Y" if premature_screen else "N",
            "Y" if large_insertion_analyzable else "N", large_insertion_support,
            len(large_sites), large_insertion_sites_text,
            args.large_insertion_min_bp, args.large_insertion_min_qnames,
            len(recs), len(circular_recs), len(internal), len(terminal),
            len(positives), false_closure,
            len(positives), len(internal_pos), len(terminal_pos),
            positions, regions, target_contigs, target_mags,
            "Y" if internal_pos else "N", "Y" if terminal_pos else "N",
            "Y" if evidence else "N", false_closure,
            args.min_linked_reads_per_side, reason,
        ])

    atomic_tsv(args.output, summary_header, summary_rows)
    print(
        f"{_detect_VERSION}: validated_breakpoints={len(detail_rows)}; "
        f"false_closure_structure_Y={sum(record[1] for record in records)}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies-fai", required=True)
    parser.add_argument("--mapping-manifest", required=True)
    parser.add_argument("--clipping-all-final-tsv", required=True)
    parser.add_argument("--contig-circularity", required=True)
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--scg-tsv", required=True)
    parser.add_argument("--premature-max-length-bp", type=int, default=500000)
    parser.add_argument("--premature-min-ribosomal-families", type=int, default=3)
    parser.add_argument("--min-linked-reads-per-side", type=int, default=2)
    parser.add_argument("--large-insertion-min-bp", type=int, default=100)
    parser.add_argument("--large-insertion-min-qnames", type=int, default=2)
    parser.add_argument("--large-insertion-details", default="")
    parser.add_argument("--candidate-read-output", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--details", required=True)
    return parser


def _detect_main() -> None:
    args = build_parser().parse_args()
    if args.premature_max_length_bp < 1:
        fail("premature-circularization maximum length must be at least 1 bp")
    if args.premature_min_ribosomal_families < 1:
        fail("minimum ribosomal-protein family count must be at least 1")
    if args.min_linked_reads_per_side < 1:
        fail("minimum linked reads per side must be at least 1")
    if args.large_insertion_min_bp < 1:
        fail("large-insertion minimum length must be at least 1 bp")
    if args.large_insertion_min_qnames < 1:
        fail("large-insertion minimum QNAME count must be at least 1")
    run(args)


def cmd_detect(argv):
    """Detect false-closure breakpoints from validated clipping and split HiFi alignments"""
    sys.argv = [sys.argv[0] + " detect", *argv]
    return _detect_main()


_mag_VERSION = "summarize_false_closure"


def read_table(path: str) -> tuple[list[str], Dict[str, Dict[str, str]]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "ID" not in fields:
            fail(f"table lacks ID column: {path}")
        rows: Dict[str, Dict[str, str]] = {}
        for row in reader:
            mag = (row.get("ID") or "").strip()
            if not mag or mag in rows:
                fail(f"blank/duplicate ID in {path}: {mag!r}")
            rows[mag] = row
    return fields, rows


def atomic_dict_table(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{os.path.basename(path)}.", suffix=".tmp",
                               dir=os.path.dirname(os.path.abspath(path)), text=True)
    os.close(fd)
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({f: row.get(f, "NA") for f in fields})
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def call(value: str, allowed: set[str], label: str, mag: str) -> str:
    v = (value or "NA").strip().upper()
    if v not in allowed:
        fail(f"invalid {label} for {mag}: {v!r}")
    return v


def intval(row: Mapping[str, str], field: str, mag: str, minimum: int = 0) -> int:
    try:
        value = int(row.get(field, "0") or 0)
    except ValueError as exc:
        fail(f"invalid integer {field} for {mag}: {exc}")
    if value < minimum:
        fail(f"{field} below {minimum} for {mag}: {value}")
    return value


def _mag_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--false-closure-tsv", required=True)
    parser.add_argument("--false-final", required=True)
    parser.add_argument("--false-audit", required=True)
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()

    source_fields, source = read_table(args.false_closure_tsv)
    required = {
        "ID", "Contig_count", "Contig", "Genome_length_bp",
        "Assembler_reported_circularity", "MAG_level_input_assembler_circularity",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "False_closure_analyzable",
        "Paper_ribosomal_protein_family_count", "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N",
        "Paper_small_genome_ribosomal_criterion_Y_N", "Paper_small_genome_proxy",
        "Paper_small_genome_proxy_contigs", "Premature_circularization_candidate",
        "Premature_circularization_screen",
        "Large_insertion_sensitivity_analyzable", "Large_insertion_sensitivity_support",
        "Large_insertion_sensitivity_site_count", "Large_insertion_sensitivity_sites",
        "Large_insertion_sensitivity_min_bp", "Large_insertion_sensitivity_min_QNAMEs",
        "All_MAG_validated_clipping_breakpoint_count",
        "Validated_clipping_breakpoint_count",
        "Validated_internal_clipping_breakpoint_count",
        "Validated_terminal_clipping_breakpoint_count",
        "False_closure_breakpoint_count", "Internal_false_closure_breakpoint_count",
        "Terminal_false_closure_breakpoint_count", "False_closure_breakpoints",
        "False_closure_regions", "Candidate_additional_target_contigs",
        "Candidate_additional_target_MAGs", "Internal_false_closure_evidence",
        "Terminal_false_closure_evidence", "False_closure_evidence", "False_closure",
        "Min_connected_QNAMEs_per_side", "False_closure_reason",
    }
    missing = required - set(source_fields)
    if missing:
        fail(f"false-closure table lacks columns {sorted(missing)}")

    final_fields = [
        "ID", "MAG_CONTIG_COUNT", "FALSE_CLOSURE_ANALYZABLE", "ASSEMBLER_REPORTED_CIRCULARITY",
        "MAG_LEVEL_INPUT_ASSEMBLER_CIRCULARITY",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "Paper_ribosomal_protein_family_count", "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N", "Paper_small_genome_ribosomal_criterion_Y_N",
        "Paper_small_genome_proxy", "Paper_small_genome_proxy_contigs",
        "Premature_circularization_candidate", "Premature_circularization_screen",
        "LARGE_INSERTION_SENSITIVITY_ANALYZABLE", "LARGE_INSERTION_SENSITIVITY_SUPPORT",
        "LARGE_INSERTION_SENSITIVITY_SITE_COUNT", "LARGE_INSERTION_SENSITIVITY_SITES",
        "Candidate_additional_target_contigs", "Candidate_additional_target_MAGs",
        "VALIDATED_CLIPPING_PRESENT", "INTERNAL_VALIDATED_CLIPPING_PRESENT",
        "TERMINAL_VALIDATED_CLIPPING_PRESENT", "INTERNAL_FALSE_CLOSURE_EVIDENCE",
        "TERMINAL_FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE",
        "FALSE_CLOSURE_BREAKPOINTS", "FALSE_CLOSURE_REGIONS",
    ]
    audit_fields = [
        "ID", "Contig_count", "Contig", "Genome_length_bp",
        "ASSEMBLER_REPORTED_CIRCULARITY", "MAG_LEVEL_INPUT_ASSEMBLER_CIRCULARITY",
        "FALSE_CLOSURE_ANALYZABLE",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT", "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "Paper_ribosomal_protein_family_count", "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N", "Paper_small_genome_ribosomal_criterion_Y_N",
        "Paper_small_genome_proxy", "Paper_small_genome_proxy_contigs",
        "Premature_circularization_candidate", "Premature_circularization_screen",
        "LARGE_INSERTION_SENSITIVITY_ANALYZABLE", "LARGE_INSERTION_SENSITIVITY_SUPPORT",
        "LARGE_INSERTION_SENSITIVITY_SITE_COUNT", "LARGE_INSERTION_SENSITIVITY_SITES",
        "Large_insertion_sensitivity_min_bp", "Large_insertion_sensitivity_min_QNAMEs",
        "All_MAG_validated_clipping_breakpoint_count", "Validated_clipping_breakpoint_count",
        "Validated_internal_clipping_breakpoint_count", "Validated_terminal_clipping_breakpoint_count",
        "False_closure_breakpoint_count", "Internal_false_closure_breakpoint_count",
        "Terminal_false_closure_breakpoint_count", "FALSE_CLOSURE_BREAKPOINTS",
        "FALSE_CLOSURE_REGIONS", "Candidate_additional_target_contigs",
        "Candidate_additional_target_MAGs", "INTERNAL_FALSE_CLOSURE_EVIDENCE",
        "TERMINAL_FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE",
        "Min_connected_QNAMEs_per_side", "False_closure_reason",
    ]

    final_rows: List[Dict[str, object]] = []
    audit_rows: List[Dict[str, object]] = []
    for mag, row in source.items():
        contig_count = intval(row, "Contig_count", mag, 1)
        circular_count = intval(row, "SOURCE_CIRCULAR_CONTIG_COUNT", mag)
        resolved_count = intval(row, "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT", mag)
        unresolved_count = intval(row, "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT", mag)
        if resolved_count + unresolved_count != contig_count:
            fail(f"assembler circularity counts do not sum to contig count for {mag}")
        if circular_count > resolved_count:
            fail(f"circular contig count exceeds resolved count for {mag}")

        assembler = call(row["Assembler_reported_circularity"], {"Y", "N", "NA"}, "source-contig aggregate assembler circularity", mag)
        mag_level_input = call(row["MAG_level_input_assembler_circularity"], {"Y", "N", "NA"}, "MAG-level input assembler circularity", mag)
        complete = call(row["SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N"], {"Y", "N"}, "source circularity completeness", mag)
        expected_aggregate = "Y" if circular_count > 0 else ("N" if complete == "Y" else "NA")
        if assembler != expected_aggregate:
            fail(f"contig aggregate circularity inconsistent for {mag}: {assembler} vs {expected_aggregate}")
        if complete != ("Y" if unresolved_count == 0 else "N"):
            fail(f"assembler circularity completeness inconsistent for {mag}")
        analyzable = call(row["False_closure_analyzable"], {"Y", "N"}, "false-closure analyzability", mag)
        false_call = call(row["False_closure"], {"Y", "N", "NA"}, "false-closure call", mag)
        if (analyzable == "Y") != (false_call in {"Y", "N"}):
            fail(f"false-closure analyzability/call mismatch for {mag}: {analyzable}/{false_call}")

        evidence = call(row["False_closure_evidence"], {"Y", "N"}, "false-closure evidence", mag)
        if false_call == "Y" and evidence != "Y":
            fail(f"FALSE_CLOSURE=Y without structural evidence for {mag}")
        premature = call(row["Premature_circularization_candidate"], {"Y", "N"}, "premature candidate", mag)
        proxy = call(row["Paper_small_genome_proxy"], {"Y", "N"}, "paper proxy", mag)
        screen = call(row["Premature_circularization_screen"], {"Y", "N"}, "premature screen", mag)
        expected_premature = "Y" if contig_count == 1 and assembler == "Y" and evidence == "Y" else "N"
        expected_screen = "Y" if expected_premature == "Y" or proxy == "Y" else "N"
        if premature != expected_premature or screen != expected_screen:
            fail(f"premature-circularization logic inconsistent for {mag}")

        validated = intval(row, "Validated_clipping_breakpoint_count", mag)
        internal_validated = intval(row, "Validated_internal_clipping_breakpoint_count", mag)
        terminal_validated = intval(row, "Validated_terminal_clipping_breakpoint_count", mag)
        false_n = intval(row, "False_closure_breakpoint_count", mag)
        internal_false_n = intval(row, "Internal_false_closure_breakpoint_count", mag)
        terminal_false_n = intval(row, "Terminal_false_closure_breakpoint_count", mag)
        if internal_validated + terminal_validated != validated:
            fail(f"validated clipping regional counts inconsistent for {mag}")
        if internal_false_n + terminal_false_n != false_n or false_n > validated:
            fail(f"false-closure breakpoint counts inconsistent for {mag}")
        if (false_n > 0) != (evidence == "Y"):
            fail(f"false-closure evidence/count mismatch for {mag}")

        large_an = call(row["Large_insertion_sensitivity_analyzable"], {"Y", "N"}, "large insertion analyzability", mag)
        large_call = call(row["Large_insertion_sensitivity_support"], {"Y", "N", "NA"}, "large insertion support", mag)
        if (large_an == "Y") != (large_call in {"Y", "N"}):
            fail(f"large-insertion analyzability/call mismatch for {mag}")

        final = {
            "ID": mag,
            "MAG_CONTIG_COUNT": contig_count,
            "FALSE_CLOSURE_ANALYZABLE": analyzable,
            "ASSEMBLER_REPORTED_CIRCULARITY": assembler,
            "MAG_LEVEL_INPUT_ASSEMBLER_CIRCULARITY": mag_level_input,
            "SOURCE_CIRCULAR_CONTIG_COUNT": circular_count,
            "SOURCE_CIRCULAR_CONTIGS": row["SOURCE_CIRCULAR_CONTIGS"],
            "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N": complete,
            "Paper_ribosomal_protein_family_count": row["Paper_ribosomal_protein_family_count"],
            "Paper_ribosomal_protein_families": row["Paper_ribosomal_protein_families"],
            "Paper_small_genome_length_criterion_Y_N": row["Paper_small_genome_length_criterion_Y_N"],
            "Paper_small_genome_ribosomal_criterion_Y_N": row["Paper_small_genome_ribosomal_criterion_Y_N"],
            "Paper_small_genome_proxy": proxy,
            "Paper_small_genome_proxy_contigs": row["Paper_small_genome_proxy_contigs"],
            "Premature_circularization_candidate": premature,
            "Premature_circularization_screen": screen,
            "LARGE_INSERTION_SENSITIVITY_ANALYZABLE": large_an,
            "LARGE_INSERTION_SENSITIVITY_SUPPORT": large_call,
            "LARGE_INSERTION_SENSITIVITY_SITE_COUNT": row["Large_insertion_sensitivity_site_count"],
            "LARGE_INSERTION_SENSITIVITY_SITES": row["Large_insertion_sensitivity_sites"],
            "Candidate_additional_target_contigs": row["Candidate_additional_target_contigs"],
            "Candidate_additional_target_MAGs": row["Candidate_additional_target_MAGs"],
            "VALIDATED_CLIPPING_PRESENT": "Y" if validated else "N",
            "INTERNAL_VALIDATED_CLIPPING_PRESENT": "Y" if internal_validated else "N",
            "TERMINAL_VALIDATED_CLIPPING_PRESENT": "Y" if terminal_validated else "N",
            "INTERNAL_FALSE_CLOSURE_EVIDENCE": row["Internal_false_closure_evidence"],
            "TERMINAL_FALSE_CLOSURE_EVIDENCE": row["Terminal_false_closure_evidence"],
            "FALSE_CLOSURE_EVIDENCE": evidence,
            "FALSE_CLOSURE": false_call,
            "FALSE_CLOSURE_BREAKPOINTS": row["False_closure_breakpoints"],
            "FALSE_CLOSURE_REGIONS": row["False_closure_regions"],
        }
        audit = {
            **final,
            "Contig_count": row["Contig_count"], "Contig": row["Contig"],
            "Genome_length_bp": row["Genome_length_bp"],
            "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT": resolved_count,
            "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT": unresolved_count,
            "Large_insertion_sensitivity_min_bp": row["Large_insertion_sensitivity_min_bp"],
            "Large_insertion_sensitivity_min_QNAMEs": row["Large_insertion_sensitivity_min_QNAMEs"],
            "All_MAG_validated_clipping_breakpoint_count": row["All_MAG_validated_clipping_breakpoint_count"],
            "Validated_clipping_breakpoint_count": validated,
            "Validated_internal_clipping_breakpoint_count": internal_validated,
            "Validated_terminal_clipping_breakpoint_count": terminal_validated,
            "False_closure_breakpoint_count": false_n,
            "Internal_false_closure_breakpoint_count": internal_false_n,
            "Terminal_false_closure_breakpoint_count": terminal_false_n,
            "Min_connected_QNAMEs_per_side": row["Min_connected_QNAMEs_per_side"],
            "False_closure_reason": row["False_closure_reason"],
        }
        final_rows.append(final)
        audit_rows.append(audit)

    atomic_dict_table(args.false_final, final_fields, final_rows)
    atomic_dict_table(args.false_audit, audit_fields, audit_rows)
    with open(args.metadata, "w", encoding="utf-8") as handle:
        handle.write(f"summarizer_version\t{_mag_VERSION}\n")
        handle.write("false_closure_main_reference\tsource_assembler_original_total_assembly\n")
        handle.write("circularity_gate\tsource_contig_assembler_reported_circularity\n")
        handle.write("final_assembler_reported_circularity\tsource_contig_aggregate_Y_if_any_source_contig_Y_N_if_all_resolved_N_else_NA\n")
        handle.write("MAG_level_input_circularity\taudit_and_single_contig_fallback_only\n")
        handle.write("multi_contig_MAG_policy\tFALSE_CLOSURE_Y_if_any_source_circular_contig_supports_false_closure; premature-circularization candidate limited to single-contig MAGs\n")
        handle.write("large_insertion_evidence\tsensitivity_only_no_change_to_false_closure_call\n")
        handle.write("paper_small_genome_proxy_changes_false_closure_call\tN\n")
    print(f"{_mag_VERSION}: MAGs={len(source)}")
    print(f"Final: {args.false_final}")
    print(f"Audit: {args.false_audit}")


def cmd_mag(argv):
    """Write the MAG table and its audit"""
    sys.argv = [sys.argv[0] + " mag", *argv]
    return _mag_main()


_contig_VERSION = "summarize_false_closure_by_contig"


def _contig_fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def pick_field(fieldnames: List[str], candidates: Tuple[str, ...], label: str) -> str:
    lower = {name.strip().lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    _contig_fail(f"missing {label}; observed columns: {fieldnames}")


def _contig_load_mag_map(path: str) -> Tuple[Dict[str, str], Dict[str, int]]:
    contig_to_mag: Dict[str, str] = {}
    counts: Dict[str, int] = defaultdict(int)
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        mag_col = pick_field(fields, ("MAG_ID", "ID"), "MAG ID column")
        contig_col = pick_field(fields, ("Contig", "Original_contig"), "contig column")
        for line_no, row in enumerate(reader, start=2):
            mag = (row.get(mag_col) or "").strip()
            contig = (row.get(contig_col) or "").strip()
            if not mag or not contig:
                _contig_fail(f"blank MAG/contig at {path}:{line_no}")
            previous = contig_to_mag.get(contig)
            if previous is not None and previous != mag:
                _contig_fail(f"contig assigned to multiple MAGs: {contig}")
            if previous is None:
                contig_to_mag[contig] = mag
                counts[mag] += 1
    if not contig_to_mag:
        _contig_fail(f"MAG contig map contains no rows: {path}")
    return contig_to_mag, counts


def load_fai(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                _contig_fail(f"malformed FAI row at {path}:{line_no}")
            contig = fields[0]
            try:
                length = int(fields[1])
            except ValueError as exc:
                _contig_fail(f"invalid FAI length at {path}:{line_no}: {exc}")
            if not contig or length < 1 or contig in lengths:
                _contig_fail(f"blank/duplicate FAI contig: {contig}")
            lengths[contig] = length
    if not lengths:
        _contig_fail(f"FAI contains no contigs: {path}")
    return lengths


def load_circularity(path: str) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "MAG_ID", "Contig", "Source_assembler",
            "SOURCE_CONTIG_ASSEMBLER_CIRCULARITY", "Circularity_resolution_method",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            _contig_fail(f"resolved circularity table lacks {sorted(missing)}")
        for line_no, row in enumerate(reader, start=2):
            contig = (row.get("Contig") or "").strip()
            call = (row.get("SOURCE_CONTIG_ASSEMBLER_CIRCULARITY") or "NA").strip().upper()
            if not contig or contig in out or call not in {"Y", "N", "NA"}:
                _contig_fail(f"invalid/duplicate circularity row at {path}:{line_no}: {contig}/{call}")
            out[contig] = row
    if not out:
        _contig_fail(f"resolved circularity table contains no contigs: {path}")
    return out


def load_clipping_state(path: str) -> Tuple[Dict[str, bool], Dict[str, bool]]:
    """Return validated-clipping presence and unresolved-candidate presence."""
    validated: Dict[str, bool] = defaultdict(bool)
    unresolved: Dict[str, bool] = defaultdict(bool)
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contig", "Clipping_analyzable", "Clipping"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            _contig_fail(f"clipping table lacks columns {sorted(missing)}")
        for line_no, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            analyzable = (row.get("Clipping_analyzable") or "").strip().upper()
            call = (row.get("Clipping") or "NA").strip().upper()
            if not contig or analyzable not in {"Y", "N"}:
                _contig_fail(f"invalid clipping row at {path}:{line_no}")
            if analyzable == "N":
                if call != "NA":
                    _contig_fail(f"Clipping_analyzable=N requires Clipping=NA at {path}:{line_no}")
                unresolved[contig] = True
            elif call == "Y":
                validated[contig] = True
            elif call != "N":
                _contig_fail(f"Clipping_analyzable=Y requires Clipping=Y/N at {path}:{line_no}")
    return validated, unresolved


def load_details(path: str) -> Dict[str, List[dict]]:
    by_contig: Dict[str, List[dict]] = defaultdict(list)
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Contig", "Region_type", "Breakpoint_1based", "False_closure_evidence",
            "Left_cross_contig_connected_QNAMEs", "Right_cross_contig_connected_QNAMEs",
            "Left_target_contigs", "Right_target_contigs",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            _contig_fail(f"false-closure breakpoint table lacks columns {sorted(missing)}")
        for row in reader:
            contig = (row.get("Contig") or "").strip()
            if contig:
                by_contig[contig].append(row)
    return by_contig


def _contig_yes(row: dict, field: str) -> bool:
    return (row.get(field) or "N").strip().upper() == "Y"


def join_unique(values) -> str:
    clean = sorted({value for value in values if value and value != "NA"})
    return ";".join(clean) if clean else "NA"


def target_names(text: str) -> List[str]:
    if not text or text == "NA":
        return []
    return [token.split(":", 1)[0] for token in text.split(";") if token]


def atomic_write(path: str, header: List[str], rows: List[List[object]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for line_no, row in enumerate(rows, start=2):
                if len(row) != len(header):
                    _contig_fail(f"schema mismatch at output line {line_no}: {len(row)} != {len(header)}")
                writer.writerow(row)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def breakpoint_details(rows: List[dict]) -> str:
    """Every validated clipping breakpoint of a contig in one cell (';' between breakpoints)."""
    def field(row: dict, name: str) -> str:
        return (row.get(name) or "NA").strip() or "NA"

    def targets(text: str) -> str:
        return "" if text in ("", "NA") else "[" + text.replace(";", ",") + "]"

    parts = []
    for row in sorted(rows, key=lambda r: int(field(r, "Breakpoint_1based")) if field(r, "Breakpoint_1based").isdigit() else 0):
        parts.append(
            f"{field(row, 'Region_type').upper()}@{field(row, 'Breakpoint_1based')}"
            f"({field(row, 'Distance_to_nearest_end_bp')})"
            f"|L={field(row, 'Left_cross_contig_connected_QNAMEs')}/{field(row, 'Left_clipped_QNAMEs')}"
            f"{targets(field(row, 'Left_target_contigs'))}"
            f"|R={field(row, 'Right_cross_contig_connected_QNAMEs')}/{field(row, 'Right_clipped_QNAMEs')}"
            f"{targets(field(row, 'Right_target_contigs'))}"
            f"|cov={field(row, 'Coverage')},clip={field(row, 'Clipping_count')},frac={field(row, 'Clipping_fraction')}"
            f"|spanning={field(row, 'Validated_clean_spanning_QNAMEs')}/{field(row, 'Validated_all_clipped_QNAMEs')}"
            f"|evidence={field(row, 'False_closure_evidence')}|{field(row, 'False_closure_reason')}"
        )
    return ";".join(parts) if parts else "NA"


def _contig_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--assemblies-fai", required=True)
    parser.add_argument("--contig-circularity", required=True)
    parser.add_argument("--clipping-all-final-tsv", required=True)
    parser.add_argument("--details-tsv", required=True)
    parser.add_argument("--scg-tsv", required=True)
    parser.add_argument("--premature-max-length-bp", type=int, default=500000)
    parser.add_argument("--premature-min-ribosomal-families", type=int, default=3)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.premature_max_length_bp < 1 or args.premature_min_ribosomal_families < 1:
        _contig_fail("premature-circularization thresholds must be positive")

    contig_to_mag, mag_counts = _contig_load_mag_map(args.mag_contig_map)
    lengths = load_fai(args.assemblies_fai)
    circularity = load_circularity(args.contig_circularity)
    if set(circularity) != set(lengths):
        missing = sorted(set(lengths) - set(circularity))[:10]
        extra = sorted(set(circularity) - set(lengths))[:10]
        _contig_fail(f"circularity/all-assemblies FAI contig mismatch: missing={missing}, extra={extra}")
    validated_clipping, unresolved_clipping = load_clipping_state(args.clipping_all_final_tsv)
    details = load_details(args.details_tsv)
    ribosomal_families = load_ribosomal_families(args.scg_tsv)

    header = [
        "Contig", "Source_assembler", "Contig_length_bp",
        "MAG_ID", "IN_FINAL_MAG_Y_N", "MAG_contig_count",
        "SOURCE_CONTIG_ASSEMBLER_CIRCULARITY", "Circularity_resolution_method",
        "FALSE_CLOSURE_ANALYZABLE", "FALSE_CLOSURE",
        "VALIDATED_CLIPPING_PRESENT", "UNRESOLVED_CLIPPING_CANDIDATE_PRESENT",
        "FALSE_CLOSURE_BREAKPOINTS", "INTERNAL_FALSE_CLOSURE_BREAKPOINTS",
        "TERMINAL_FALSE_CLOSURE_BREAKPOINTS", "FALSE_CLOSURE_POSITIONS",
        "Max_left_cross_contig_connected_QNAMEs", "Max_right_cross_contig_connected_QNAMEs",
        "Candidate_additional_target_contigs", "INTERNAL_FALSE_CLOSURE_EVIDENCE",
        "TERMINAL_FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE_EVIDENCE",
        "Ribosomal_protein_family_count", "Ribosomal_protein_families",
        "Small_genome_length_criterion_Y_N", "Small_genome_ribosomal_criterion_Y_N",
        "Paper_small_genome_proxy", "Premature_circularization_candidate",
        "Premature_circularization_screen",
        "Validated_clipping_breakpoints", "Breakpoint_details",
    ]

    output_rows: List[List[object]] = []
    circular_contigs = [
        contig for contig, row in circularity.items()
        if (row.get("SOURCE_CONTIG_ASSEMBLER_CIRCULARITY") or "NA").strip().upper() == "Y"
    ]
    for contig in sorted(circular_contigs, key=lambda c: ((circularity[c].get("Source_assembler") or ""), c)):
        circ_row = circularity[contig]
        rows = details.get(contig, [])
        positive = [r for r in rows if _contig_yes(r, "False_closure_evidence")]
        internal_pos = [r for r in positive if (r.get("Region_type") or "").strip().upper() == "INTERNAL"]
        terminal_pos = [r for r in positive if (r.get("Region_type") or "").strip().upper() == "TERMINAL"]
        unresolved = bool(unresolved_clipping.get(contig, False))
        if positive:
            false_call = "Y"
        elif unresolved:
            false_call = "NA"
        else:
            false_call = "N"
        analyzable = "Y" if false_call in {"Y", "N"} else "N"

        positions = join_unique(
            f"{(r.get('Region_type') or 'NA').strip().upper()}:{(r.get('Breakpoint_1based') or '').strip()}"
            for r in positive
        )
        targets: List[str] = []
        for r in positive:
            targets.extend(target_names((r.get("Left_target_contigs") or "NA").strip()))
            targets.extend(target_names((r.get("Right_target_contigs") or "NA").strip()))
        max_left = max((int(r.get("Left_cross_contig_connected_QNAMEs") or 0) for r in rows), default=0)
        max_right = max((int(r.get("Right_cross_contig_connected_QNAMEs") or 0) for r in rows), default=0)

        families = ribosomal_families.get(contig, ())
        length_ok = lengths[contig] < args.premature_max_length_bp
        ribosomal_ok = len(families) >= args.premature_min_ribosomal_families
        proxy = length_ok and ribosomal_ok
        premature = false_call == "Y"
        screen = premature or proxy
        mag = contig_to_mag.get(contig, "NA")
        in_mag = mag != "NA"

        output_rows.append([
            contig, circ_row["Source_assembler"], lengths[contig],
            mag, "Y" if in_mag else "N", mag_counts.get(mag, 0) if in_mag else 0,
            "Y", circ_row["Circularity_resolution_method"], analyzable, false_call,
            "Y" if validated_clipping.get(contig, False) else "N", "Y" if unresolved else "N",
            len(positive), len(internal_pos), len(terminal_pos), positions,
            max_left, max_right, join_unique(targets),
            "Y" if internal_pos else "N", "Y" if terminal_pos else "N",
            "Y" if positive else "N",
            len(families), ";".join(families) if families else "NA",
            "Y" if length_ok else "N", "Y" if ribosomal_ok else "N",
            "Y" if proxy else "N", "Y" if premature else "N", "Y" if screen else "N",
            len(rows), breakpoint_details(rows),
        ])

    atomic_write(args.output, header, output_rows)
    print(
        f"{_contig_VERSION}: circular_contigs={len(output_rows)}; "
        f"false_closure_Y={sum(row[9] == 'Y' for row in output_rows)}; "
        f"small_genome_proxy_Y={sum(row[header.index('Paper_small_genome_proxy')] == 'Y' for row in output_rows)}"
    )
    print(f"Output: {args.output}")


def cmd_contig(argv):
    """Write the contig table"""
    sys.argv = [sys.argv[0] + " contig", *argv]
    return _contig_main()


COMMANDS = {
    'resolve': cmd_resolve,
    'detect': cmd_detect,
    'mag': cmd_mag,
    'contig': cmd_contig,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
