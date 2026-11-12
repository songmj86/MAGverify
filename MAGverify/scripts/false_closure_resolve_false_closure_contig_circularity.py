#!/usr/bin/env python3
"""Resolve assembler-reported circularity at source-contig level.

Priority:
  1. assembler-native contig metadata
     - metaMDBG FASTA header: circular=yes/no
     - hifiasm-meta contig suffix: ...c / ...l
     - metaFlye assembly_info.txt: circ. Y/N
  2. MAG-level circularity only for single-contig MAGs (identity unambiguous)

Multi-contig MAGs are never assigned a contig circularity from a MAG-level Y/N
alone. This prevents a catalogue assembled from multiple source assemblers from
projecting one MAG-level circularity call onto every contig in the bin.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, NoReturn, Sequence, Tuple

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from classify_circularity import resolve_call_table

VERSION = "resolve_false_closure_contig_circularity"
HIFIASM_META_RE = re.compile(r"^s[0-9]+\.[uc]tg[0-9]{6}([lc])$")


def fail(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


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
        # metaMDBG (and any compatible header): circular=yes/no
        match = re.search(r"(?:^|\s)circular\s*=\s*([^\s]+)", header, flags=re.I)
        if match:
            call = normalize_call(match.group(1))
            if call is not None:
                out[contig] = (call, "SOURCE_FASTA_HEADER_CIRCULAR_FIELD")
                continue
        # hifiasm-meta official contig convention ends in l/c.
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
                # Standard Flye assembly_info columns, even if header was removed.
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
    # deduplicate while preserving order
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
                fail(f"source assembly FASTA missing for {assembler}: {fasta}")
            fasta_of[assembler] = fasta
            header_calls = parse_header_circularity(fasta, assembler)
            names = []
            for contig, _header in fasta_headers(fasta):
                names.append(contig)
                prev = source_of.get(contig)
                if prev is not None and prev != assembler:
                    fail(f"contig occurs in multiple source assemblers: {contig} ({prev}, {assembler})")
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
                            "multiple Flye assembly_info tables tie for source-contig overlap: "
                            + ", ".join(x[1] for x in best)
                        )
                    for contig, value in best[0][2].items():
                        if contig in name_set:
                            calls[contig] = value
    if not source_of:
        fail("mapping manifest/source FASTAs contain no contigs")
    return source_of, calls, fasta_of


def atomic_tsv(path: str, header: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
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


def main() -> None:
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
        fail(f"MAG contigs absent from all source assemblies: {missing_mag_contigs[:10]}")

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

    # The MAG-level table is compatibility/audit metadata only. For a multi-contig
    # MAG its semantics cannot be assumed to mean that every (or any particular)
    # source contig is circular, so no aggregate MAG-vs-contig consistency rule is
    # enforced. The only hard conflict check is for a single-contig MAG, where the
    # MAG and source-contig identities are unambiguous.

    atomic_tsv(
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
            "source-contig circularity conflicts with MAG-level circularity for single-contig MAGs: "
            + ", ".join(single_conflicts[:10])
        )
    resolved = sum(1 for row in rows if row[5] in {"Y", "N"})
    circular = sum(1 for row in rows if row[5] == "Y")
    print(f"{VERSION}: contigs={len(rows)} resolved={resolved} circular={circular} unresolved={len(rows)-resolved}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
