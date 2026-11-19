#!/usr/bin/env python3
"""Create the assembly-wide contig-level false-closure summary.

All assembler-reported circular source contigs are listed, irrespective of MAG
membership. False closure is called from HiFi-supported clipping breakpoints;
MAG membership is added only as a downstream annotation.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from false_closure_scg_utils import load_ribosomal_families

VERSION = "summarize_false_closure_by_contig"


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def pick_field(fieldnames: List[str], candidates: Tuple[str, ...], label: str) -> str:
    lower = {name.strip().lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    fail(f"missing {label}; observed columns: {fieldnames}")


def load_mag_map(path: str) -> Tuple[Dict[str, str], Dict[str, int]]:
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
                fail(f"blank MAG/contig at {path}:{line_no}")
            previous = contig_to_mag.get(contig)
            if previous is not None and previous != mag:
                fail(f"contig assigned to multiple MAGs: {contig}")
            if previous is None:
                contig_to_mag[contig] = mag
                counts[mag] += 1
    if not contig_to_mag:
        fail(f"MAG contig map contains no rows: {path}")
    return contig_to_mag, counts


def load_fai(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                fail(f"malformed FAI row at {path}:{line_no}")
            contig = fields[0]
            try:
                length = int(fields[1])
            except ValueError as exc:
                fail(f"invalid FAI length at {path}:{line_no}: {exc}")
            if not contig or length < 1 or contig in lengths:
                fail(f"blank/duplicate FAI contig: {contig}")
            lengths[contig] = length
    if not lengths:
        fail(f"FAI contains no contigs: {path}")
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
            fail(f"resolved circularity table lacks {sorted(missing)}")
        for line_no, row in enumerate(reader, start=2):
            contig = (row.get("Contig") or "").strip()
            call = (row.get("SOURCE_CONTIG_ASSEMBLER_CIRCULARITY") or "NA").strip().upper()
            if not contig or contig in out or call not in {"Y", "N", "NA"}:
                fail(f"invalid/duplicate circularity row at {path}:{line_no}: {contig}/{call}")
            out[contig] = row
    if not out:
        fail(f"resolved circularity table contains no contigs: {path}")
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
            fail(f"clipping table lacks columns {sorted(missing)}")
        for line_no, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            analyzable = (row.get("Clipping_analyzable") or "").strip().upper()
            call = (row.get("Clipping") or "NA").strip().upper()
            if not contig or analyzable not in {"Y", "N"}:
                fail(f"invalid clipping row at {path}:{line_no}")
            if analyzable == "N":
                if call != "NA":
                    fail(f"Clipping_analyzable=N requires Clipping=NA at {path}:{line_no}")
                unresolved[contig] = True
            elif call == "Y":
                validated[contig] = True
            elif call != "N":
                fail(f"Clipping_analyzable=Y requires Clipping=Y/N at {path}:{line_no}")
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
            fail(f"false-closure breakpoint table lacks columns {sorted(missing)}")
        for row in reader:
            contig = (row.get("Contig") or "").strip()
            if contig:
                by_contig[contig].append(row)
    return by_contig


def yes(row: dict, field: str) -> bool:
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
                    fail(f"schema mismatch at output line {line_no}: {len(row)} != {len(header)}")
                writer.writerow(row)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--combined-fai", required=True)
    parser.add_argument("--contig-circularity", required=True)
    parser.add_argument("--clipping-all-final-tsv", required=True)
    parser.add_argument("--details-tsv", required=True)
    parser.add_argument("--scg-tsv", required=True)
    parser.add_argument("--premature-max-length-bp", type=int, default=500000)
    parser.add_argument("--premature-min-ribosomal-families", type=int, default=3)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.premature_max_length_bp < 1 or args.premature_min_ribosomal_families < 1:
        fail("premature-circularization thresholds must be positive")

    contig_to_mag, mag_counts = load_mag_map(args.mag_contig_map)
    lengths = load_fai(args.combined_fai)
    circularity = load_circularity(args.contig_circularity)
    if set(circularity) != set(lengths):
        missing = sorted(set(lengths) - set(circularity))[:10]
        extra = sorted(set(circularity) - set(lengths))[:10]
        fail(f"circularity/combined-FAI contig mismatch: missing={missing}, extra={extra}")
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
    ]

    output_rows: List[List[object]] = []
    circular_contigs = [
        contig for contig, row in circularity.items()
        if (row.get("SOURCE_CONTIG_ASSEMBLER_CIRCULARITY") or "NA").strip().upper() == "Y"
    ]
    for contig in sorted(circular_contigs, key=lambda c: ((circularity[c].get("Source_assembler") or ""), c)):
        circ_row = circularity[contig]
        rows = details.get(contig, [])
        positive = [r for r in rows if yes(r, "False_closure_evidence")]
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
        ])

    atomic_write(args.output, header, output_rows)
    print(
        f"{VERSION}: circular_contigs={len(output_rows)}; "
        f"false_closure_Y={sum(row[9] == 'Y' for row in output_rows)}; "
        f"small_genome_proxy_Y={sum(row[header.index('Paper_small_genome_proxy')] == 'Y' for row in output_rows)}"
    )
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
