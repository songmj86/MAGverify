#!/usr/bin/env python3
"""Refresh MAG-level circularity compatibility outputs from canonical HiFi results.

This script performs no read mapping and no junction re-analysis. It only
rebuilds MAG_circularity.final.tsv, and optionally MAG_circularity.audit.tsv,
from MAG_END_START_closure.final.tsv plus assembler-reported circularity calls.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, NoReturn, Optional, Sequence, Set, Tuple


def fail(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


TRUE_VALUES = {"Y", "YES", "TRUE", "T", "1", "CIRCULAR", "CIRCULAR=Y"}
FALSE_VALUES = {"N", "NO", "FALSE", "F", "0", "LINEAR", "NOT_CIRCULAR"}
ID_COLUMNS = ("id", "mag_id", "genome", "genome_id", "bin", "bin_id")
CALL_COLUMNS = ("assembler_circularity", "circularity", "circular", "is_circular")


def canonical_id(value: str) -> str:
    value = str(value).strip().strip('"').strip("'")
    value = re.split(r"[/\\]", value)[-1]
    lower = value.lower()
    for suffix in (".fasta", ".fna", ".fa"):
        if lower.endswith(suffix):
            return value[: -len(suffix)]
    return value


def normalize_call(value: str) -> Optional[str]:
    upper = str(value).strip().upper()
    if upper in TRUE_VALUES:
        return "Y"
    if upper in FALSE_VALUES:
        return "N"
    return None


def table_header(path: str) -> Optional[Tuple[list[str], int, int]]:
    try:
        with open(path, encoding="utf-8-sig") as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                header = raw.rstrip("\n\r").split("\t")
                lowered = [field.strip().lower() for field in header]
                id_index = next((lowered.index(x) for x in ID_COLUMNS if x in lowered), None)
                call_index = next((lowered.index(x) for x in CALL_COLUMNS if x in lowered), None)
                if id_index is None:
                    return None
                if call_index is not None:
                    return header, id_index, call_index
                if "assembler" not in lowered:
                    return None
                assembler_index = lowered.index("assembler")
                saw_binary = False
                for data_raw in handle:
                    if not data_raw.strip() or data_raw.startswith("#"):
                        continue
                    fields = data_raw.rstrip("\n\r").split("\t")
                    if len(fields) <= assembler_index or not fields[assembler_index].strip():
                        continue
                    if normalize_call(fields[assembler_index]) is None:
                        return None
                    saw_binary = True
                return (header, id_index, assembler_index) if saw_binary else None
    except (OSError, UnicodeDecodeError):
        return None
    return None


def parse_call_table(path: str) -> Dict[str, str]:
    header_info = table_header(path)
    if header_info is None:
        fail(f"no recognized ID/circularity header in {path}")
    _header, id_index, call_index = header_info
    calls: Dict[str, str] = {}
    header_seen = False
    with open(path, encoding="utf-8-sig") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if not header_seen:
                header_seen = True
                continue
            if len(fields) <= max(id_index, call_index):
                continue
            mag_id = canonical_id(fields[id_index])
            call = normalize_call(fields[call_index])
            if not mag_id or call is None:
                continue
            if mag_id in calls and calls[mag_id] != call:
                fail(f"conflicting circularity calls for {mag_id} in {path}")
            calls[mag_id] = call
    return calls


def resolve_call_table(path: str, mag_ids: Set[str]) -> Tuple[str, Dict[str, str]]:
    candidate_path = Path(path)
    if candidate_path.is_file():
        return str(candidate_path), parse_call_table(str(candidate_path))
    if not candidate_path.is_dir():
        fail(f"CIRCULARITY_TSV is neither a file nor directory: {path}")
    scored: list[Tuple[int, int, str, Dict[str, str]]] = []
    for candidate in sorted(candidate_path.rglob("*")):
        if not candidate.is_file():
            continue
        if table_header(str(candidate)) is None:
            continue
        calls = parse_call_table(str(candidate))
        scored.append((len(mag_ids & set(calls)), len(calls), str(candidate), calls))
    if not scored:
        fail(f"no table with recognized circularity calls under {path}")
    scored.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
    best_overlap, best_count = scored[0][0], scored[0][1]
    if best_overlap == 0:
        fail(f"circularity tables under {path} do not overlap target MAG IDs")
    tied = [x for x in scored if x[0] == best_overlap and x[1] == best_count]
    if len(tied) > 1:
        fail("multiple circularity tables tie for best MAG-ID overlap")
    return scored[0][2], scored[0][3]

def read_tsv(path: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not fields:
        fail(f"empty TSV header: {path}")
    return fields, rows


def atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        prefix=os.path.basename(path) + ".", suffix=".tmp", dir=directory, text=True
    )
    os.close(fd)
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def load_assembler_calls(circularity_source: str, mag_ids: set[str]) -> Dict[str, str]:
    _selected_call_table, assembler_calls = resolve_call_table(circularity_source, mag_ids)
    missing = sorted(mag_ids - set(assembler_calls))
    if missing:
        preview = ", ".join(missing[:10])
        fail(
            f"assembler circularity call is missing for {len(missing)} MAGs; "
            f"examples: {preview}"
        )
    invalid = {mag: call for mag, call in assembler_calls.items() if mag in mag_ids and call not in {"Y", "N"}}
    if invalid:
        fail(f"invalid assembler circularity calls: {invalid}")
    return assembler_calls


def build_rows(
    canonical_rows: Sequence[Mapping[str, str]], assembler_calls: Mapping[str, str]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    final_rows: List[dict[str, str]] = []
    audit_rows: List[dict[str, str]] = []
    seen: set[str] = set()

    for row in canonical_rows:
        mag = row["MAG_ID"]
        if mag in seen:
            fail(f"duplicate MAG_ID in canonical circularity output: {mag}")
        seen.add(mag)

        try:
            n = int(row["Number_of_contigs"])
        except ValueError:
            fail(f"invalid Number_of_contigs for {mag}: {row['Number_of_contigs']}")
        if n < 1:
            fail(f"Number_of_contigs must be >=1 for {mag}: {n}")

        eligibility = "Y_single_contig" if n == 1 else "N_not_single_contig"
        canonical_analyzable = row.get("CIRCULARITY_ANALYZABLE", "")
        canonical_topology = row.get("HIFI_END_START_SUPPORTED", "")

        if canonical_analyzable not in {"Y", "N"}:
            fail(f"invalid CIRCULARITY_ANALYZABLE for {mag}: {canonical_analyzable}")
        if canonical_topology not in {"Y", "N", "NA"}:
            fail(f"invalid HIFI_END_START_SUPPORTED for {mag}: {canonical_topology}")

        if n != 1:
            if canonical_analyzable != "N" or canonical_topology != "NA":
                fail(
                    f"multi-contig MAG must be non-analyzable with topology NA: {mag} "
                    f"({canonical_analyzable}, {canonical_topology})"
                )
            evaluable = "N"
            topology = "NA"
        else:
            evaluable = canonical_analyzable
            if evaluable == "Y":
                if canonical_topology not in {"Y", "N"}:
                    fail(f"evaluable single-contig MAG lacks Y/N topology call: {mag}")
                topology = canonical_topology
            else:
                if canonical_topology != "NA":
                    fail(f"non-evaluable single-contig MAG must have topology NA: {mag}")
                topology = "NA"

        assembler = assembler_calls[mag]
        if assembler == "NA":
            high_confidence = additional_candidate = "NA"
            circular = topology
        else:
            high_confidence = "Y" if assembler == "Y" and topology == "Y" else "N"
            additional_candidate = "Y" if assembler == "N" and topology == "Y" else "N"
            circular = high_confidence

        final_rows.append(
            {
                "ID": mag,
                "Assembler_reported_circularity": assembler,
                "Circularity_evaluable": evaluable,
                "Circularity_eligibility": eligibility,
                "Read_supported_circular_topology": topology,
                "High_confidence_circular": high_confidence,
                "Additional_circular_candidate": additional_candidate,
                "Circular": circular,
            }
        )
        audit_rows.append(
            {"ID": mag, "Contig_count": str(n), "Eligibility": eligibility}
        )

    return final_rows, audit_rows


def validate_expected_structure(rows: Sequence[Mapping[str, str]]) -> None:
    high_conf = Counter(row["High_confidence_circular"] for row in rows)
    candidates = Counter(row["Additional_circular_candidate"] for row in rows)
    both = [row["ID"] for row in rows if row["High_confidence_circular"] == "Y" and row["Additional_circular_candidate"] == "Y"]
    if both:
        fail(f"MAG cannot be both high-confidence circular and additional candidate: {both[:10]}")
    if high_conf.get("Y", 0) + candidates.get("Y", 0) > sum(row["Read_supported_circular_topology"] == "Y" for row in rows):
        fail("derived circular classes exceed read-supported topology calls")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical", required=True, help="MAG_END_START_closure.final.tsv")
    parser.add_argument("--assembler-circularity", default="",
                        help="Assembler circularity TSV or directory; if empty, calls are NA")
    parser.add_argument("--final", required=True, help="MAG_circularity.final.tsv")
    parser.add_argument("--audit", help="Optional MAG_circularity.audit.tsv")
    args = parser.parse_args()

    fields, canonical_rows = read_tsv(args.canonical)
    required = {
        "MAG_ID", "Number_of_contigs", "CIRCULARITY_ANALYZABLE",
        "HIFI_END_START_SUPPORTED",
    }
    missing = required - set(fields)
    if missing:
        fail(f"canonical circularity output lacks {sorted(missing)}: {args.canonical}")
    if not canonical_rows:
        fail(f"canonical circularity output has no data rows: {args.canonical}")

    mag_ids = {row["MAG_ID"] for row in canonical_rows}
    if len(mag_ids) != len(canonical_rows):
        fail("canonical circularity output contains duplicate MAG_ID values")

    if args.assembler_circularity:
        assembler_calls = load_assembler_calls(args.assembler_circularity, mag_ids)
    else:
        assembler_calls = {mag: "NA" for mag in mag_ids}
    final_rows, audit_rows = build_rows(canonical_rows, assembler_calls)
    validate_expected_structure(final_rows)

    final_fields = [
        "ID", "Assembler_reported_circularity", "Circularity_evaluable",
        "Circularity_eligibility", "Read_supported_circular_topology",
        "High_confidence_circular", "Additional_circular_candidate", "Circular",
    ]
    atomic_write(args.final, final_fields, final_rows)
    if args.audit:
        atomic_write(args.audit, ["ID", "Contig_count", "Eligibility"], audit_rows)

    counts = {
        "total": len(final_rows),
        "single_contig": sum(row["Circularity_eligibility"] == "Y_single_contig" for row in final_rows),
        "evaluable": sum(row["Circularity_evaluable"] == "Y" for row in final_rows),
        "topology_Y": sum(row["Read_supported_circular_topology"] == "Y" for row in final_rows),
        "topology_N": sum(row["Read_supported_circular_topology"] == "N" for row in final_rows),
        "topology_NA": sum(row["Read_supported_circular_topology"] == "NA" for row in final_rows),
        "assembler_Y": sum(row["Assembler_reported_circularity"] == "Y" for row in final_rows),
        "high_confidence": sum(row["High_confidence_circular"] == "Y" for row in final_rows),
        "additional_candidates": sum(row["Additional_circular_candidate"] == "Y" for row in final_rows),
    }
    print("Circularity compatibility output refreshed:")
    for key, value in counts.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
