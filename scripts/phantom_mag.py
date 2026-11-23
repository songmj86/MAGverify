#!/usr/bin/env python3
"""Create MAG-associated phantom-region, MAG-level and contig-level tables.

This script does not recalculate coverage, BLAST, or k-mer evidence. It joins an
existing ``phantom.region.tsv`` table to the MAG-contig ownership map
and summarizes the retained phantom candidates across every mapped MAG.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


PHANTOM_REQUIRED = (
    "region_id",
    "contig",
    "start_0based",
    "end_0based",
    "length_bp",
    "raw_read_blast_hsps",
    "best_hit_pident",
    "best_hit_aligned_bp",
    "raw_read_blast_query_covered_bp",
    "raw_read_blast_query_coverage",
    "phantom_sequence_call",
)

SITE_HEADER = ("MAG_ID",) + PHANTOM_REQUIRED
SUMMARY_HEADER = (
    "MAG_ID",
    "n_contigs",
    "total_length_bp",
    "n_phantom_regions",
    "n_contigs_with_phantom",
    "phantom_bp",
    "phantom_fraction",
    "max_phantom_region_bp",
    "phantom_contigs",
    "phantom_region_ids",
    "LR_PHANTOM_ANALYZABLE",
    "LR_PHANTOM",
)
CONTIG_HEADER = (
    "Contig",
    "MAG_ID",
    "length_bp",
    "n_phantom_regions",
    "phantom_bp",
    "phantom_fraction",
    "max_phantom_region_bp",
    "phantom_regions_0based",
    "phantom_region_ids",
    "LR_PHANTOM",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate MAG phantom outputs from an existing "
            "phantom.region.tsv. No mapping, BLAST, or Meryl analysis "
            "is rerun."
        )
    )
    parser.add_argument("--phantom", required=True, help="Existing phantom candidate TSV")
    parser.add_argument("--mag-contig-map", required=True, help="MAG-to-contig TSV")
    parser.add_argument("--fai", required=True, help="Assembly FASTA index (.fai)")
    parser.add_argument("--site-output", required=True, help="MAG-associated region TSV")
    parser.add_argument("--summary-output", required=True, help="MAG-level summary TSV")
    parser.add_argument("--contig-output", help="Contig-level table: one row per assembly contig")
    return parser.parse_args()


def require_nonempty(path: str, label: str) -> None:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"ERROR: {label} not found or empty: {path}")


def resolve_column(fieldnames: Sequence[str], candidates: Iterable[str], label: str) -> str:
    normalized = {name.strip().lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    raise SystemExit(
        f"ERROR: MAG_CONTIG_MAP needs a {label} column; observed columns: {fieldnames}"
    )


def load_fai(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                raise SystemExit(f"ERROR: malformed FAI row at line {line_number}: {path}")
            contig = fields[0].strip()
            try:
                length = int(fields[1])
            except ValueError as error:
                raise SystemExit(
                    f"ERROR: invalid FAI length at line {line_number}: {fields[1]}"
                ) from error
            if not contig or length < 1 or contig in lengths:
                raise SystemExit(f"ERROR: blank, duplicate, or invalid FAI contig: {contig}")
            lengths[contig] = length
    if not lengths:
        raise SystemExit(f"ERROR: FAI contains no contigs: {path}")
    return lengths


def load_mag_map(
    path: str, lengths: Mapping[str, int]
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    contig_to_mag: Dict[str, str] = {}
    mag_to_contigs: Dict[str, List[str]] = defaultdict(list)
    seen_rows = set()

    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        mag_column = resolve_column(
            fields,
            ("mag_id", "entity_id", "genome_id", "bin_id", "genome", "bin", "id"),
            "MAG ID",
        )
        contig_column = resolve_column(fields, ("contig", "original_contig"), "Contig")
        for line_number, row in enumerate(reader, start=2):
            mag_id = (row.get(mag_column) or "").strip()
            contig = (row.get(contig_column) or "").strip()
            if not mag_id or not contig:
                raise SystemExit(
                    f"ERROR: blank MAG/contig in MAG_CONTIG_MAP at line {line_number}"
                )
            key = (mag_id, contig)
            if key in seen_rows:
                raise SystemExit(f"ERROR: duplicate MAG-map row: {mag_id}/{contig}")
            seen_rows.add(key)
            if contig not in lengths:
                raise SystemExit(f"ERROR: MAG-map contig absent from FAI: {contig}")
            previous = contig_to_mag.get(contig)
            if previous is not None and previous != mag_id:
                raise SystemExit(
                    f"ERROR: contig assigned to multiple MAGs: {contig} ({previous}, {mag_id})"
                )
            if previous is None:
                contig_to_mag[contig] = mag_id
                mag_to_contigs[mag_id].append(contig)

    if not mag_to_contigs:
        raise SystemExit(f"ERROR: MAG_CONTIG_MAP contains no assignments: {path}")
    for mag_id in mag_to_contigs:
        mag_to_contigs[mag_id].sort()
    return contig_to_mag, dict(mag_to_contigs)


def parse_int(value: str, field: str, line_number: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise SystemExit(
            f"ERROR: invalid integer in {field} at phantom line {line_number}: {value}"
        ) from error


def parse_float(value: str, field: str, line_number: int) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise SystemExit(
            f"ERROR: invalid number in {field} at phantom line {line_number}: {value}"
        ) from error
    if not math.isfinite(result):
        raise SystemExit(
            f"ERROR: non-finite number in {field} at phantom line {line_number}: {value}"
        )
    return result


def load_phantom(path: str, lengths: Mapping[str, int]) -> List[Dict[str, str]]:
    records: List[Dict[str, str]] = []
    seen_region_ids = set()
    seen_intervals = set()

    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [column for column in PHANTOM_REQUIRED if column not in fields]
        if missing:
            raise SystemExit(f"ERROR: phantom table lacks columns: {missing}")

        for line_number, row in enumerate(reader, start=2):
            clean = {column: (row.get(column) or "").strip() for column in PHANTOM_REQUIRED}
            region_id = clean["region_id"]
            contig = clean["contig"]
            if not region_id or not contig:
                raise SystemExit(f"ERROR: blank region_id/contig at phantom line {line_number}")
            if region_id in seen_region_ids:
                raise SystemExit(f"ERROR: duplicate phantom region_id: {region_id}")
            seen_region_ids.add(region_id)
            if contig not in lengths:
                raise SystemExit(
                    f"ERROR: phantom contig absent from assembly FAI at line "
                    f"{line_number}: {contig}"
                )

            start = parse_int(clean["start_0based"], "start_0based", line_number)
            end = parse_int(clean["end_0based"], "end_0based", line_number)
            length = parse_int(clean["length_bp"], "length_bp", line_number)
            blast_hsps = parse_int(
                clean["raw_read_blast_hsps"], "raw_read_blast_hsps", line_number
            )
            best_hit_aligned = parse_int(
                clean["best_hit_aligned_bp"], "best_hit_aligned_bp", line_number
            )
            query_covered = parse_int(
                clean["raw_read_blast_query_covered_bp"],
                "raw_read_blast_query_covered_bp",
                line_number,
            )
            query_coverage = parse_float(
                clean["raw_read_blast_query_coverage"],
                "raw_read_blast_query_coverage",
                line_number,
            )

            if not (0 <= start < end <= lengths[contig]):
                raise SystemExit(
                    f"ERROR: phantom coordinates out of range at line {line_number}: "
                    f"{contig}:{start}-{end}"
                )
            if end - start != length or length < 1:
                raise SystemExit(
                    f"ERROR: phantom length disagrees with coordinates at line {line_number}"
                )
            if blast_hsps != 0 or best_hit_aligned != 0 or query_covered != 0:
                raise SystemExit(
                    f"ERROR: candidate-only phantom table contains raw-read alignment "
                    f"support at line {line_number}"
                )
            if not math.isclose(query_coverage, 0.0, rel_tol=0.0, abs_tol=5e-7):
                raise SystemExit(
                    f"ERROR: candidate-only phantom query coverage is nonzero at line "
                    f"{line_number}"
                )
            if clean["phantom_sequence_call"] != "No_raw_read_BLAST_hit_phantom_candidate":
                raise SystemExit(
                    f"ERROR: unexpected phantom_sequence_call at line {line_number}: "
                    f"{clean['phantom_sequence_call']}"
                )

            interval = (contig, start, end)
            if interval in seen_intervals:
                raise SystemExit(f"ERROR: duplicate phantom interval: {contig}:{start}-{end}")
            seen_intervals.add(interval)
            records.append(clean)

    records.sort(
        key=lambda row: (
            row["contig"],
            int(row["start_0based"]),
            int(row["end_0based"]),
            row["region_id"],
        )
    )
    return records


def atomic_write_rows(path: str, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    output = Path(path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            delete=False,
            dir=str(output.parent),
            prefix=output.name + ".",
            suffix=".tmp",
        ) as handle:
            temp_name = handle.name
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
        os.replace(temp_name, output)
        temp_name = None
    finally:
        if temp_name and os.path.exists(temp_name):
            os.unlink(temp_name)


def main() -> None:
    args = parse_args()
    require_nonempty(args.phantom, "phantom candidate table")
    require_nonempty(args.mag_contig_map, "MAG-contig map")
    require_nonempty(args.fai, "assembly FAI")

    lengths = load_fai(args.fai)
    contig_to_mag, mag_to_contigs = load_mag_map(args.mag_contig_map, lengths)
    phantom_records = load_phantom(args.phantom, lengths)

    mag_sites: List[Tuple[object, ...]] = []
    by_mag: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for record in phantom_records:
        mag_id = contig_to_mag.get(record["contig"])
        if mag_id is None:
            continue
        by_mag[mag_id].append(record)
        mag_sites.append((mag_id,) + tuple(record[column] for column in PHANTOM_REQUIRED))

    mag_sites.sort(
        key=lambda row: (row[0], row[2], int(row[3]), int(row[4]), row[1])
    )
    atomic_write_rows(args.site_output, SITE_HEADER, mag_sites)

    summary_rows: List[Tuple[object, ...]] = []
    for mag_id in sorted(mag_to_contigs):
        members = mag_to_contigs[mag_id]
        total_length = sum(lengths[contig] for contig in members)
        records = sorted(
            by_mag.get(mag_id, []),
            key=lambda row: (
                row["contig"],
                int(row["start_0based"]),
                int(row["end_0based"]),
                row["region_id"],
            ),
        )
        phantom_bp = sum(int(row["length_bp"]) for row in records)
        phantom_contigs = sorted({row["contig"] for row in records})
        region_ids = [row["region_id"] for row in records]
        maximum = max((int(row["length_bp"]) for row in records), default=0)
        fraction = phantom_bp / total_length if total_length else 0.0
        summary_rows.append(
            (
                mag_id,
                len(members),
                total_length,
                len(records),
                len(phantom_contigs),
                phantom_bp,
                f"{fraction:.6f}",
                maximum,
                ",".join(phantom_contigs) if phantom_contigs else "NA",
                ",".join(region_ids) if region_ids else "NA",
                "Y",
                "Y" if records else "N",
            )
        )

    atomic_write_rows(args.summary_output, SUMMARY_HEADER, summary_rows)

    if args.contig_output:
        by_contig: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        for record in phantom_records:
            by_contig[record["contig"]].append(record)
        contig_rows: List[Tuple[object, ...]] = []
        for contig, length in lengths.items():
            records = by_contig.get(contig, [])
            phantom_bp = sum(int(row["length_bp"]) for row in records)
            contig_rows.append((
                contig,
                contig_to_mag.get(contig, "NA"),
                length,
                len(records),
                phantom_bp,
                f"{phantom_bp / length:.6f}" if length else "0.000000",
                max((int(row["length_bp"]) for row in records), default=0),
                ",".join(f"{row['start_0based']}-{row['end_0based']}" for row in records) or "NA",
                ",".join(row["region_id"] for row in records) or "NA",
                "Y" if records else "N",
            ))
        atomic_write_rows(args.contig_output, CONTIG_HEADER, contig_rows)
        print(f"Contigs with phantom regions: {sum(1 for row in contig_rows if row[-1] == 'Y')}")

    print(f"Input phantom regions: {len(phantom_records)}")
    print(f"MAG-associated phantom regions: {len(mag_sites)}")
    print(f"MAGs summarized: {len(summary_rows)}")
    print(f"MAGs with phantom regions: {sum(1 for row in summary_rows if row[11] == 'Y')}")
    print(f"MAG-associated phantom bp: {sum(int(row[5]) for row in summary_rows)}")


if __name__ == "__main__":
    main()
