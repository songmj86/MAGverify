#!/usr/bin/env python3
"""Merge source-assembler-specific assembly-variant results.

Each assembler is profiled independently against its complete co-assembly using
its own HiFi_all_BAM. This helper combines those independent classifications
without changing their reference populations and writes assembly-wide and
assembler-level statistics.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping

POSITIVE = {"Minority_assembly_variant", "Unsupported_assembly_variant"}
COUNT_FIELDS = (
    "minority_assembly_snv_count",
    "unsupported_assembly_snv_count",
    "minority_assembly_indel_count",
    "unsupported_assembly_indel_count",
)
EXAMINED_FIELDS = (
    "anvio_variable_nucleotide_rows_examined",
    "anvio_indel_rows_examined",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--detail-out", required=True)
    p.add_argument("--contig-out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--assembler-summary-out", required=True)
    return p.parse_args()


def read_tsv(path: str) -> tuple[List[str], List[Dict[str, str]]]:
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: TSV lacks a header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def parse_int(value: str, label: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"ERROR: invalid integer for {label}: {value!r}") from exc
    if result < 0:
        raise SystemExit(f"ERROR: negative integer for {label}: {value!r}")
    return result


def atomic_write(path: str, fields: Iterable[str], rows: Iterable[Mapping[str, object]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    field_list = list(fields)
    with open(tmp, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=field_list,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def per_mbp(value: int, bp: int) -> str:
    return f"{value / (bp / 1_000_000.0):.6f}" if bp > 0 else "0.000000"


def summarize_contigs(
    rows: List[Dict[str, str]],
    detail_rows: List[Dict[str, str]],
) -> Dict[str, object]:
    total_bp = sum(parse_int(row["length"], "length") for row in rows)
    totals = {field: sum(parse_int(row[field], field) for row in rows) for field in COUNT_FIELDS}
    examined = {
        field: sum(parse_int(row[field], field) for row in rows) for field in EXAMINED_FIELDS
    }
    total_variants = sum(totals.values())
    assembly_variant_result = "Y" if total_variants > 0 else "N"
    contigs_with = sum((row.get("assembly_allele_support_call") or "") in POSITIVE for row in rows)
    contigs_minority = sum(
        (row.get("assembly_allele_support_call") or "") == "Minority_assembly_variant"
        for row in rows
    )
    contigs_unsupported = sum(
        (row.get("assembly_allele_support_call") or "") == "Unsupported_assembly_variant"
        for row in rows
    )

    minority_genes = set()
    unsupported_genes = set()
    for row in detail_rows:
        gene = (row.get("corresponding_gene_call") or "").strip()
        if gene in {"", "NA", "None", "-1"}:
            continue
        gene_key = ((row.get("Source_assembler") or "").strip(), gene)
        cls = (row.get("assembly_support_class") or "").strip()
        if cls == "Minority_assembly_variant":
            minority_genes.add(gene_key)
        elif cls == "Unsupported_assembly_variant":
            unsupported_genes.add(gene_key)

    return {
        "analyzed_contig_count": len(rows),
        "Assembly_variant": assembly_variant_result,
        "assembly_length_bp": total_bp,
        **examined,
        "minority_assembly_snv_count": totals["minority_assembly_snv_count"],
        "minority_assembly_snv_per_Mbp": per_mbp(totals["minority_assembly_snv_count"], total_bp),
        "unsupported_assembly_snv_count": totals["unsupported_assembly_snv_count"],
        "unsupported_assembly_snv_per_Mbp": per_mbp(totals["unsupported_assembly_snv_count"], total_bp),
        "minority_assembly_indel_count": totals["minority_assembly_indel_count"],
        "minority_assembly_indel_per_Mbp": per_mbp(totals["minority_assembly_indel_count"], total_bp),
        "unsupported_assembly_indel_count": totals["unsupported_assembly_indel_count"],
        "unsupported_assembly_indel_per_Mbp": per_mbp(totals["unsupported_assembly_indel_count"], total_bp),
        "minority_impacted_gene_count": len(minority_genes),
        "unsupported_impacted_gene_count": len(unsupported_genes),
        "contigs_with_assembly_allele_disagreement": contigs_with,
        "contigs_classified_as_minority_assembly_variant": contigs_minority,
        "contigs_classified_as_unsupported_assembly_variant": contigs_unsupported,
        "assembly_variant_count": total_variants,
        "assembly_variant_per_Mbp": per_mbp(total_variants, total_bp),
    }


def main() -> None:
    args = parse_args()
    manifest_fields, manifest_rows = read_tsv(args.manifest)
    required_manifest = {
        "Assembler",
        "Assembly_FASTA",
        "Detail_TSV",
        "Contig_TSV",
        "Summary_TSV",
        "SNV_raw_TSV",
        "INDEL_raw_TSV",
    }
    missing = sorted(required_manifest - set(manifest_fields))
    if missing:
        raise SystemExit("ERROR: per-assembler result manifest lacks columns: " + ", ".join(missing))
    if not manifest_rows:
        raise SystemExit("ERROR: per-assembler result manifest contains no rows")

    all_detail: List[Dict[str, str]] = []
    all_contigs: List[Dict[str, str]] = []
    detail_fields: List[str] | None = None
    contig_fields: List[str] | None = None
    seen_assemblers = set()
    seen_contigs: Dict[str, str] = {}
    seen_events = set()
    per_assembler: List[Dict[str, object]] = []

    for manifest_row in manifest_rows:
        assembler = (manifest_row.get("Assembler") or "").strip()
        if not assembler or assembler in seen_assemblers:
            raise SystemExit(f"ERROR: blank or duplicate assembler in manifest: {assembler!r}")
        seen_assemblers.add(assembler)
        detail_path = manifest_row["Detail_TSV"]
        contig_path = manifest_row["Contig_TSV"]
        summary_path = manifest_row["Summary_TSV"]
        for path in (detail_path, contig_path, summary_path):
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise SystemExit(f"ERROR: missing/empty per-assembler result for {assembler}: {path}")

        d_fields, d_rows = read_tsv(detail_path)
        c_fields, c_rows = read_tsv(contig_path)
        if detail_fields is None:
            detail_fields = d_fields
        elif d_fields != detail_fields:
            raise SystemExit(f"ERROR: detail schema differs for assembler {assembler}")
        if contig_fields is None:
            contig_fields = c_fields
        elif c_fields != contig_fields:
            raise SystemExit(f"ERROR: contig schema differs for assembler {assembler}")

        assembler_detail = []
        for row in d_rows:
            raw_event = (row.get("event_id") or "").strip()
            if not raw_event:
                raise SystemExit(f"ERROR: blank event_id in {detail_path}")
            event_id = f"{assembler}::{raw_event}"
            if event_id in seen_events:
                raise SystemExit(f"ERROR: duplicate merged event_id: {event_id}")
            seen_events.add(event_id)
            merged = {"Source_assembler": assembler, **row, "event_id": event_id}
            all_detail.append(merged)
            assembler_detail.append(merged)

        for row in c_rows:
            contig = (row.get("contig") or "").strip()
            if not contig:
                raise SystemExit(f"ERROR: blank contig in {contig_path}")
            previous = seen_contigs.get(contig)
            if previous is not None:
                raise SystemExit(
                    "ERROR: a contig ID occurs in more than one source assembly and cannot be "
                    f"merged unambiguously: {contig} ({previous}, {assembler})"
                )
            seen_contigs[contig] = assembler
            all_contigs.append({"Source_assembler": assembler, **row})

        stats = summarize_contigs(
            [{"Source_assembler": assembler, **row} for row in c_rows], assembler_detail
        )
        mag_rows = [row for row in c_rows if (row.get("MAG_ID") or "").strip() not in {"", "NA", "None"}]
        mag_bp = sum(parse_int(row["length"], "length") for row in mag_rows)
        mag_variants = sum(
            sum(parse_int(row[field], field) for field in COUNT_FIELDS) for row in mag_rows
        )
        mag_contigs_with = sum(
            (row.get("assembly_allele_support_call") or "") in POSITIVE for row in mag_rows
        )
        per_assembler.append(
            {
                "Assembler": assembler,
                "Assembly_FASTA": os.path.realpath(manifest_row["Assembly_FASTA"]),
                **stats,
                "MAG_assigned_contig_count": len(mag_rows),
                "MAG_assigned_bp": mag_bp,
                "MAG_assigned_contigs_with_assembly_allele_disagreement": mag_contigs_with,
                "MAG_assigned_assembly_variant_count": mag_variants,
                "MAG_assigned_assembly_variant_per_Mbp": per_mbp(mag_variants, mag_bp),
            }
        )

    assert detail_fields is not None and contig_fields is not None
    all_detail.sort(
        key=lambda row: (
            row["Source_assembler"],
            row["contig"],
            int(row["pos_1based"]),
            row["variant_type"],
            row["event_id"],
        )
    )
    all_contigs.sort(key=lambda row: (row["Source_assembler"], row["contig"]))
    per_assembler.sort(key=lambda row: str(row["Assembler"]))

    atomic_write(args.detail_out, ["Source_assembler", *detail_fields], all_detail)
    atomic_write(args.contig_out, ["Source_assembler", *contig_fields], all_contigs)

    overall = summarize_contigs(all_contigs, all_detail)
    summary_fields = list(overall)
    atomic_write(args.summary_out, summary_fields, [overall])

    assembler_fields = [
        "Assembler",
        "Assembly_FASTA",
        *summary_fields,
        "MAG_assigned_contig_count",
        "MAG_assigned_bp",
        "MAG_assigned_contigs_with_assembly_allele_disagreement",
        "MAG_assigned_assembly_variant_count",
        "MAG_assigned_assembly_variant_per_Mbp",
    ]
    atomic_write(args.assembler_summary_out, assembler_fields, per_assembler)

    print(f"Assemblers merged: {len(per_assembler)}")
    print(f"Assembly-wide contigs: {len(all_contigs)}")
    print(f"Assembly-wide disagreement events: {len(all_detail)}")
    print(f"Assembly-wide sequence length: {overall['assembly_length_bp']} bp")


if __name__ == "__main__":
    main()
