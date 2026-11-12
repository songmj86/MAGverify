#!/usr/bin/env python3
"""Filter a protein FASTA to genes from contigs with at least N predicted genes.

The coordinate table is expected to be the four-column table emitted by
chimera_link_prodigal_accessions.py and to preserve protein-FASTA record order. The
script validates accession order while streaming so it does not need to keep
all gene IDs in memory.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def count_genes_by_contig(coord_path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with coord_path.open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) != 4:
                fail(f"malformed coordinate row {line_no}: expected 4 columns")
            gene, contig, start, end = fields
            if not gene or not contig:
                fail(f"blank gene/contig in coordinate row {line_no}")
            try:
                start_i, end_i = int(start), int(end)
            except ValueError as exc:
                fail(f"invalid coordinates in row {line_no}: {start}-{end}")
            if start_i < 1 or end_i < start_i:
                fail(f"invalid coordinates in row {line_no}: {start}-{end}")
            counts[contig] += 1
    if not counts:
        fail(f"coordinate table is empty: {coord_path}")
    return counts


def next_coordinate(handle, line_no_ref: list[int]):
    for raw in handle:
        line_no_ref[0] += 1
        if not raw.strip():
            continue
        fields = raw.rstrip("\n\r").split("\t")
        if len(fields) != 4:
            fail(f"malformed coordinate row {line_no_ref[0]} during streaming")
        return fields
    return None


def filter_faa(faa: Path, coords: Path, output: Path, min_genes: int) -> dict[str, int]:
    if min_genes < 1:
        fail("min_genes must be >=1")
    counts = count_genes_by_contig(coords)
    eligible = {contig for contig, n in counts.items() if n >= min_genes}

    output.parent.mkdir(parents=True, exist_ok=True)
    total_genes = selected_genes = 0
    total_contigs = len(counts)
    selected_contigs = len(eligible)

    with coords.open(encoding="utf-8") as coord_handle, faa.open(encoding="utf-8") as faa_handle, output.open("w", encoding="utf-8", newline="") as out:
        coord_line_no = [0]
        write_record = False
        for faa_line_no, raw in enumerate(faa_handle, 1):
            if raw.startswith(">"):
                fields = next_coordinate(coord_handle, coord_line_no)
                if fields is None:
                    fail(f"FASTA contains more records than coordinate table; first extra FASTA record at line {faa_line_no}")
                gene, contig, _start, _end = fields
                accession = raw[1:].split(None, 1)[0].strip()
                if accession != gene:
                    fail(
                        "FASTA/coordinate accession order mismatch: "
                        f"FASTA={accession!r}, coordinates={gene!r} at coordinate row {coord_line_no[0]}"
                    )
                total_genes += 1
                write_record = contig in eligible
                if write_record:
                    selected_genes += 1
                    out.write(raw)
            elif write_record:
                out.write(raw)

        extra = next_coordinate(coord_handle, coord_line_no)
        if extra is not None:
            fail(f"coordinate table contains more genes than FASTA; first extra gene={extra[0]!r}")

    if total_genes != sum(counts.values()):
        fail(f"FASTA/coordinate gene counts disagree: {total_genes} versus {sum(counts.values())}")

    return {
        "min_genes_per_contig": min_genes,
        "total_contigs": total_contigs,
        "selected_contigs": selected_contigs,
        "skipped_contigs": total_contigs - selected_contigs,
        "total_genes": total_genes,
        "selected_genes": selected_genes,
        "skipped_genes": total_genes - selected_genes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("faa", type=Path)
    parser.add_argument("coordinates", type=Path)
    parser.add_argument("output_faa", type=Path)
    parser.add_argument("stats_json", type=Path)
    parser.add_argument("min_genes", type=int)
    args = parser.parse_args()

    stats = filter_faa(args.faa, args.coordinates, args.output_faa, args.min_genes)
    args.stats_json.parent.mkdir(parents=True, exist_ok=True)
    args.stats_json.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    print(
        "Taxonomy prefilter: "
        f"{stats['selected_contigs']}/{stats['total_contigs']} contigs and "
        f"{stats['selected_genes']}/{stats['total_genes']} genes retained "
        f"(minimum {stats['min_genes_per_contig']} predicted genes per contig)"
    )


if __name__ == "__main__":
    main()
