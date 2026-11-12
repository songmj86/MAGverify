#!/usr/bin/env python3
"""Merge non-overlapping nine-column gene-taxonomy tables."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

HEADER = [
    "contig", "start", "end", "domain", "phylum", "class", "order",
    "family", "gene_id",
]


def load(path: Path):
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != HEADER:
            raise SystemExit(
                f"ERROR: unexpected taxonomy header in {path}: {reader.fieldnames!r}; "
                f"expected {HEADER!r}"
            )
        for line_no, row in enumerate(reader, 2):
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError as error:
                raise SystemExit(f"ERROR: invalid coordinates in {path}:{line_no}") from error
            if not row["contig"] or not row["gene_id"] or start < 1 or end < start:
                raise SystemExit(f"ERROR: invalid taxonomy row in {path}:{line_no}")
            rows.append(
                (
                    row["contig"], start, end, row["domain"], row["phylum"],
                    row["class"], row["order"], row["family"], row["gene_id"],
                )
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()

    rows = []
    genes = set()
    for path in args.inputs:
        for row in load(path):
            gene = row[8]
            if gene in genes:
                raise SystemExit(f"ERROR: duplicate gene_id across taxonomy tables: {gene}")
            genes.add(gene)
            rows.append(row)

    rows.sort(key=lambda row: (row[0], row[1], row[2], row[8]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as out:
        out.write("\t".join(HEADER) + "\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")
    temporary.replace(args.output)
    print(f"Merged {len(rows)} genes from {len(args.inputs)} tables: {args.output}")


if __name__ == "__main__":
    main()
