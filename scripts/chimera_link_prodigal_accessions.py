#!/usr/bin/env python3
"""Link Prodigal protein FASTA accessions to exact GFF CDS coordinates."""
from __future__ import annotations

import argparse
from typing import List, Optional, Tuple


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def load_gff(path: str) -> List[Tuple[str, str, int, int, str]]:
    rows: List[Tuple[str, str, int, int, str]] = []
    with open(path, encoding="utf-8") as src:
        for line_no, line in enumerate(src, 1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value
            gff_id = attrs.get("ID", "").strip()
            contig = fields[0].strip()
            if not gff_id or not contig:
                fail(f"missing CDS ID/contig in GFF line {line_no}")
            start, end = int(fields[3]), int(fields[4])
            strand = fields[6].strip()
            if start < 1 or end < start:
                fail(f"invalid coordinates for {gff_id}: {start}-{end}")
            rows.append((gff_id, contig, start, end, strand))
    if not rows:
        fail("no CDS features were recovered from Prodigal GFF")
    return rows


def parse_fasta_header_coordinates(
    header: str,
) -> Optional[Tuple[int, int, str]]:
    parts = [part.strip() for part in header.split(" # ")]
    if len(parts) < 4:
        return None
    try:
        start = int(parts[1])
        end = int(parts[2])
        strand_code = int(parts[3])
    except ValueError:
        return None
    strand = "+" if strand_code == 1 else "-" if strand_code == -1 else "?"
    return start, end, strand


def load_faa(path: str) -> List[Tuple[str, Optional[Tuple[int, int, str]], str]]:
    rows = []
    with open(path, encoding="utf-8") as src:
        for line_no, line in enumerate(src, 1):
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n\r")
            accession = header.split(None, 1)[0].strip()
            if not accession:
                fail(f"blank protein accession in FASTA line {line_no}")
            rows.append((accession, parse_fasta_header_coordinates(header), header))
    if not rows:
        fail("no protein accessions were recovered from Prodigal FASTA")
    return rows


def link(faa: str, gff: str, output: str) -> None:
    faa_rows = load_faa(faa)
    gff_rows = load_gff(gff)
    if len(faa_rows) != len(gff_rows):
        fail(
            "Prodigal FASTA/GFF feature counts differ: "
            f"{len(faa_rows)} versus {len(gff_rows)}"
        )
    seen = set()
    with open(output, "w", encoding="utf-8", newline="") as out:
        for index, ((accession, parsed, header), gff_row) in enumerate(
            zip(faa_rows, gff_rows), 1
        ):
            gff_id, contig, start, end, strand = gff_row
            if accession in seen:
                fail(f"duplicate Prodigal protein accession: {accession}")
            if any(ch in accession + contig for ch in "\t\n\r"):
                fail(f"tab/newline in protein or contig identifier: {accession}")
            if parsed is not None and parsed != (start, end, strand):
                fail(
                    "Prodigal FASTA/GFF order or coordinates disagree at "
                    f"feature {index}: FASTA={accession} {parsed}; "
                    f"GFF={gff_id} {(start, end, strand)}; header={header!r}"
                )
            seen.add(accession)
            out.write(f"{accession}\t{contig}\t{start}\t{end}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("faa")
    parser.add_argument("gff")
    parser.add_argument("output")
    args = parser.parse_args()
    link(args.faa, args.gff, args.output)


if __name__ == "__main__":
    main()
