#!/usr/bin/env python3
"""Select assembly FASTAs or precomputed genes using MAG-contig ownership."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, Iterator, Tuple, TextIO


ENTITY_COLUMNS = (
    "mag_id", "entity_id", "genome_id", "bin_id", "genome", "bin", "id"
)
CONTIG_COLUMNS = ("contig", "original_contig", "contig_id", "sequence", "seqid")


def choose_column(fieldnames: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {name.strip().lower(): name for name in fieldnames}
    return next((lookup[name] for name in candidates if name in lookup), None)


def load_map(path: Path) -> Dict[str, str]:
    ownership: Dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: empty MAG contig map: {path}")
        entity_col = choose_column(reader.fieldnames, ENTITY_COLUMNS)
        contig_col = choose_column(reader.fieldnames, CONTIG_COLUMNS)
        if entity_col is None or contig_col is None:
            raise SystemExit(
                "ERROR: MAG contig map needs a MAG/entity column and a contig column"
            )
        for line_no, row in enumerate(reader, 2):
            entity = row[entity_col].strip()
            contig = row[contig_col].strip()
            if not entity or not contig:
                raise SystemExit(
                    f"ERROR: blank MAG/contig value in map line {line_no}: {path}"
                )
            previous = ownership.get(contig)
            if previous is not None and previous != entity:
                raise SystemExit(
                    f"ERROR: contig assigned to multiple MAGs: {contig} ({previous}, {entity})"
                )
            ownership[contig] = entity
    if not ownership:
        raise SystemExit(f"ERROR: MAG contig map contains no contigs: {path}")
    return ownership


def fasta_records(path: Path) -> Iterator[Tuple[str, str, list[str]]]:
    header: str | None = None
    name: str | None = None
    sequence: list[str] = []
    observed: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if header is not None and name is not None:
                    yield name, header, sequence
                header = raw.rstrip("\n\r")
                text = header[1:].strip()
                name = text.split(None, 1)[0] if text else ""
                if not name:
                    raise SystemExit(f"ERROR: blank FASTA identifier at line {line_no}")
                if name in observed:
                    raise SystemExit(f"ERROR: duplicate FASTA identifier: {name}")
                observed.add(name)
                sequence = []
            else:
                if header is None:
                    if raw.strip():
                        raise SystemExit("ERROR: sequence data found before first FASTA header")
                    continue
                sequence.append(raw.rstrip("\n\r"))
    if header is not None and name is not None:
        yield name, header, sequence


def fasta_main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("include", "exclude"))
    parser.add_argument("assembly_fasta", type=Path)
    parser.add_argument("mag_contig_map", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("stats_json", type=Path)
    args = parser.parse_args(argv)

    ownership = load_map(args.mag_contig_map)
    found_mapped: set[str] = set()
    n_input = 0
    n_output = 0
    bp_input = 0
    bp_output = 0

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.output_fasta.open("w", encoding="utf-8", newline="") as out:
        for name, header, sequence_lines in fasta_records(args.assembly_fasta):
            n_input += 1
            length = sum(len(line.strip()) for line in sequence_lines)
            bp_input += length
            is_mag = name in ownership
            if is_mag:
                found_mapped.add(name)
            keep = is_mag if args.mode == "include" else not is_mag
            if not keep:
                continue
            n_output += 1
            bp_output += length
            out.write(header + "\n")
            for line in sequence_lines:
                out.write(line + "\n")

    missing = sorted(set(ownership) - found_mapped)
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise SystemExit(
            f"ERROR: {len(missing)} mapped contigs are absent from assembly: {preview}{suffix}"
        )

    stats = {
        "mode": args.mode,
        "input_contigs": n_input,
        "input_bp": bp_input,
        "mapped_contigs": len(ownership),
        "output_contigs": n_output,
        "output_bp": bp_output,
    }
    args.stats_json.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(stats, sort_keys=True))




def protein_records(path: Path) -> Iterator[Tuple[str, str, list[str]]]:
    header: str | None = None
    gene_id: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if header is not None and gene_id is not None:
                    yield gene_id, header, sequence
                header = raw.rstrip("\n\r")
                text = header[1:].strip()
                gene_id = text.split(None, 1)[0] if text else ""
                if not gene_id:
                    raise SystemExit(
                        f"ERROR: blank protein identifier in {path}:{line_no}"
                    )
                sequence = []
            else:
                if header is None:
                    if raw.strip():
                        raise SystemExit(
                            f"ERROR: protein sequence before first header in {path}:{line_no}"
                        )
                    continue
                sequence.append(raw.rstrip("\n\r"))
    if header is not None and gene_id is not None:
        yield gene_id, header, sequence


def write_protein(out: TextIO, header: str, sequence: list[str]) -> None:
    out.write(header + "\n")
    for line in sequence:
        out.write(line + "\n")


def genes_main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("include", "exclude"))
    parser.add_argument("assembly_fasta", type=Path)
    parser.add_argument("mag_contig_map", type=Path)
    parser.add_argument("all_genes_faa", type=Path)
    parser.add_argument("all_gene_coordinates", type=Path)
    parser.add_argument("output_genes_faa", type=Path)
    parser.add_argument("output_gene_coordinates", type=Path)
    parser.add_argument("stats_json", type=Path)
    args = parser.parse_args(argv)

    ownership = load_map(args.mag_contig_map)
    assembly_contigs: set[str] = set()
    n_input_contigs = n_output_contigs = 0
    bp_input = bp_output = 0
    for name, _header, sequence_lines in fasta_records(args.assembly_fasta):
        assembly_contigs.add(name)
        n_input_contigs += 1
        length = sum(len(line.strip()) for line in sequence_lines)
        bp_input += length
        is_mag = name in ownership
        keep = is_mag if args.mode == "include" else not is_mag
        if keep:
            n_output_contigs += 1
            bp_output += length

    missing = sorted(set(ownership) - assembly_contigs)
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise SystemExit(
            f"ERROR: {len(missing)} mapped contigs are absent from assembly: "
            f"{preview}{suffix}"
        )

    args.output_genes_faa.parent.mkdir(parents=True, exist_ok=True)
    args.output_gene_coordinates.parent.mkdir(parents=True, exist_ok=True)
    proteins = protein_records(args.all_genes_faa)
    n_input_genes = n_output_genes = 0
    with (
        args.all_gene_coordinates.open(encoding="utf-8") as coordinates,
        args.output_genes_faa.open("w", encoding="utf-8", newline="") as faa_out,
        args.output_gene_coordinates.open(
            "w", encoding="utf-8", newline=""
        ) as coordinate_out,
    ):
        for line_no, raw in enumerate(coordinates, 1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) != 4:
                raise SystemExit(
                    f"ERROR: invalid gene-coordinate row in "
                    f"{args.all_gene_coordinates}:{line_no}"
                )
            gene_id, contig, start_text, end_text = fields
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as error:
                raise SystemExit(
                    f"ERROR: invalid gene coordinates in "
                    f"{args.all_gene_coordinates}:{line_no}"
                ) from error
            if not gene_id or not contig or start < 1 or end < start:
                raise SystemExit(
                    f"ERROR: invalid gene-coordinate values in "
                    f"{args.all_gene_coordinates}:{line_no}"
                )

            try:
                faa_gene_id, header, sequence = next(proteins)
            except StopIteration as error:
                raise SystemExit(
                    "ERROR: gene-coordinate table contains more records than protein FASTA"
                ) from error
            if faa_gene_id != gene_id:
                raise SystemExit(
                    "ERROR: protein FASTA and coordinate order disagree at record "
                    f"{line_no}: {faa_gene_id!r} versus {gene_id!r}"
                )

            n_input_genes += 1
            is_mag = contig in ownership
            keep = is_mag if args.mode == "include" else not is_mag
            if keep:
                write_protein(faa_out, header, sequence)
                coordinate_out.write(raw if raw.endswith("\n") else raw + "\n")
                n_output_genes += 1

        try:
            extra_gene_id, _header, _sequence = next(proteins)
        except StopIteration:
            pass
        else:
            raise SystemExit(
                "ERROR: protein FASTA contains records absent from the coordinate table: "
                f"first extra gene={extra_gene_id}"
            )

    stats = {
        "mode": args.mode,
        "input_contigs": n_input_contigs,
        "input_bp": bp_input,
        "mapped_contigs": len(ownership),
        "output_contigs": n_output_contigs,
        "output_bp": bp_output,
        "input_genes": n_input_genes,
        "output_genes": n_output_genes,
    }
    args.stats_json.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(stats, sort_keys=True))




def main():
    import sys
    if len(sys.argv) < 2 or sys.argv[1] not in ("fasta", "genes"):
        raise SystemExit("Usage: chimera_select_sequences.py {fasta|genes} ARGS...")
    (fasta_main if sys.argv[1] == "fasta" else genes_main)(sys.argv[2:])


if __name__ == "__main__":
    main()
