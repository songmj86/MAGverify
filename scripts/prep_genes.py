#!/usr/bin/env python3
"""Shared gene predictions: protein-to-coordinate table and anvi'o external gene calls.

Subcommands:
  link               Link predicted proteins to contig coordinates (genes.coordinates.tsv)
  anvio              Write the shared genes as an anvi'o external gene calls file
"""
from __future__ import annotations

import argparse
import os
import sys
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


def _link_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("faa")
    parser.add_argument("gff")
    parser.add_argument("output")
    args = parser.parse_args()
    link(args.faa, args.gff, args.output)


def cmd_link(argv):
    """Link predicted proteins to contig coordinates (genes.coordinates.tsv)"""
    sys.argv = [sys.argv[0] + " link", *argv]
    return _link_main()


def read_gff(path):
    genes = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            partial = 0 if attrs.get('partial', '00') == '00' else 1
            genes.append((fields[0], int(fields[3]), int(fields[4]), fields[6], partial))
    return genes


def read_proteins(path):
    proteins, name, chunks = [], None, []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    proteins.append(''.join(chunks).rstrip('*'))
                name, chunks = line[1:].split(None, 1)[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        proteins.append(''.join(chunks).rstrip('*'))
    return proteins


def read_name_map(path, known):
    """Return original contig -> anvi'o name; column order is detected from GFF contigs."""
    mapping = {}
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) < 2:
                continue
            a, b = fields[0].split()[0] if fields[0] else '', fields[1].split()[0] if fields[1] else ''
            if b in known:
                mapping[b] = a
            elif a in known:
                mapping[a] = b
    return mapping


def _anvio_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gff')
    parser.add_argument('faa')
    parser.add_argument('name_map')
    parser.add_argument('output')
    parser.add_argument('--version', default='NA', help='Pyrodigal version recorded in the table')
    args = parser.parse_args()

    genes = read_gff(args.gff)
    proteins = read_proteins(args.faa)
    if len(genes) != len(proteins):
        sys.exit(f'ERROR: GFF and protein FASTA gene counts differ: {len(genes)} vs {len(proteins)}')
    mapping = read_name_map(args.name_map, {gene[0] for gene in genes})

    temp = f'{args.output}.tmp.{os.getpid()}'
    written = 0
    with open(temp, 'w', encoding='utf-8') as out:
        out.write('gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type'
                  '\tsource\tversion\taa_sequence\n')
        for (contig, start, end, strand, partial), protein in zip(genes, proteins):
            name = mapping.get(contig)
            if name is None:
                continue
            direction = 'f' if strand == '+' else 'r'
            out.write(f'{written}\t{name}\t{start - 1}\t{end}\t{direction}\t{partial}\t1'
                      f'\tpyrodigal\t{args.version}\t{protein}\n')
            written += 1
    if written == 0:
        os.remove(temp)
        sys.exit('ERROR: no genes matched the anvi\'o name map')
    os.replace(temp, args.output)
    print(f'External gene calls: {written} genes -> {args.output}')


def cmd_anvio(argv):
    """Write the shared genes as an anvi'o external gene calls file"""
    sys.argv = [sys.argv[0] + " anvio", *argv]
    return _anvio_main()


COMMANDS = {
    'link': cmd_link,
    'anvio': cmd_anvio,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
