#!/usr/bin/env python3
"""Write shared Pyrodigal genes as an anvi'o external gene calls file.

Only contigs listed in the anvi-script-reformat-fasta name map are written, under
their simplified anvi'o names. Start is 0-based and stop is 1-based (anvi'o rule).
"""
import argparse
import os
import sys


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


def main():
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


if __name__ == '__main__':
    main()
