#!/usr/bin/env python3
"""Produce the compact complete-genome table and validate integrated QC tables."""
import csv
import sys
from pathlib import Path

def complete_genome_qc(source, output):
    required = {
        'MAG_ID', 'Circular', 'Taxonomic_chimera_Y_N', 'Phantom_sequence_Y_N',
        'Excessive_repeat_Y_N', 'False_duplication_Y_N', 'False_closure_Y_N',
        'Premature_circularization_Y_N',
    }

    with open(source, newline='') as handle, open(output, 'w', newline='') as out:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f'ERROR: integrated MAG summary lacks complete-genome QC columns: {sorted(missing)}'
            )

        fields = [
            'ID',
            'Circular',
            'Taxonomic_chimera',
            'Phantom_sequence',
            'Excessive_repeat',
            'False_duplication',
            'False_closure_or_premature_circularization',
        ]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()

        for row in reader:
            values = {
                'Circular': row['Circular'],
                'Taxonomic_chimera_Y_N': row['Taxonomic_chimera_Y_N'],
                'Phantom_sequence_Y_N': row['Phantom_sequence_Y_N'],
                'Excessive_repeat_Y_N': row['Excessive_repeat_Y_N'],
                'False_duplication_Y_N': row['False_duplication_Y_N'],
                'False_closure_Y_N': row['False_closure_Y_N'],
                'Premature_circularization_Y_N': row['Premature_circularization_Y_N'],
            }
            invalid = {k: v for k, v in values.items()
                       if k != 'Circular' and v not in {'Y', 'N'}}
            if invalid:
                raise SystemExit(
                    f'ERROR: invalid complete-genome QC field for {row["MAG_ID"]}: {invalid}'
                )
            circular = row['Circular']
            if circular not in {'Y', 'N'}:
                circular = 'NA'

            false_closure_or_premature = (
                'Y' if row['False_closure_Y_N'] == 'Y'
                or row['Premature_circularization_Y_N'] == 'Y' else 'N'
            )
            writer.writerow({
                'ID': row['MAG_ID'],
                'Circular': circular,
                'Taxonomic_chimera': row['Taxonomic_chimera_Y_N'],
                'Phantom_sequence': row['Phantom_sequence_Y_N'],
                'Excessive_repeat': row['Excessive_repeat_Y_N'],
                'False_duplication': row['False_duplication_Y_N'],
                'False_closure_or_premature_circularization': false_closure_or_premature,
            })


def rows(path, required):
    with open(path, newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(required).issubset(header):
            raise ValueError(f"Invalid TSV header: {path}")
        result = list(reader)
    if not result:
        raise ValueError(f"Empty QC table: {path}")
    return result


def validate(contig_path, mag_path, filtered_path, flags_path, validity_path):
    contigs = rows(contig_path, ('MAG_ID', 'contig', 'length'))
    mags = rows(mag_path, ('MAG_ID', 'n_contigs', 'Any_detected_QC_signal_Y_N'))
    filtered = rows(filtered_path, ('MAG_ID', 'contig', 'length'))
    flags = rows(flags_path, ('MAG_ID', 'n_contigs'))
    validity = rows(validity_path, ('ID',))
    ids = [row['MAG_ID'] for row in mags]
    if any(not mag or mag == 'NA' for mag in ids) or len(set(ids)) != len(ids):
        raise ValueError('Invalid or duplicated MAG_ID in MAG summary')
    counts = {}
    keys = set()
    for row in contigs:
        key = row['contig']
        if not key or key in keys:
            raise ValueError(f'Duplicated or blank contig in contig summary: {key}')
        keys.add(key)
        if int(row['length']) < 1:
            raise ValueError(f'Invalid contig length: {key}')
        mag = row['MAG_ID']
        if mag != 'NA':
            counts[mag] = counts.get(mag, 0) + 1
    mag_ids = set(ids)
    if set(counts) != mag_ids:
        raise ValueError('MAG membership differs between MAG and contig summaries')
    for row in mags:
        mag = row['MAG_ID']
        if int(row['n_contigs']) != counts[mag]:
            raise ValueError(f'Contig count mismatch for {mag}')
    if len(filtered) != sum(v for v in counts.values()):
        raise ValueError('MAG contig view row count differs from MAG membership')
    if {row['contig'] for row in filtered} != {row['contig'] for row in contigs if row['MAG_ID'] != 'NA'}:
        raise ValueError('MAG contig view omits or duplicates a contig')
    if len(flags) != len(mags) or {row['MAG_ID'] for row in flags} != mag_ids:
        raise ValueError('QC flags MAG IDs do not match the MAG summary')
    if len(validity) != len(mags) or {row['ID'] for row in validity} != mag_ids:
        raise ValueError('Complete-genome QC IDs do not match the MAG summary')
    print(f'Validated QC summary: {len(mags)} MAGs, {len(contigs)} contigs')




def main():
    if len(sys.argv) == 4 and sys.argv[1] == 'complete-genome':
        complete_genome_qc(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 7 and sys.argv[1] == 'validate':
        try:
            validate(*(Path(arg) for arg in sys.argv[2:]))
        except (ValueError, OSError, KeyError) as exc:
            raise SystemExit(f'ERROR: {exc}') from None
    else:
        raise SystemExit('Usage: summary_finalize.py complete-genome MAG OUTPUT | validate CONTIG MAG MAG_CONTIG FLAGS VALIDITY')


if __name__ == '__main__':
    main()
