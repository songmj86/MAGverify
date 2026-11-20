#!/usr/bin/env python3
"""Combine the repeat HSP and region tables into one table each.

repeats.HSP.tsv     INTRA (same contig) and INTER (same MAG, different contig) BLAST HSPs
repeats.region.tsv  merged INTRA, INTER and TOTAL (INTRA+INTER union) repeat regions

Every input row is kept; a Repeat_type column records which table it came from.
The BED files are not needed: repeats.region.tsv holds the same intervals
(start_1based/end_1based).
"""

import csv
import os
import sys
import tempfile

HSP_COLUMNS = ['Repeat_type', 'MAG_ID', 'qseqid', 'sseqid', 'pident', 'length',
               'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'orientation']


def read_rows(path):
    with open(path, newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def contig_mags(path):
    mags = {}
    for row in read_rows(path):
        contig = (row.get('Contig') or row.get('contig') or '').strip()
        mag = (row.get('MAG_ID') or '').strip()
        if contig and mag:
            mags[contig] = mag
    return mags


def write_atomic(path, columns, rows):
    fd, temp = tempfile.mkstemp(prefix=f'.{os.path.basename(path)}.', dir=os.path.dirname(path) or '.')
    with os.fdopen(fd, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def main():
    if len(sys.argv) != 9:
        raise SystemExit('usage: repeats_consolidate_outputs.py MAG_MAP INTRA_HSP INTER_HSP '
                         'INTRA_REGIONS INTER_REGIONS TOTAL_REGIONS HSP_OUT REGION_OUT')
    mag_map, intra_hsp, inter_hsp, intra_reg, inter_reg, total_reg, hsp_out, region_out = sys.argv[1:]
    mags = contig_mags(mag_map)

    hsps = []
    for row in read_rows(intra_hsp):
        qstart, qend, sstart, send = (int(row[k]) for k in ('qstart', 'qend', 'sstart', 'send'))
        same_strand = (qend >= qstart) == (send >= sstart)
        hsps.append({**row, 'Repeat_type': 'INTRA', 'MAG_ID': mags.get(row['qseqid'], 'NA'),
                     'orientation': 'direct' if same_strand else 'inverted'})
    for row in read_rows(inter_hsp):
        hsps.append({**row, 'Repeat_type': 'INTER'})
    write_atomic(hsp_out, HSP_COLUMNS, hsps)

    regions, columns = [], None
    for path in (intra_reg, inter_reg, total_reg):
        with open(path, newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if columns is None:
                columns = list(reader.fieldnames or [])
            elif list(reader.fieldnames or []) != columns:
                raise SystemExit(f'ERROR: repeat region columns differ in {path}')
            regions.extend(reader)
    write_atomic(region_out, columns, regions)
    print(f'repeats.HSP.tsv: {len(hsps)} HSPs; repeats.region.tsv: {len(regions)} regions')


if __name__ == '__main__':
    main()
