import csv
import os
import sys
from collections import defaultdict

contig_path, mag_path, dataset_path, contig_high, genome_high = sys.argv[1:]
contig_high = float(contig_high)
genome_high = float(genome_high)


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if not fields:
        raise SystemExit(f'ERROR: missing TSV header: {path}')
    return fields, rows


def write_tsv_atomic(path, fields, rows):
    tmp = f'{path}.tmp.{os.getpid()}'
    with open(tmp, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter='\t',
            lineterminator='\n',
            extrasaction='ignore',
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


contig_fields_in, contig_rows_in = read_tsv(contig_path)
required_contig = {
    'contig', 'MAG_ID', 'MAG_structure', 'length',
    'INTRA_REPEAT', 'intra_repeat_bp', 'intra_repeat_fraction',
    'n_intra_repeat_hsps',
}
missing = required_contig - set(contig_fields_in)
if missing:
    raise SystemExit(
        f'ERROR: contig repeat summary lacks required INTRA columns: {sorted(missing)}'
    )

contig_fields = [
    'contig', 'MAG_ID', 'MAG_structure', 'length',
    'INTRA_REPEAT', 'intra_repeat_bp', 'intra_repeat_fraction',
    'n_intra_repeat_hsps', 'REPEAT_ANALYZABLE', 'EXCESSIVE_REPEAT_CONTIG',
]
contig_rows = []
mag_to_excessive_contigs = defaultdict(int)
for row in contig_rows_in:
    length = int(row['length'])
    intra_bp = int(row['intra_repeat_bp'])
    if not (0 <= intra_bp <= length):
        raise SystemExit(
            f'ERROR: invalid INTRA repeat breadth for {row["contig"]}: '
            f'{intra_bp}/{length}'
        )
    intra_fraction = intra_bp / length if length else 0.0
    excessive = intra_fraction >= contig_high

    out = {
        'contig': row['contig'],
        'MAG_ID': row.get('MAG_ID', 'NA') or 'NA',
        'MAG_structure': row.get('MAG_structure', 'unbinned') or 'unbinned',
        'length': length,
        'INTRA_REPEAT': 'Y' if intra_bp > 0 else 'N',
        'intra_repeat_bp': intra_bp,
        'intra_repeat_fraction': f'{intra_fraction:.6f}',
        'n_intra_repeat_hsps': int(row.get('n_intra_repeat_hsps') or 0),
        'REPEAT_ANALYZABLE': 'Y',
        'EXCESSIVE_REPEAT_CONTIG': 'Y' if excessive else 'N',
    }
    contig_rows.append(out)
    mag = out['MAG_ID']
    if mag not in {'', 'NA'} and excessive:
        mag_to_excessive_contigs[mag] += 1

write_tsv_atomic(contig_path, contig_fields, contig_rows)


mag_fields_in, mag_rows_in = read_tsv(mag_path)
required_mag = {
    'MAG_ID', 'MAG_structure', 'n_contigs', 'total_length_bp',
    'INTRA_REPEAT', 'intra_repeat_bp',
    'INTER_REPEAT', 'inter_repeat_bp',
    'TOTAL_REPEAT', 'total_repeat_bp',
}
missing = required_mag - set(mag_fields_in)
if missing:
    raise SystemExit(
        f'ERROR: MAG repeat summary lacks required INTRA/INTER/TOTAL columns: {sorted(missing)}'
    )

mag_fields = [
    'MAG_ID', 'MAG_structure', 'n_contigs', 'total_length_bp',
    'INTRA_REPEAT', 'intra_repeat_bp', 'intra_repeat_fraction',
    'INTER_REPEAT', 'inter_repeat_bp', 'inter_repeat_fraction',
    'TOTAL_REPEAT', 'total_repeat_bp', 'total_repeat_fraction',
    'n_excessive_repeat_contigs', 'HAS_EXCESSIVE_INTRA_REPEAT_CONTIG',
    'REPEAT_ANALYZABLE', 'EXCESSIVE_REPEAT_MAG',
]
mag_rows = []
for row in mag_rows_in:
    mag = row['MAG_ID']
    length = int(row['total_length_bp'])
    intra_bp = int(row['intra_repeat_bp'])
    inter_bp = int(row['inter_repeat_bp'])
    total_bp = int(row['total_repeat_bp'])
    if not (
        0 <= intra_bp <= length
        and 0 <= inter_bp <= length
        and 0 <= total_bp <= length
    ):
        raise SystemExit(f'ERROR: invalid MAG repeat breadth for {mag}')

    intra_fraction = intra_bp / length if length else 0.0
    inter_fraction = inter_bp / length if length else 0.0
    total_fraction = total_bp / length if length else 0.0

    out = {
        'MAG_ID': mag,
        'MAG_structure': row['MAG_structure'],
        'n_contigs': int(row['n_contigs']),
        'total_length_bp': length,
        'INTRA_REPEAT': 'Y' if intra_bp > 0 else 'N',
        'intra_repeat_bp': intra_bp,
        'intra_repeat_fraction': f'{intra_fraction:.6f}',
        'INTER_REPEAT': 'Y' if inter_bp > 0 else 'N',
        'inter_repeat_bp': inter_bp,
        'inter_repeat_fraction': f'{inter_fraction:.6f}',
        'TOTAL_REPEAT': 'Y' if total_bp > 0 else 'N',
        'total_repeat_bp': total_bp,
        'total_repeat_fraction': f'{total_fraction:.6f}',
        # Count/flag use the INTRA-only contig call and do not define the MAG call.
        'n_excessive_repeat_contigs': mag_to_excessive_contigs.get(mag, 0),
        'HAS_EXCESSIVE_INTRA_REPEAT_CONTIG': (
            'Y' if mag_to_excessive_contigs.get(mag, 0) > 0 else 'N'
        ),
        'REPEAT_ANALYZABLE': 'Y',
        # MAG-level call uses the INTRA+INTER union fraction.
        'EXCESSIVE_REPEAT_MAG': 'Y' if total_fraction >= genome_high else 'N',
    }
    mag_rows.append(out)

write_tsv_atomic(mag_path, mag_fields, mag_rows)


# Dataset summary is intentionally contig-level and therefore INTRA-only.
dataset_fields = [
    'scope_id', 'n_contigs', 'total_length_bp',
    'intra_repeat_bp', 'intra_repeat_fraction',
    'n_excessive_repeat_contigs',
]
dataset_rows = []
for scope, members in (
    ('ALL_INPUT_CONTIGS', contig_rows),
    (
        'FINAL_MAG_CONTIGS',
        [r for r in contig_rows if r['MAG_ID'] not in {'', 'NA'}],
    ),
):
    total_length = sum(int(r['length']) for r in members)
    intra_bp = sum(int(r['intra_repeat_bp']) for r in members)
    n_excessive = sum(r['EXCESSIVE_REPEAT_CONTIG'] == 'Y' for r in members)
    dataset_rows.append({
        'scope_id': scope,
        'n_contigs': len(members),
        'total_length_bp': total_length,
        'intra_repeat_bp': intra_bp,
        'intra_repeat_fraction': f'{intra_bp / total_length if total_length else 0.0:.6f}',
        'n_excessive_repeat_contigs': n_excessive,
    })

write_tsv_atomic(dataset_path, dataset_fields, dataset_rows)

print(
    'Final contig call (INTRA-only): '
    f'{sum(r["EXCESSIVE_REPEAT_CONTIG"] == "Y" for r in contig_rows)} excessive contigs'
)
print(
    'Final MAG call (INTRA+INTER union): '
    f'{sum(r["EXCESSIVE_REPEAT_MAG"] == "Y" for r in mag_rows)} excessive MAGs'
)
