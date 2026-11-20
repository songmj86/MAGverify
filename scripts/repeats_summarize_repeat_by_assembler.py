import csv
import os
import sys
from collections import Counter, defaultdict

assemblies_fai, contig_path, mag_path, dataset_path, *source_fastas = sys.argv[1:]

if not source_fastas:
    raise SystemExit('ERROR: no source assembler FASTAs supplied')


def read_fai(path):
    result = {}
    with open(path, encoding='utf-8') as handle:
        for raw in handle:
            f = raw.rstrip('\n\r').split('\t')
            if len(f) < 2:
                raise SystemExit(f'ERROR: malformed FAI row: {path}')
            name, length = f[0], int(f[1])
            if name in result:
                raise SystemExit(f'ERROR: duplicate contig ID in FAI: {name}')
            result[name] = length
    return result


def assembler_label(path):
    name = os.path.basename(path)
    for suffix in ('.fasta', '.fna', '.fa'):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return os.path.splitext(name)[0]


combined = read_fai(assemblies_fai)
contig_source = {}
assembler_order = []
for fasta in source_fastas:
    assembler = assembler_label(fasta)
    if assembler in assembler_order:
        raise SystemExit(f'ERROR: duplicate assembler label: {assembler}')
    assembler_order.append(assembler)
    source = read_fai(fasta + '.fai')
    for contig, length in source.items():
        if contig in contig_source:
            raise SystemExit(
                f'ERROR: contig ID occurs in multiple source assemblers: {contig}'
            )
        contig_source[contig] = (assembler, length)

if set(contig_source) != set(combined):
    missing = sorted(set(combined) - set(contig_source))[:20]
    extra = sorted(set(contig_source) - set(combined))[:20]
    raise SystemExit(
        'ERROR: source-assembler contig union differs from concatenated assembly; '
        f'missing_from_sources={missing}, extra_in_sources={extra}'
    )
for contig, (_, length) in contig_source.items():
    if combined[contig] != length:
        raise SystemExit(
            f'ERROR: contig length mismatch for {contig}: '
            f'source={length}, combined={combined[contig]}'
        )


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
            handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore'
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


contig_fields, contig_rows = read_tsv(contig_path)
if 'contig' not in contig_fields or 'EXCESSIVE_REPEAT_CONTIG' not in contig_fields:
    raise SystemExit('ERROR: contig repeat summary lacks required columns')
contig_fields = [f for f in contig_fields if f != 'source_assembler']
contig_fields.insert(contig_fields.index('contig') + 1, 'source_assembler')
mag_sources = defaultdict(set)
for row in contig_rows:
    contig = row['contig']
    if contig not in contig_source:
        raise SystemExit(f'ERROR: contig summary ID lacks source assembler: {contig}')
    assembler = contig_source[contig][0]
    row['source_assembler'] = assembler
    mag = row.get('MAG_ID', '').strip()
    if mag and mag != 'NA':
        mag_sources[mag].add(assembler)
write_tsv_atomic(contig_path, contig_fields, contig_rows)

mag_fields, mag_rows = read_tsv(mag_path)
if 'MAG_ID' not in mag_fields:
    raise SystemExit('ERROR: MAG repeat summary lacks MAG_ID')
mag_fields = [f for f in mag_fields if f != 'source_assembler']
mag_fields.insert(mag_fields.index('MAG_ID') + 1, 'source_assembler')
for row in mag_rows:
    mag = row['MAG_ID']
    sources = sorted(mag_sources.get(mag, []), key=lambda x: assembler_order.index(x))
    row['source_assembler'] = '|'.join(sources) if sources else 'NA'
write_tsv_atomic(mag_path, mag_fields, mag_rows)

dataset_fields, dataset_rows = read_tsv(dataset_path)
count_columns = [f'n_excessive_repeat_contigs_{a}' for a in assembler_order]
dataset_fields = [f for f in dataset_fields if f not in count_columns]
dataset_fields.extend(count_columns)

all_rows = contig_rows
mag_contig_rows = [
    row for row in contig_rows
    if row.get('MAG_ID', '').strip() not in {'', 'NA'}
]
rows_by_scope = {
    'ALL_INPUT_CONTIGS': all_rows,
    'FINAL_MAG_CONTIGS': mag_contig_rows,
}
for row in dataset_rows:
    members = rows_by_scope.get(row.get('scope_id', ''))
    if members is None:
        for col in count_columns:
            row[col] = 'NA'
        continue
    counts = Counter(
        r['source_assembler'] for r in members
        if r.get('EXCESSIVE_REPEAT_CONTIG') == 'Y'
    )
    for assembler, col in zip(assembler_order, count_columns):
        row[col] = str(counts.get(assembler, 0))
write_tsv_atomic(dataset_path, dataset_fields, dataset_rows)

all_counts = Counter(
    row['source_assembler'] for row in contig_rows
    if row.get('EXCESSIVE_REPEAT_CONTIG') == 'Y'
)
print('Excessive-repeat contigs by assembler: ' + ', '.join(
    f'{assembler}={all_counts.get(assembler, 0)}' for assembler in assembler_order
))
