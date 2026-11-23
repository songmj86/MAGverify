#!/usr/bin/env python3
"""Repeat tables: evidence summary, normalization, per-assembler columns and consolidated region/HSP tables.

Subcommands:
  evidence           Summarize INTRA/INTER repeat evidence per contig, MAG and dataset
  normalize          Normalize the repeat tables and calls
  assembler          Add assembler columns and the per-assembler dataset rows
  consolidate        Write repeats.region.tsv and repeats.HSP.tsv
"""

import csv
import os
import sys
import tempfile
from collections import Counter, defaultdict


def cmd_evidence(argv):
    """Summarize INTRA/INTER repeat evidence per contig, MAG and dataset"""
    (
        fai_path, map_path,
        intra_bed, inter_bed, total_bed,
        intra_hsp, inter_hsp,
        contig_out, mag_out,
        intra_region_bed, intra_region_tsv,
        inter_region_bed, inter_region_tsv,
        total_region_bed, total_region_tsv,
        dataset_out, contig_high, genome_high,
    ) = argv

    contig_high = float(contig_high)

    genome_high = float(genome_high)

    lengths = {}

    with open(fai_path) as handle:
        for line in handle:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2:
                lengths[f[0]] = int(f[1])

    with open(map_path, newline='') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))

    if rows and len(rows[0]) >= 2 and rows[0][0].strip().lower() in {'mag_id','mag','id'}:
        rows = rows[1:]

    contig_to_mag = {}

    mag_to_contigs = defaultdict(list)

    for row in rows:
        if len(row) < 2:
            continue
        mag, contig = row[0].strip(), row[1].strip()
        if not mag or not contig:
            continue
        contig_to_mag[contig] = mag
        mag_to_contigs[mag].append(contig)

    mag_structure = {
        mag: ('single_contig' if len(members) == 1 else 'multi_contig')
        for mag, members in mag_to_contigs.items()
    }

    def read_bed(path):
        intervals = defaultdict(list)
        bp = defaultdict(int)
        with open(path) as handle:
            for line in handle:
                if not line.strip():
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 3:
                    continue
                c, s, e = f[0], int(f[1]), int(f[2])
                if c not in lengths:
                    raise SystemExit(f'ERROR: repeat BED contig absent from FASTA: {c}')
                if s < 0 or e <= s or e > lengths[c]:
                    raise SystemExit(f'ERROR: invalid interval {c}:{s}-{e}')
                intervals[c].append((s,e))
                bp[c] += e - s
        return intervals, bp

    intra_iv, intra_bp = read_bed(intra_bed)

    inter_iv, inter_bp = read_bed(inter_bed)

    total_iv, total_bp = read_bed(total_bed)

    intra_hsp_count = defaultdict(int)

    with open(intra_hsp, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            if row.get('qseqid'):
                intra_hsp_count[row['qseqid']] += 1

    inter_hsp_count = defaultdict(int)

    with open(inter_hsp, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            q, s = row.get('qseqid',''), row.get('sseqid','')
            if q:
                inter_hsp_count[q] += 1
            if s:
                inter_hsp_count[s] += 1

    summary = {}

    fields = [
        'contig','MAG_ID','MAG_structure','length',
        'INTRA_REPEAT','intra_repeat_bp','intra_repeat_fraction','n_intra_repeat_hsps','EXCESSIVE_INTRA_REPEAT',
        'INTER_REPEAT','inter_repeat_bp','inter_repeat_fraction','n_inter_repeat_hsps',
        'TOTAL_REPEAT','total_repeat_bp','total_repeat_fraction','REPEAT_ANALYZABLE','EXCESSIVE_REPEAT_CONTIG',
    ]

    with open(contig_out, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for c in sorted(lengths):
            L = lengths[c]
            ib, eb, tb = intra_bp.get(c,0), inter_bp.get(c,0), total_bp.get(c,0)
            if not (0 <= ib <= L and 0 <= eb <= L and 0 <= tb <= L):
                raise SystemExit(f'ERROR: invalid repeat breadth for {c}')
            fi = ib/L if L else 0.0
            fe = eb/L if L else 0.0
            ft = tb/L if L else 0.0
            excessive_intra = fi >= contig_high
            excessive_contig = excessive_intra
            mag = contig_to_mag.get(c, 'NA')
            structure = mag_structure.get(mag, 'unbinned')
            row = {
                'contig': c, 'MAG_ID': mag, 'MAG_structure': structure, 'length': L,
                'INTRA_REPEAT': 'Y' if ib > 0 else 'N',
                'intra_repeat_bp': ib, 'intra_repeat_fraction': f'{fi:.6f}',
                'n_intra_repeat_hsps': intra_hsp_count.get(c,0),
                'EXCESSIVE_INTRA_REPEAT': 'Y' if excessive_intra else 'N',
                'INTER_REPEAT': 'Y' if eb > 0 else 'N',
                'inter_repeat_bp': eb, 'inter_repeat_fraction': f'{fe:.6f}',
                'n_inter_repeat_hsps': inter_hsp_count.get(c,0),
                'TOTAL_REPEAT': 'Y' if tb > 0 else 'N',
                'total_repeat_bp': tb, 'total_repeat_fraction': f'{ft:.6f}',
                'REPEAT_ANALYZABLE': 'Y',
                'EXCESSIVE_REPEAT_CONTIG': 'Y' if excessive_contig else 'N',
            }
            writer.writerow(row)
            summary[c] = row

    mag_fields = [
        'MAG_ID','MAG_structure','n_contigs','total_length_bp',
        'INTRA_REPEAT','intra_repeat_bp','intra_repeat_fraction','n_excessive_repeat_contigs','HAS_EXCESSIVE_INTRA_REPEAT_CONTIG',
        'INTER_REPEAT','inter_repeat_bp','inter_repeat_fraction',
        'TOTAL_REPEAT','total_repeat_bp','total_repeat_fraction',
        'REPEAT_ANALYZABLE','EXCESSIVE_REPEAT_MAG',
    ]

    mag_summary = {}

    with open(mag_out, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=mag_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for mag in sorted(mag_to_contigs):
            members = mag_to_contigs[mag]
            L = sum(lengths[c] for c in members)
            ib = sum(int(summary[c]['intra_repeat_bp']) for c in members)
            eb = sum(int(summary[c]['inter_repeat_bp']) for c in members)
            tb = sum(int(summary[c]['total_repeat_bp']) for c in members)
            fi, fe, ft = (ib/L if L else 0.0), (eb/L if L else 0.0), (tb/L if L else 0.0)
            n_excessive_contigs = sum(summary[c]['EXCESSIVE_REPEAT_CONTIG'] == 'Y' for c in members)
            excessive_mag = ft >= genome_high
            row = {
                'MAG_ID': mag, 'MAG_structure': mag_structure[mag], 'n_contigs': len(members), 'total_length_bp': L,
                'INTRA_REPEAT': 'Y' if ib > 0 else 'N', 'intra_repeat_bp': ib, 'intra_repeat_fraction': f'{fi:.6f}',
                'n_excessive_repeat_contigs': n_excessive_contigs,
                'HAS_EXCESSIVE_INTRA_REPEAT_CONTIG': 'Y' if n_excessive_contigs > 0 else 'N',
                'INTER_REPEAT': 'Y' if eb > 0 else 'N', 'inter_repeat_bp': eb, 'inter_repeat_fraction': f'{fe:.6f}',
                'TOTAL_REPEAT': 'Y' if tb > 0 else 'N', 'total_repeat_bp': tb, 'total_repeat_fraction': f'{ft:.6f}',
                'REPEAT_ANALYZABLE': 'Y',
                'EXCESSIVE_REPEAT_MAG': 'Y' if excessive_mag else 'N',
            }
            writer.writerow(row)
            mag_summary[mag] = row

    def write_regions(path_bed, path_tsv, intervals, repeat_type):
        with open(path_bed, 'w') as bed, open(path_tsv, 'w', newline='') as tsv:
            cols = ['contig','MAG_ID','MAG_structure','repeat_type','repeat_region_id','start_1based','end_1based','repeat_region_bp']
            writer = csv.DictWriter(tsv, fieldnames=cols, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for c in sorted(intervals):
                for i, (s,e) in enumerate(intervals[c],1):
                    rid = f'{c}.{repeat_type.lower()}_repeat_region_{i:06d}'
                    bed.write(f'{c}\t{s}\t{e}\t{rid}\t0\t.\n')
                    mag = contig_to_mag.get(c,'NA')
                    writer.writerow({
                        'contig':c, 'MAG_ID':mag, 'MAG_structure':mag_structure.get(mag,'unbinned'),
                        'repeat_type':repeat_type, 'repeat_region_id':rid,
                        'start_1based':s+1, 'end_1based':e, 'repeat_region_bp':e-s,
                    })

    write_regions(intra_region_bed, intra_region_tsv, intra_iv, 'INTRA')

    write_regions(inter_region_bed, inter_region_tsv, inter_iv, 'INTER')

    write_regions(total_region_bed, total_region_tsv, total_iv, 'TOTAL')

    dataset_fields = [
        'scope_id','n_contigs','total_length_bp',
        'intra_repeat_bp','intra_repeat_fraction','n_excessive_repeat_contigs',
    ]

    with open(dataset_out, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=dataset_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for scope, members in [
            ('ALL_INPUT_CONTIGS', sorted(lengths)),
            ('FINAL_MAG_CONTIGS', sorted(contig_to_mag)),
        ]:
            L = sum(lengths[c] for c in members)
            ib = sum(int(summary[c]['intra_repeat_bp']) for c in members)
            writer.writerow({
                'scope_id':scope,
                'n_contigs':len(members), 'total_length_bp':L,
                'intra_repeat_bp':ib, 'intra_repeat_fraction':f'{ib/L if L else 0:.6f}',
                'n_excessive_repeat_contigs':sum(summary[c]['EXCESSIVE_REPEAT_CONTIG']=='Y' for c in members),
            })

    print(f'Contigs summarized: {len(lengths)}')

    print(f'MAGs summarized: {len(mag_to_contigs)}')

    print(f'Contigs with INTRA_REPEAT=Y: {sum(r["INTRA_REPEAT"]=="Y" for r in summary.values())}')

    print(f'Contigs with INTER_REPEAT=Y: {sum(r["INTER_REPEAT"]=="Y" for r in summary.values())}')

    print(f'Contigs with EXCESSIVE_REPEAT_CONTIG=Y: {sum(r["EXCESSIVE_REPEAT_CONTIG"]=="Y" for r in summary.values())}')

    print(f'MAGs with EXCESSIVE_REPEAT_MAG=Y: {sum(r["EXCESSIVE_REPEAT_MAG"]=="Y" for r in mag_summary.values())}')


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


def cmd_normalize(argv):
    """Normalize the repeat tables and calls"""
    contig_path, mag_path, dataset_path, contig_high, genome_high = argv

    contig_high = float(contig_high)

    genome_high = float(genome_high)

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
            'n_excessive_repeat_contigs': mag_to_excessive_contigs.get(mag, 0),
            'HAS_EXCESSIVE_INTRA_REPEAT_CONTIG': (
                'Y' if mag_to_excessive_contigs.get(mag, 0) > 0 else 'N'
            ),
            'REPEAT_ANALYZABLE': 'Y',
            'EXCESSIVE_REPEAT_MAG': 'Y' if total_fraction >= genome_high else 'N',
        }
        mag_rows.append(out)

    write_tsv_atomic(mag_path, mag_fields, mag_rows)

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


def cmd_assembler(argv):
    """Add assembler columns and the per-assembler dataset rows"""
    assemblies_fai, contig_path, mag_path, dataset_path, *source_fastas = argv

    if not source_fastas:
        raise SystemExit('ERROR: no assembler FASTAs supplied')

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
                    f'ERROR: contig ID occurs in multiple assemblers: {contig}'
                )
            contig_source[contig] = (assembler, length)

    if set(contig_source) != set(combined):
        missing = sorted(set(combined) - set(contig_source))[:20]
        extra = sorted(set(contig_source) - set(combined))[:20]
        raise SystemExit(
            'ERROR: assembler contig union differs from concatenated assembly; '
            f'missing_from_sources={missing}, extra_in_sources={extra}'
        )

    for contig, (_, length) in contig_source.items():
        if combined[contig] != length:
            raise SystemExit(
                f'ERROR: contig length mismatch for {contig}: '
                f'assembler={length}, combined={combined[contig]}'
            )

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
            raise SystemExit(f'ERROR: contig summary ID lacks assembler: {contig}')
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


def _consolidate_main():
    if len(sys.argv) != 9:
        raise SystemExit('usage: repeats_tables.py consolidate MAG_MAP INTRA_HSP INTER_HSP '
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


def cmd_consolidate(argv):
    """Write repeats.region.tsv and repeats.HSP.tsv"""
    sys.argv = [sys.argv[0] + " consolidate", *argv]
    return _consolidate_main()


COMMANDS = {
    'evidence': cmd_evidence,
    'normalize': cmd_normalize,
    'assembler': cmd_assembler,
    'consolidate': cmd_consolidate,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
