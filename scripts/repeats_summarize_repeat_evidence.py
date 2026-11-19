import csv
import sys
from collections import defaultdict

(
    fai_path, map_path,
    intra_bed, inter_bed, total_bed,
    intra_hsp, inter_hsp,
    contig_out, mag_out,
    intra_region_bed, intra_region_tsv,
    inter_region_bed, inter_region_tsv,
    total_region_bed, total_region_tsv,
    dataset_out, contig_high, genome_high,
) = sys.argv[1:]
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
        # Contig-level excessive-repeat classification is INTRA-only.
        # INTER/TOTAL values are retained here solely for MAG-level aggregation.
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

# Dataset summary is contig-level and therefore INTRA-only.
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
