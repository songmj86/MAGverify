#!/usr/bin/env bash
set -Eeuo pipefail

INTEGRATED_SCRIPT_ID="05_phantom"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"

run_hifi_phantom_phase() (
SCRIPT_ID="05_phantom"
: "${OUT:?}" "${ASSEMBLY_FA:?}" "${READS:?}" "${ASSEMBLER_MAPPING_MANIFEST:?}"
: "${ZERO_MIN_BP:?}" "${PHANTOM_ZERO_BED:?}" "${PHANTOM_ZERO_TSV:?}"
: "${PHANTOM_ZERO_FA:?}" "${PHANTOM_BLAST_DB:?}" "${PHANTOM_BLAST:?}"
: "${PHANTOM_FINAL:?}" "${PHANTOM_BLAST_THREADS:?}"
: "${PHANTOM_KMER_K:?}" "${PHANTOM_MERYL_READ_DB:?}"
: "${PHANTOM_MERYL_READ_CACHE:?}" "${PHANTOM_KMER_FINAL:?}"
: "${PHANTOM_KMER_CACHE_META:?}"
: "${MAG_CONTIG_MAP:?}" "${PHANTOM_MAG_SITES:?}" "${PHANTOM_MAG_SUMMARY:?}"
: "${PHANTOM_MAG_MERYL_FINAL:?}" "${PHANTOM_MAG_MERYL_CATALOGUE:?}"
: "${PHANTOM_MAG_MERYL_CACHE_META:?}" "${PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META:?}"

PHANTOM_MAG_SUMMARIZER="$SCRIPT_DIR/../scripts/phantom_mag.py"
PHANTOM_CACHE_META="${PHANTOM_CACHE_META:-${PHANTOM_DIR}/.cache/phantom.region.tsv.inputs.sha256}"
PHANTOM_DB_CACHE_META="${PHANTOM_BLAST_DB}.inputs.sha256"
PHANTOM_SOURCE_MAPPING_AUDIT="${PHANTOM_CACHE_DIR}/mapping_manifest.tsv"
FORCE_RERUN_PHANTOM="${FORCE_RERUN_PHANTOM:-N}"

migrate_phantom_layout() {
    local pair old new
    mkdir -p "$PHANTOM_CACHE_DIR"
    if [[ -e "${PHANTOM_DIR}/MAG_phantom_meryl.final.tsv" && -e "${PHANTOM_DIR}/phantom.MAG.tsv" ]]; then
        mv -f "${PHANTOM_DIR}/phantom.MAG.tsv" "$PHANTOM_MAG_SUMMARY"
    fi
    for pair in \
        "MAG_phantom_meryl.final.tsv|$PHANTOM_MAG_MERYL_FINAL" \
        "MAG_phantom_meryl.final.tsv.inputs.sha256|$PHANTOM_MAG_MERYL_CACHE_META" \
        "MAG_phantom_meryl.catalogue.final.tsv|$PHANTOM_MAG_MERYL_CATALOGUE" \
        "MAG_phantom_meryl.catalogue.final.tsv.inputs.sha256|$PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META" \
        "phantom_kmer_confirmation.final.tsv|$PHANTOM_KMER_FINAL" \
        "phantom_kmer_confirmation.final.tsv.inputs.sha256|$PHANTOM_KMER_CACHE_META" \
        "phantom.MAG_region.tsv|$PHANTOM_MAG_SITES" \
        "zero_regions.bed|$PHANTOM_ZERO_BED" "zero_regions.fa|$PHANTOM_ZERO_FA" \
        "zero_regions.tsv|$PHANTOM_ZERO_TSV" "zero_regions_vs_raw_reads.tsv|$PHANTOM_BLAST" \
        "raw_reads_blastdb|$(dirname "$PHANTOM_BLAST_DB")" \
        "meryl|$(dirname "$PHANTOM_MERYL_READ_DB")" \
        "phantom.mapping_manifest.tsv|$PHANTOM_SOURCE_MAPPING_AUDIT"; do
        old="${PHANTOM_DIR}/${pair%%|*}"
        new="${pair#*|}"
        [[ -e "$old" && "$old" != "$new" ]] || continue
        if [[ ! -e "$new" ]]; then
            mkdir -p "$(dirname "$new")"
            mv "$old" "$new"
        elif [[ -f "$old" ]]; then
            rm -f "$old"
        fi
    done
}

write_phantom_dataset() {
    python3 - "$PHANTOM_KMER_FINAL" "$PHANTOM_MAG_MERYL_CATALOGUE" "$PHANTOM_DATASET" <<'PY_DATASET'
import csv
import os
import sys
kmer_path, catalogue_path, output = sys.argv[1:]
fields = ['scope', 'k', 'n_MAGs_total', 'n_MAGs_with_phantom', 'n_regions', 'phantom_bp',
          'candidate_distinct_kmers', 'absent_from_reads_distinct_kmers', 'absent_kmer_fraction',
          'interpretation']
rows = []
with open(kmer_path, encoding='utf-8', newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        rows.append({'scope': 'BLAST_negative_regions_all_contigs', 'n_MAGs_total': 'NA',
                     'n_MAGs_with_phantom': 'NA', 'n_regions': row['n_blast_negative_regions'],
                     'phantom_bp': 'NA', **{k: row[k] for k in fields if k in row}})
if os.path.isfile(catalogue_path):
    with open(catalogue_path, encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            rows.append({'scope': 'phantom_regions_in_MAGs', 'n_regions': row['n_phantom_regions'],
                         **{k: row[k] for k in fields if k in row}})
temp = f'{output}.tmp.{os.getpid()}'
with open(temp, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
os.replace(temp, output)
PY_DATASET
}

migrate_phantom_layout
mkdir -p "$OUT/00_log" "${PHANTOM_DIR}" \
    "$(dirname "$PHANTOM_BLAST_DB")" "$(dirname "$PHANTOM_FINAL")" \
    "$(dirname "$PHANTOM_MAG_SITES")" "$(dirname "$PHANTOM_MAG_SUMMARY")" \
    "$(dirname "$PHANTOM_MAG_MERYL_FINAL")"
LOG="$OUT/00_log/05_phantom.log"
exec > >(tee -a "$LOG") 2>&1

for tool in samtools bedtools blastn makeblastdb blastdbcmd meryl python3 awk; do
    command -v "$tool" >/dev/null 2>&1 \
        || { echo "ERROR: missing tool: $tool"; exit 1; }
done
[[ -s "$ASSEMBLY_FA" ]] || { echo "ERROR: missing ASSEMBLY_FA: $ASSEMBLY_FA"; exit 1; }
[[ -s "$READS" ]] || { echo "ERROR: missing READS: $READS"; exit 1; }
[[ -s "$MAG_CONTIG_MAP" ]] || { echo "ERROR: missing MAG_CONTIG_MAP: $MAG_CONTIG_MAP"; exit 1; }
[[ -f "$PHANTOM_MAG_SUMMARIZER" && -r "$PHANTOM_MAG_SUMMARIZER" ]] || {
    echo "ERROR: MAG phantom summarizer is missing or unreadable: $PHANTOM_MAG_SUMMARIZER"
    exit 1
}
[[ "$FORCE_RERUN_PHANTOM" == "Y" || "$FORCE_RERUN_PHANTOM" == "N" ]] \
    || { echo "ERROR: FORCE_RERUN_PHANTOM must be Y or N"; exit 1; }

python3 - "$ZERO_MIN_BP" "$PHANTOM_BLAST_THREADS" "$PHANTOM_KMER_K" <<'PY_VALIDATE'
import sys
zero, threads, kmer_k = map(int, sys.argv[1:])
if zero < 501:
    raise SystemExit('ERROR: ZERO_MIN_BP must be >=501 for the paper >500-bp rule')
if threads < 1:
    raise SystemExit('ERROR: PHANTOM_BLAST_THREADS must be >=1')
if kmer_k != 21:
    raise SystemExit('ERROR: paper confirmation requires PHANTOM_KMER_K=21')
PY_VALIDATE

samtools faidx "$ASSEMBLY_FA"
[[ -s "$ASSEMBLER_MAPPING_MANIFEST" ]] || {
    echo "ERROR: missing mapping manifest: $ASSEMBLER_MAPPING_MANIFEST"
    echo "       rerun 01_mapping.sh with per-assembler HiFi mapping enabled"
    exit 1
}

while IFS=$'\t' read -r assembler source_fasta hifi_all; do
    [[ -s "$source_fasta" ]] || { echo "ERROR: assembly missing: $source_fasta"; exit 1; }
    samtools faidx "$source_fasta"
    [[ -s "$hifi_all" ]] || { echo "ERROR: HiFi_all_BAM missing for $assembler: $hifi_all"; exit 1; }
    samtools quickcheck -v "$hifi_all" || { echo "ERROR: samtools quickcheck failed: $hifi_all"; exit 1; }
    qc_bam_matches_fai "$hifi_all" "${source_fasta}.fai" || {
        echo "ERROR: HiFi BAM/assembly mismatch: $assembler"
        exit 1
    }
done < <(qc_manifest_columns "$ASSEMBLER_MAPPING_MANIFEST" Assembler Assembly_FASTA HiFi_all_BAM)

python3 - "$ASSEMBLER_MAPPING_MANIFEST" "${ASSEMBLY_FA}.fai" <<'PY_SOURCE_UNION'
import csv
import sys
manifest, assemblies_fai = sys.argv[1:]

def load_fai(path):
    result = {}
    with open(path, encoding='utf-8') as handle:
        for raw in handle:
            fields = raw.rstrip('\n\r').split('\t')
            if len(fields) < 2:
                raise SystemExit(f'ERROR: malformed FAI row: {path}')
            name, length = fields[0], int(fields[1])
            if name in result:
                raise SystemExit(f'ERROR: duplicate contig ID in FAI: {name}')
            result[name] = length
    return result

combined = load_fai(assemblies_fai)
union = {}
with open(manifest, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'Assembler', 'Assembly_FASTA', 'HiFi_all_BAM'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'ERROR: assembler mapping manifest lacks columns: {sorted(missing)}')
    rows = list(reader)
if not rows:
    raise SystemExit('ERROR: assembler mapping manifest contains no rows')
for row in rows:
    assembler = row['Assembler'].strip()
    source = load_fai(row['Assembly_FASTA'].strip() + '.fai')
    for contig, length in source.items():
        if contig in union:
            raise SystemExit(f'ERROR: contig ID occurs in multiple assemblers: {contig}')
        union[contig] = (assembler, length)
if set(union) != set(combined):
    missing = sorted(set(combined) - set(union))[:20]
    extra = sorted(set(union) - set(combined))[:20]
    raise SystemExit(
        'ERROR: assembler contig union differs from concatenated assembly; '
        f'missing_from_sources={missing}, extra_in_sources={extra}'
    )
for contig, (_, length) in union.items():
    if combined[contig] != length:
        raise SystemExit(
            f'ERROR: contig length mismatch for {contig}: assembler={length}, combined={combined[contig]}'
        )
PY_SOURCE_UNION

cp -f -- "$ASSEMBLER_MAPPING_MANIFEST" "${PHANTOM_SOURCE_MAPPING_AUDIT}.tmp.$$"
mv -f -- "${PHANTOM_SOURCE_MAPPING_AUDIT}.tmp.$$" "$PHANTOM_SOURCE_MAPPING_AUDIT"

SOURCE_MAPPING_FP="$(python3 - "$ASSEMBLER_MAPPING_MANIFEST" "${ASSEMBLY_FA}.fai" <<'PY_SOURCE_FP'
import csv
import hashlib
import json
import os
import sys
manifest, assemblies_fai = sys.argv[1:]

def file_stat(path):
    st = os.stat(path)
    return [os.path.realpath(path), st.st_size, st.st_mtime_ns]

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def manifest_digest(path):
    with open(path, 'rb') as handle:
        lines = handle.read().split(b'\n')
    return hashlib.sha256(b'\n'.join(b'\t'.join(line.split(b'\t')[:7]) for line in lines)).hexdigest()

rows = []
with open(manifest, encoding='utf-8-sig', newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        bam = row['HiFi_all_BAM'].strip()
        meta = os.path.join(os.path.dirname(bam), '.cache', os.path.basename(bam) + '.inputs.sha256')
        rows.append({
            'assembler': row['Assembler'].strip(),
            'assembly': file_stat(row['Assembly_FASTA'].strip()),
            'assembly_fai_sha256': digest(row['Assembly_FASTA'].strip() + '.fai'),
            'hifi_all_bam': file_stat(bam),
            'hifi_all_bam_cache': (open(meta, encoding='utf-8').read().strip() if os.path.isfile(meta) else 'MISSING'),
        })
payload = {
    'design': 'source_assembler_specific_all_mapped_zero_depth_v1',
    'assemblies_fai_sha256': digest(assemblies_fai),
    'manifest_sha256': manifest_digest(manifest),
    'sources': rows,
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_SOURCE_FP
)"
SAMTOOLS_VERSION="$(samtools --version | awk 'NR==1{print}')"

write_mag_phantom_outputs() {
    [[ -s "$PHANTOM_FINAL" ]] || {
        echo "ERROR: phantom table required for MAG outputs is missing: $PHANTOM_FINAL"
        exit 1
    }
    echo
    echo "Writing MAG-associated phantom outputs"
    python3 "$PHANTOM_MAG_SUMMARIZER" \
        --phantom "$PHANTOM_FINAL" \
        --mag-contig-map "$MAG_CONTIG_MAP" \
        --fai "${ASSEMBLY_FA}.fai" \
        --site-output "$PHANTOM_MAG_SITES" \
        --summary-output "$PHANTOM_MAG_SUMMARY" \
        --contig-output "$PHANTOM_CONTIG"
}

DB_FP="$(python3 - "$READS" <<'PY_DB_FP'
import hashlib, json, os, subprocess, sys
reads = sys.argv[1]
st = os.stat(reads)
payload = {
    'reads': [os.path.realpath(reads), st.st_size, st.st_mtime_ns],
    'makeblastdb': subprocess.check_output(
        ['makeblastdb', '-version'], text=True
    ).splitlines()[0],
    'converter': 'strict_four_line_fastq_to_unique_fasta_v1',
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_DB_FP
)"
BLAST_VERSION="$(blastn -version | awk 'NR==1{print}')"
MERYL_VERSION="$(meryl --version 2>&1 | awk 'NR==1{print}')"
MERYL_READ_FP="$(qc_hash_text "$DB_FP" "$MERYL_VERSION" \
    "k=$PHANTOM_KMER_K" "raw_read_distinct_kmers_v1")"
PHANTOM_FP="$(qc_hash_text "$SOURCE_MAPPING_FP" "$DB_FP" \
    "$ZERO_MIN_BP" "$BLAST_VERSION" "$MERYL_VERSION" "$SAMTOOLS_VERSION" \
    "k=$PHANTOM_KMER_K" \
    "source_assembler_specific_samtools_depth_aa_include_secondary_v1" \
    "dust_no_max_target_1_query_union_meryl_difference_v3")"
MAG_MAP_SHA="$(qc_file_sha256 "$MAG_CONTIG_MAP")"
PHANTOM_MAG_MERYL_FP="$(qc_hash_text "$PHANTOM_FP" "$MAG_MAP_SHA" \
    "mag_level_meryl_union_per_mag_and_catalogue_v1")"

echo "[05_phantom] Started: $(date)"
echo "SCRIPT_ID: $SCRIPT_ID"
echo "Zero-coverage reference: each contig's original assembly"
echo "Zero-coverage population: ALL mapped alignments (primary+secondary+supplementary)"
echo "Minimum zero-region length: $ZERO_MIN_BP bp"
echo "Raw-read confirmation: BLASTN -dust no -max_target_seqs 1"
echo "K-mer confirmation: Meryl difference, k=$PHANTOM_KMER_K"

TMP_BED="${PHANTOM_ZERO_BED}.tmp.$$"
TMP_TSV="${PHANTOM_ZERO_TSV}.tmp.$$"
TMP_FA="${PHANTOM_ZERO_FA}.tmp.$$"
TMP_BLAST="${PHANTOM_BLAST}.tmp.$$"
TMP_FINAL="${PHANTOM_FINAL}.tmp.$$"
TMP_NOHIT_FA="${PHANTOM_CACHE_DIR}/no_blast_hit_regions.fa.tmp.$$"
TMP_KMER_FINAL="${PHANTOM_KMER_FINAL}.tmp.$$"
TMP_CANDIDATE_DB="$(dirname "$PHANTOM_MERYL_READ_DB")/candidates.k${PHANTOM_KMER_K}.tmp.$$.meryl"
TMP_DIFFERENCE_DB="$(dirname "$PHANTOM_MERYL_READ_DB")/absent.k${PHANTOM_KMER_K}.tmp.$$.meryl"
TMP_READ_DB="${PHANTOM_MERYL_READ_DB}.tmp.$$"
TMP_MAG_MERYL_DIR="${PHANTOM_CACHE_DIR}/mag_meryl.tmp.$$"

remove_meryl_db() {
    local path="$1"
    [[ -e "$path" ]] || return 0
    python3 - "$OUT" "$path" <<'PY_REMOVE_DB'
import os
import shutil
import sys
root = os.path.realpath(sys.argv[1])
path = os.path.realpath(sys.argv[2])
if path == root or not path.startswith(root + os.sep):
    raise SystemExit(f'ERROR: refusing to remove Meryl DB outside OUT: {path}')
if os.path.isdir(path):
    shutil.rmtree(path)
else:
    os.unlink(path)
PY_REMOVE_DB
}

ensure_read_meryl_db() {
    mkdir -p "$(dirname "$PHANTOM_MERYL_READ_DB")"
    if ! qc_cache_matches "$PHANTOM_MERYL_READ_CACHE" "$MERYL_READ_FP" \
        || ! meryl statistics "$PHANTOM_MERYL_READ_DB" >/dev/null 2>&1; then
        echo "Building cached raw-read Meryl database (k=$PHANTOM_KMER_K)."
        remove_meryl_db "$PHANTOM_MERYL_READ_DB"
        remove_meryl_db "$TMP_READ_DB"
        meryl count "k=$PHANTOM_KMER_K" "$READS" output "$TMP_READ_DB"
        meryl statistics "$TMP_READ_DB" >/dev/null
        mv "$TMP_READ_DB" "$PHANTOM_MERYL_READ_DB"
        qc_write_cache_atomic "$PHANTOM_MERYL_READ_CACHE" "$MERYL_READ_FP"
    else
        echo "Reusing fingerprinted raw-read Meryl database."
    fi
}

write_mag_meryl_outputs() {
    [[ -s "$PHANTOM_MAG_SUMMARY" ]] || {
        echo "ERROR: missing MAG phantom summary: $PHANTOM_MAG_SUMMARY"
        exit 1
    }
    [[ -s "$PHANTOM_MAG_SITES" ]] || {
        [[ -f "$PHANTOM_MAG_SITES" ]] || {
            echo "ERROR: missing MAG phantom site table: $PHANTOM_MAG_SITES"
            exit 1
        }
    }

    rm -rf "$TMP_MAG_MERYL_DIR"
    mkdir -p "$TMP_MAG_MERYL_DIR"
    local manifest="$TMP_MAG_MERYL_DIR/manifest.tsv"
    local stats="$TMP_MAG_MERYL_DIR/per_mag.stats.tsv"
    local catalogue_bed="$TMP_MAG_MERYL_DIR/catalogue.bed"

    python3 - "$PHANTOM_MAG_SUMMARY" "$PHANTOM_MAG_SITES" \
        "$TMP_MAG_MERYL_DIR" "$manifest" "$catalogue_bed" <<'PY_MAG_PREP'
import csv
import os
import sys
from collections import defaultdict

summary_path, sites_path, workdir, manifest_path, catalogue_bed = sys.argv[1:]
summary_required = [
    'MAG_ID', 'n_contigs', 'total_length_bp', 'n_phantom_regions',
    'phantom_bp', 'phantom_fraction', 'phantom_contigs', 'phantom_region_ids'
]
site_required = ['MAG_ID', 'region_id', 'contig', 'start_0based', 'end_0based', 'length_bp']

with open(summary_path, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    missing = [c for c in summary_required if c not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit(f'ERROR: MAG phantom summary lacks columns: {missing}')
    summaries = []
    seen = set()
    for line_no, row in enumerate(reader, start=2):
        mag = (row.get('MAG_ID') or '').strip()
        if not mag or mag in seen:
            raise SystemExit(f'ERROR: blank or duplicate MAG_ID at summary line {line_no}: {mag}')
        seen.add(mag)
        try:
            n_regions = int(row['n_phantom_regions'])
            phantom_bp = int(row['phantom_bp'])
        except ValueError as exc:
            raise SystemExit(f'ERROR: invalid MAG phantom counts at summary line {line_no}') from exc
        if n_regions < 0 or phantom_bp < 0:
            raise SystemExit(f'ERROR: negative MAG phantom count at summary line {line_no}')
        summaries.append(row)

sites = defaultdict(list)
with open(sites_path, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    missing = [c for c in site_required if c not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit(f'ERROR: MAG phantom site table lacks columns: {missing}')
    for line_no, row in enumerate(reader, start=2):
        mag = (row.get('MAG_ID') or '').strip()
        if mag not in seen:
            raise SystemExit(f'ERROR: site MAG_ID absent from summary at line {line_no}: {mag}')
        try:
            start = int(row['start_0based'])
            end = int(row['end_0based'])
            length = int(row['length_bp'])
        except ValueError as exc:
            raise SystemExit(f'ERROR: invalid site coordinates at line {line_no}') from exc
        if start < 0 or end <= start or end - start != length:
            raise SystemExit(f'ERROR: inconsistent site coordinates at line {line_no}')
        sites[mag].append((row['contig'], start, end, row['region_id'], length))

for row in summaries:
    mag = row['MAG_ID']
    expected = int(row['n_phantom_regions'])
    if len(sites.get(mag, [])) != expected:
        raise SystemExit(
            f'ERROR: {mag} has {expected} phantom regions in summary but '
            f'{len(sites.get(mag, []))} in site table'
        )

with open(manifest_path, 'w', encoding='utf-8', newline='') as out:
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerow([
        'token', 'MAG_ID', 'n_contigs', 'total_length_bp', 'n_phantom_regions',
        'phantom_bp', 'phantom_fraction', 'phantom_contigs', 'phantom_region_ids'
    ])
    for idx, row in enumerate(summaries, start=1):
        token = f'M{idx:06d}'
        writer.writerow([token] + [row[c] for c in summary_required])
        if int(row['n_phantom_regions']) > 0:
            bed_path = os.path.join(workdir, token + '.bed')
            with open(bed_path, 'w', encoding='utf-8') as bed:
                for contig, start, end, region_id, _ in sites[row['MAG_ID']]:
                    bed.write(f'{contig}\t{start}\t{end}\t{region_id}\n')

with open(catalogue_bed, 'w', encoding='utf-8') as bed:
    for row in summaries:
        for contig, start, end, region_id, _ in sites.get(row['MAG_ID'], []):
            bed.write(f'{contig}\t{start}\t{end}\t{row["MAG_ID"]}|{region_id}\n')
PY_MAG_PREP

    printf 'token\tcandidate_distinct_kmers\tabsent_from_reads_distinct_kmers\tabsent_kmer_fraction\n' > "$stats"

    local token n_regions candidate_fa candidate_db absent_db total absent fraction
    while IFS=$'\t' read -r token n_regions; do
        (( n_regions > 0 )) || continue
        candidate_fa="$TMP_MAG_MERYL_DIR/${token}.fa"
        candidate_db="$TMP_MAG_MERYL_DIR/${token}.candidate.meryl"
        absent_db="$TMP_MAG_MERYL_DIR/${token}.absent.meryl"
        bedtools getfasta -fi "$ASSEMBLY_FA" \
            -bed "$TMP_MAG_MERYL_DIR/${token}.bed" -nameOnly -fo "$candidate_fa"
        meryl count "k=$PHANTOM_KMER_K" "$candidate_fa" output "$candidate_db"
        meryl difference "$candidate_db" "$PHANTOM_MERYL_READ_DB" output "$absent_db"
        total="$(meryl threads=1 print "$candidate_db" | awk 'END{print NR+0}')"
        absent="$(meryl threads=1 print "$absent_db" | awk 'END{print NR+0}')"
        [[ "$total" =~ ^[0-9]+$ && "$absent" =~ ^[0-9]+$ ]] || {
            echo "ERROR: invalid Meryl counts for $token"
            exit 1
        }
        (( absent <= total )) || {
            echo "ERROR: absent k-mers exceed candidate k-mers for $token"
            exit 1
        }
        fraction="$(awk -v a="$absent" -v t="$total" \
            'BEGIN{if(t>0) printf "%.6f",a/t; else printf "NA"}')"
        printf '%s\t%s\t%s\t%s\n' "$token" "$total" "$absent" "$fraction" >> "$stats"
        remove_meryl_db "$candidate_db"
        remove_meryl_db "$absent_db"
        rm -f "$candidate_fa"
    done < <(awk -F'\t' 'NR>1{print $1"\t"$5}' "$manifest")

    python3 - "$manifest" "$stats" "$PHANTOM_KMER_K" \
        "$PHANTOM_MAG_MERYL_FINAL" <<'PY_MAG_MERGE'
import csv
import os
import sys
import tempfile

manifest_path, stats_path, k, output_path = sys.argv[1:]
stats = {}
with open(stats_path, encoding='utf-8', newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        stats[row['token']] = row

header = [
    'MAG_ID', 'n_contigs', 'total_length_bp', 'n_phantom_regions', 'phantom_bp',
    'phantom_fraction', 'phantom_contigs', 'phantom_region_ids', 'k',
    'candidate_distinct_kmers', 'absent_from_reads_distinct_kmers',
    'absent_kmer_fraction', 'meryl_interpretation'
]
outdir = os.path.dirname(os.path.abspath(output_path))
os.makedirs(outdir, exist_ok=True)
fd, tmp = tempfile.mkstemp(prefix=os.path.basename(output_path)+'.', suffix='.tmp', dir=outdir)
os.close(fd)
try:
    with open(manifest_path, encoding='utf-8', newline='') as source, \
         open(tmp, 'w', encoding='utf-8', newline='') as out:
        reader = csv.DictReader(source, delimiter='\t')
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        for row in reader:
            n_regions = int(row['n_phantom_regions'])
            if n_regions == 0:
                total = absent = '0'
                fraction = 'NA'
                interpretation = 'No_phantom_region'
            else:
                stat = stats.get(row['token'])
                if stat is None:
                    raise SystemExit(f'ERROR: missing Meryl result for {row["MAG_ID"]}')
                total = stat['candidate_distinct_kmers']
                absent = stat['absent_from_reads_distinct_kmers']
                fraction = stat['absent_kmer_fraction']
                interpretation = 'Measurement_only_no_binary_threshold'
            writer.writerow([
                row['MAG_ID'], row['n_contigs'], row['total_length_bp'],
                row['n_phantom_regions'], row['phantom_bp'], row['phantom_fraction'],
                row['phantom_contigs'], row['phantom_region_ids'], k,
                total, absent, fraction, interpretation
            ])
    os.replace(tmp, output_path)
except Exception:
    if os.path.exists(tmp):
        os.unlink(tmp)
    raise
PY_MAG_MERGE

    local n_mags_total n_mags_with n_regions_total phantom_bp_total
    n_mags_total="$(awk 'NR>1{n++} END{print n+0}' "$manifest")"
    n_mags_with="$(awk -F'\t' 'NR>1 && $5>0{n++} END{print n+0}' "$manifest")"
    n_regions_total="$(awk -F'\t' 'NR>1{s+=$5} END{print s+0}' "$manifest")"
    phantom_bp_total="$(awk -F'\t' 'NR>1{s+=$6} END{print s+0}' "$manifest")"

    local cat_total=0 cat_absent=0 cat_fraction="NA"
    if (( n_regions_total > 0 )); then
        local cat_fa="$TMP_MAG_MERYL_DIR/catalogue.fa"
        local cat_db="$TMP_MAG_MERYL_DIR/catalogue.candidate.meryl"
        local cat_abs_db="$TMP_MAG_MERYL_DIR/catalogue.absent.meryl"
        bedtools getfasta -fi "$ASSEMBLY_FA" -bed "$catalogue_bed" \
            -nameOnly -fo "$cat_fa"
        meryl count "k=$PHANTOM_KMER_K" "$cat_fa" output "$cat_db"
        meryl difference "$cat_db" "$PHANTOM_MERYL_READ_DB" output "$cat_abs_db"
        cat_total="$(meryl threads=1 print "$cat_db" | awk 'END{print NR+0}')"
        cat_absent="$(meryl threads=1 print "$cat_abs_db" | awk 'END{print NR+0}')"
        (( cat_absent <= cat_total )) || {
            echo "ERROR: catalogue absent k-mers exceed candidate k-mers"
            exit 1
        }
        cat_fraction="$(awk -v a="$cat_absent" -v t="$cat_total" \
            'BEGIN{if(t>0) printf "%.6f",a/t; else printf "NA"}')"
        remove_meryl_db "$cat_db"
        remove_meryl_db "$cat_abs_db"
        rm -f "$cat_fa"
    fi

    {
        printf 'k\tn_MAGs_total\tn_MAGs_with_phantom\tn_phantom_regions\tphantom_bp\tcandidate_distinct_kmers\tabsent_from_reads_distinct_kmers\tabsent_kmer_fraction\tinterpretation\n'
        if (( n_regions_total > 0 )); then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tMeasurement_only_no_binary_threshold\n' \
                "$PHANTOM_KMER_K" "$n_mags_total" "$n_mags_with" "$n_regions_total" \
                "$phantom_bp_total" "$cat_total" "$cat_absent" "$cat_fraction"
        else
            printf '%s\t%s\t0\t0\t0\t0\t0\tNA\tNo_phantom_regions_in_MAG_catalogue\n' \
                "$PHANTOM_KMER_K" "$n_mags_total"
        fi
    } > "$PHANTOM_MAG_MERYL_CATALOGUE"

    echo "MAG-level Meryl rows: $n_mags_total"
    echo "MAGs with measured phantom-region Meryl evidence: $n_mags_with"
}

cleanup() {
    rm -f "$TMP_BED" "$TMP_TSV" "$TMP_FA" "$TMP_BLAST" \
          "$TMP_FINAL" "$TMP_NOHIT_FA" "$TMP_KMER_FINAL"
    remove_meryl_db "$TMP_CANDIDATE_DB" || true
    remove_meryl_db "$TMP_DIFFERENCE_DB" || true
    remove_meryl_db "$TMP_READ_DB" || true
    rm -rf "$TMP_MAG_MERYL_DIR"
}
trap cleanup EXIT

if [[ "$FORCE_RERUN_PHANTOM" != "Y" ]] \
    && qc_result_cache_matches "$PHANTOM_CACHE_META" "$PHANTOM_FP" \
        "$PHANTOM_FINAL" \
    && qc_result_cache_matches "$PHANTOM_KMER_CACHE_META" "$PHANTOM_FP" \
        "$PHANTOM_KMER_FINAL"; then
    echo "Reusing fingerprinted phantom result: $PHANTOM_FINAL"
    write_mag_phantom_outputs
    if qc_result_cache_matches "$PHANTOM_MAG_MERYL_CACHE_META" \
            "$PHANTOM_MAG_MERYL_FP" "$PHANTOM_MAG_MERYL_FINAL" \
       && qc_result_cache_matches "$PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META" \
            "$PHANTOM_MAG_MERYL_FP" "$PHANTOM_MAG_MERYL_CATALOGUE"; then
        echo "Reusing fingerprinted MAG-level Meryl outputs."
    else
        N_MAG_PHANTOM="$(awk -F'\t' 'NR>1 && $4>0{n++} END{print n+0}' \
            "$PHANTOM_MAG_SUMMARY")"
        if (( N_MAG_PHANTOM > 0 )); then
            ensure_read_meryl_db
        fi
        write_mag_meryl_outputs
        qc_write_result_cache_atomic "$PHANTOM_MAG_MERYL_CACHE_META" \
            "$PHANTOM_MAG_MERYL_FP" "$PHANTOM_MAG_MERYL_FINAL"
        qc_write_result_cache_atomic "$PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META" \
            "$PHANTOM_MAG_MERYL_FP" "$PHANTOM_MAG_MERYL_CATALOGUE"
        remove_meryl_db "$PHANTOM_MERYL_READ_DB"
        rm -f "$PHANTOM_MERYL_READ_CACHE"
    fi
    write_phantom_dataset
    echo "Outputs:"
    echo "  $PHANTOM_FINAL"
    echo "  $PHANTOM_CONTIG"
    echo "  $PHANTOM_MAG_MERYL_FINAL"
    echo "  $PHANTOM_DATASET"
    echo "[05_phantom] Finished: $(date)"
    exit 0
fi

rm -f "$PHANTOM_ZERO_BED" "$PHANTOM_ZERO_TSV" "$PHANTOM_ZERO_FA" \
      "$PHANTOM_BLAST" "$PHANTOM_FINAL" "$PHANTOM_CACHE_META" \
      "$PHANTOM_KMER_FINAL" "$PHANTOM_KMER_CACHE_META" \
      "$PHANTOM_MAG_MERYL_FINAL" "$PHANTOM_MAG_MERYL_CATALOGUE" \
      "$PHANTOM_MAG_MERYL_CACHE_META" "$PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META"

: > "$TMP_BED"
zero_depth_one_assembler() {
    local hifi_all="$1" part="$2"
    samtools depth -aa -g SECONDARY "$hifi_all" \
        | awk -v minbp="$ZERO_MIN_BP" '
            BEGIN { OFS="\t"; active=0; n=0 }
            function emit_run(    len,rid) {
                if (!active) return
                len=end0-start0
                if (len >= minbp) {
                    rid=contig"|zero_"++n"|start="start0+1"|end="end0"|len="len
                    print contig,start0,end0,rid,len
                }
                active=0
            }
            {
                c=$1; p=$2; depth=$3
                if (depth == 0) {
                    if (!active || c != contig || p != end0+1) {
                        emit_run()
                        contig=c; start0=p-1; end0=p; active=1
                    } else {
                        end0=p
                    }
                } else if (active) {
                    emit_run()
                }
            }
            END { emit_run() }
        ' > "$part.tmp"
    mv -f "$part.tmp" "$part"
}
mkdir -p "$PHANTOM_DIR"
ZERO_PARTS_DIR="$(mktemp -d "$PHANTOM_DIR/.zero_depth.XXXXXX")"
ZERO_ROWS=()
while IFS=$'\t' read -r assembler source_fasta hifi_all _rest; do
    [[ "$assembler" == "Assembler" || -z "$assembler" ]] && continue
    ZERO_ROWS+=("$assembler"$'\t'"$hifi_all")
done < "$ASSEMBLER_MAPPING_MANIFEST"
ZERO_JOBS="$(qc_parallel_jobs "${#ZERO_ROWS[@]}")"
echo "  Per-assembler zero-depth: ${#ZERO_ROWS[@]} assemblers, $ZERO_JOBS in parallel"
for row in "${ZERO_ROWS[@]}"; do
    IFS=$'\t' read -r assembler hifi_all <<< "$row"
    ( zero_depth_one_assembler "$hifi_all" "$ZERO_PARTS_DIR/$assembler.bed" ) &
    qc_wait_for_slot "$ZERO_JOBS"
done
qc_wait_all
for row in "${ZERO_ROWS[@]}"; do
    IFS=$'\t' read -r assembler hifi_all <<< "$row"
    [[ -e "$ZERO_PARTS_DIR/$assembler.bed" ]] || { echo "ERROR: zero-depth scan failed for $assembler"; exit 1; }
    cat "$ZERO_PARTS_DIR/$assembler.bed" >> "$TMP_BED"
done
rm -rf "$ZERO_PARTS_DIR"

{
    printf 'contig\tstart_0based\tend_0based\tregion_id\tlength_bp\n'
    cat "$TMP_BED"
} > "$TMP_TSV"

N_ZERO="$(awk 'END{print NR+0}' "$TMP_BED")"
echo "Exact zero-coverage regions >=${ZERO_MIN_BP} bp: $N_ZERO"

if (( N_ZERO == 0 )); then
    : > "$TMP_FA"
    : > "$TMP_BLAST"
    : > "$TMP_NOHIT_FA"
    printf 'region_id\tcontig\tstart_0based\tend_0based\tlength_bp\traw_read_blast_hsps\tbest_hit_pident\tbest_hit_aligned_bp\traw_read_blast_query_covered_bp\traw_read_blast_query_coverage\tphantom_sequence_call\n' > "$TMP_FINAL"
else
    if bedtools getfasta -h 2>&1 \
        | awk '/-nameOnly/{ok=1} END{exit(ok?0:1)}'; then
        bedtools getfasta -fi "$ASSEMBLY_FA" -bed "$TMP_BED" \
            -nameOnly -fo "$TMP_FA"
    else
        bedtools getfasta -fi "$ASSEMBLY_FA" -bed "$TMP_BED" \
            -name -fo "$TMP_FA"
    fi

    if ! qc_cache_matches "$PHANTOM_DB_CACHE_META" "$DB_FP" \
        || ! blastdbcmd -db "$PHANTOM_BLAST_DB" -info >/dev/null 2>&1; then
        echo "Building cached BLAST database from original long reads."
        case "$PHANTOM_BLAST_DB" in
            "$OUT"/*) rm -f "${PHANTOM_BLAST_DB}".* ;;
            *) echo "ERROR: BLAST DB path is outside OUT"; exit 1 ;;
        esac
        python3 - "$READS" <<'PY_FASTQ' | makeblastdb \
            -in - -dbtype nucl -out "$PHANTOM_BLAST_DB" -title raw_long_reads
import gzip
import sys

path = sys.argv[1]
opener = gzip.open if path.lower().endswith('.gz') else open
count = 0
with opener(path, 'rt', encoding='utf-8', errors='strict') as handle:
    while True:
        header = handle.readline()
        if not header:
            break
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not sequence or not plus or not quality:
            raise SystemExit('ERROR: truncated FASTQ record in READS')
        if not header.startswith('@') or not plus.startswith('+'):
            raise SystemExit('ERROR: READS must be four-line FASTQ/FASTQ.GZ')
        sequence = sequence.strip()
        quality = quality.rstrip('\n\r')
        if not sequence or len(sequence) != len(quality):
            raise SystemExit('ERROR: FASTQ sequence/quality length mismatch')
        count += 1
        sys.stdout.write(f'>read_{count}\n{sequence}\n')
if count == 0:
    raise SystemExit('ERROR: READS contains no FASTQ records')
PY_FASTQ
        blastdbcmd -db "$PHANTOM_BLAST_DB" -info >/dev/null
        qc_write_cache_atomic "$PHANTOM_DB_CACHE_META" "$DB_FP"
    else
        echo "Reusing fingerprinted raw-read BLAST database."
    fi

    blastn -query "$TMP_FA" -db "$PHANTOM_BLAST_DB" \
        -dust no -max_target_seqs 1 -num_threads "$PHANTOM_BLAST_THREADS" \
        -outfmt '6 qseqid sseqid pident length qlen qstart qend' -out "$TMP_BLAST"

    python3 - "$TMP_TSV" "$TMP_BLAST" "$TMP_FA" "$TMP_FINAL" \
        "$TMP_NOHIT_FA" <<'PY_SUMMARY'
import csv
import sys

regions_path, blast_path, query_fasta, output, nohit_fasta = sys.argv[1:]
hit_queries = set()
with open(blast_path, encoding='utf-8') as handle:
    for raw in handle:
        fields = raw.rstrip('\n').split('\t')
        if not fields or not fields[0]:
            continue
        hit_queries.add(fields[0].split('::', 1)[0])

def fasta_records(path):
    name = None
    chunks = []
    with open(path, encoding='utf-8') as handle:
        for raw in handle:
            if raw.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name = raw[1:].strip().split()[0].split('::', 1)[0]
                chunks = []
            else:
                chunks.append(raw.strip())
        if name is not None:
            yield name, ''.join(chunks)

candidate_sequences = dict(fasta_records(query_fasta))

with open(regions_path, encoding='utf-8') as source, \
     open(output, 'w', newline='', encoding='utf-8') as out:
    reader = csv.DictReader(source, delimiter='\t')
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerow([
        'region_id', 'contig', 'start_0based', 'end_0based', 'length_bp',
        'raw_read_blast_hsps', 'best_hit_pident', 'best_hit_aligned_bp',
        'raw_read_blast_query_covered_bp', 'raw_read_blast_query_coverage',
        'phantom_sequence_call',
    ])
    for row in reader:
        region_id = row['region_id']
        has_blast_hit = region_id in hit_queries
        sequence = candidate_sequences.get(region_id)
        if sequence is None:
            raise SystemExit(f'ERROR: region missing from FASTA: {region_id}')
        has_ambiguous_base = any(base not in 'ACGTacgt' for base in sequence)

        # PHANTOM_FINAL is a candidate-only table. Regions with any raw-read
        # BLAST HSP or any non-ACGT base are omitted from this output.
        if has_blast_hit or has_ambiguous_base:
            continue

        writer.writerow([
            row['region_id'], row['contig'], row['start_0based'],
            row['end_0based'], row['length_bp'], 0,
            'NA', '0', '0', '0.000000',
            'No_raw_read_BLAST_hit_phantom_candidate',
        ])

with open(nohit_fasta, 'w', encoding='utf-8') as out:
    for name, sequence in candidate_sequences.items():
        has_ambiguous_base = any(base not in 'ACGTacgt' for base in sequence)
        if name not in hit_queries and not has_ambiguous_base:
            out.write(f'>{name}\n{sequence}\n')
PY_SUMMARY

    N_CANDIDATES="$(awk 'NR>1{n++} END{print n+0}' "$TMP_FINAL")"
    echo "No-BLAST-hit phantom-sequence candidates written: $N_CANDIDATES"
fi

N_NOHIT="$(awk '/^>/{n++} END{print n+0}' "$TMP_NOHIT_FA")"
printf 'k\tn_blast_negative_regions\tcandidate_distinct_kmers\tabsent_from_reads_distinct_kmers\tabsent_kmer_fraction\tinterpretation\n' > "$TMP_KMER_FINAL"

if (( N_NOHIT == 0 )); then
    printf '%s\t0\t0\t0\tNA\tNo_BLAST_negative_candidate_regions\n' \
        "$PHANTOM_KMER_K" >> "$TMP_KMER_FINAL"
else
    ensure_read_meryl_db

    remove_meryl_db "$TMP_CANDIDATE_DB"
    remove_meryl_db "$TMP_DIFFERENCE_DB"
    meryl count "k=$PHANTOM_KMER_K" "$TMP_NOHIT_FA" \
        output "$TMP_CANDIDATE_DB"
    meryl difference "$TMP_CANDIDATE_DB" "$PHANTOM_MERYL_READ_DB" \
        output "$TMP_DIFFERENCE_DB"

    TOTAL_KMERS="$(meryl threads=1 print "$TMP_CANDIDATE_DB" | awk 'END{print NR+0}')"
    ABSENT_KMERS="$(meryl threads=1 print "$TMP_DIFFERENCE_DB" | awk 'END{print NR+0}')"
    [[ "$TOTAL_KMERS" =~ ^[0-9]+$ && "$ABSENT_KMERS" =~ ^[0-9]+$ ]] || {
        echo "ERROR: could not count distinct k-mers from Meryl databases"
        exit 1
    }
    (( ABSENT_KMERS <= TOTAL_KMERS )) || {
        echo "ERROR: absent distinct k-mers exceed candidate distinct k-mers"
        exit 1
    }
    ABSENT_FRACTION="$(awk -v absent="$ABSENT_KMERS" -v total="$TOTAL_KMERS" \
        'BEGIN{if(total>0) printf "%.6f", absent/total; else printf "NA"}')"
    printf '%s\t%s\t%s\t%s\t%s\tMeasurement_only_no_paper_binary_threshold\n' \
        "$PHANTOM_KMER_K" "$N_NOHIT" "$TOTAL_KMERS" "$ABSENT_KMERS" \
        "$ABSENT_FRACTION" >> "$TMP_KMER_FINAL"
fi

mv -f "$TMP_BED" "$PHANTOM_ZERO_BED"
mv -f "$TMP_TSV" "$PHANTOM_ZERO_TSV"
mv -f "$TMP_FA" "$PHANTOM_ZERO_FA"
mv -f "$TMP_BLAST" "$PHANTOM_BLAST"
mv -f "$TMP_FINAL" "$PHANTOM_FINAL"
mv -f "$TMP_KMER_FINAL" "$PHANTOM_KMER_FINAL"
qc_write_result_cache_atomic "$PHANTOM_CACHE_META" "$PHANTOM_FP" "$PHANTOM_FINAL"
qc_write_result_cache_atomic "$PHANTOM_KMER_CACHE_META" "$PHANTOM_FP" \
    "$PHANTOM_KMER_FINAL"

write_mag_phantom_outputs
write_mag_meryl_outputs
qc_write_result_cache_atomic "$PHANTOM_MAG_MERYL_CACHE_META" \
    "$PHANTOM_MAG_MERYL_FP" "$PHANTOM_MAG_MERYL_FINAL"
qc_write_result_cache_atomic "$PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META" \
    "$PHANTOM_MAG_MERYL_FP" "$PHANTOM_MAG_MERYL_CATALOGUE"

rm -f "$PHANTOM_ZERO_BED" "$PHANTOM_ZERO_TSV" "$PHANTOM_ZERO_FA" \
      "$PHANTOM_BLAST" "$TMP_NOHIT_FA"
remove_meryl_db "$TMP_CANDIDATE_DB"
remove_meryl_db "$TMP_DIFFERENCE_DB"

case "$PHANTOM_BLAST_DB" in
    "$OUT"/*) rm -f "${PHANTOM_BLAST_DB}".* ;;
    *) echo "ERROR: BLAST DB path is outside OUT"; exit 1 ;;
esac
remove_meryl_db "$PHANTOM_MERYL_READ_DB"
rm -f "$PHANTOM_MERYL_READ_CACHE"

write_phantom_dataset
echo "Outputs:"
echo "  $PHANTOM_FINAL"
echo "  $PHANTOM_CONTIG"
echo "  $PHANTOM_MAG_MERYL_FINAL"
echo "  $PHANTOM_DATASET"
echo "[05_phantom] Finished: $(date)"
)

echo "[05_phantom] Integrated workflow started: $(date)"
echo "SCRIPT_ID: $INTEGRATED_SCRIPT_ID"
run_hifi_phantom_phase

echo "[05_phantom] Integrated workflow finished: $(date)"
