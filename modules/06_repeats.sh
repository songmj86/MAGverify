#!/usr/bin/env bash
set -euo pipefail

SCRIPT_ID="06_repeats"
REPEAT_EVIDENCE_CACHE_VERSION="repeat_evidence_contig_and_genome"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

: "${OUT:?ERROR: OUT is not set in lib_config.sh}"
: "${ASSEMBLY_FA:?ERROR: ASSEMBLY_FA is not set in lib_config.sh}"
: "${THREADS:?ERROR: THREADS is not set in lib_config.sh}"
: "${REPEAT_SUM:?ERROR: REPEAT_SUM is not set in lib_config.sh}"
: "${MAG_CONTIG_MAP:?ERROR: MAG_CONTIG_MAP is not set in lib_config.sh}"
: "${ASSEMBLER_ASSEMBLY_DIR:?ERROR: ASSEMBLER_ASSEMBLY_DIR is not set in lib_config.sh}"
: "${ASSEMBLER_FASTAS:?ERROR: ASSEMBLER_FASTAS is not set in lib_config.sh}"

REPEAT_FASTA="$ASSEMBLY_FA"
REPEAT_DIR="${REPEAT_DIR}"
LOG="$OUT/00_log/06_repeats.log"

FORCE_RERUN_REPEAT="${FORCE_RERUN_REPEAT:-N}"
REPEAT_CACHE_META="${REPEAT_CACHE_META:-${REPEAT_DIR}/.cache/repeats.contig.tsv.inputs.sha256}"

REPEAT_MIN_ALIGN_PAPER="${REPEAT_MIN_ALIGN_PAPER:-200}"
REPEAT_PIDENT_PAPER="${REPEAT_PIDENT_PAPER:-80}"
REPEAT_CONTIG_HIGH="${REPEAT_CONTIG_HIGH:-${REPEAT_HIGH:-0.70}}"
REPEAT_GENOME_HIGH="${REPEAT_GENOME_HIGH:-${REPEAT_HIGH:-0.70}}"

REPEAT_CHUNK_BP="${REPEAT_CHUNK_BP:-10000000}"
REPEAT_BLAST_THREADS_PER_JOB=1
REPEAT_JOBS="${REPEAT_JOBS:-$(( THREADS / REPEAT_BLAST_THREADS_PER_JOB ))}"
(( REPEAT_JOBS >= 1 )) || REPEAT_JOBS=1
REPEAT_SORT_THREADS="${REPEAT_SORT_THREADS:-$THREADS}"
REPEAT_SORT_MEMORY="${REPEAT_SORT_MEMORY:-50%}"
REPEAT_MIN_CONTIG_LEN="${REPEAT_MIN_CONTIG_LEN:-0}"

(( REPEAT_JOBS > THREADS )) && REPEAT_JOBS="$THREADS"
if [[ -n "${MEMORY_GB:-}" ]]; then
    sort_gb=$(( MEMORY_GB / 2 )); (( sort_gb < 1 )) && sort_gb=1
    REPEAT_SORT_MEMORY="${sort_gb}G"
    memory_jobs=$(( (MEMORY_GB - sort_gb) / 2 )); (( memory_jobs < 1 )) && memory_jobs=1
    (( REPEAT_JOBS > memory_jobs )) && REPEAT_JOBS="$memory_jobs"
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

REPEAT_WORK_DIR="$REPEAT_DIR/.work"
SELF_BLAST="$REPEAT_WORK_DIR/self_blast.paper_same_contig.filtered.tsv"
INTER_BLAST="$REPEAT_WORK_DIR/inter_blast.paper_same_MAG_different_contig.filtered.tsv"
REPEAT_HSP="$REPEAT_DIR/.cache/repeats.HSP.tsv"
REPEAT_REGIONS="$REPEAT_DIR/.cache/repeats.region.tsv"
REPEAT_MAG_MAP="${REPEAT_MAG_MAP:-$MAG_CONTIG_MAP}"
REPEAT_MAG_SUM="${REPEAT_MAG_SUM:-$REPEAT_DIR/.cache/repeats.MAG.tsv}"
REPEAT_DATASET_SUM="${REPEAT_DATASET_SUM:-$REPEAT_DIR/.cache/repeats.dataset.tsv}"

REPEAT_REGION_BED="$REPEAT_WORK_DIR/repeat_regions.final.bed"
REPEAT_REGION_TSV="$REPEAT_WORK_DIR/repeat_regions.final.tsv"
INTRA_REGION_BED="$REPEAT_WORK_DIR/intra_repeat_regions.final.bed"
INTRA_REGION_TSV="$REPEAT_WORK_DIR/intra_repeat_regions.final.tsv"
INTER_REGION_BED="$REPEAT_WORK_DIR/inter_repeat_regions.final.bed"
INTER_REGION_TSV="$REPEAT_WORK_DIR/inter_repeat_regions.final.tsv"

for legacy in self_blast.paper_same_contig.filtered.tsv inter_blast.paper_same_MAG_different_contig.filtered.tsv \
              repeat_regions.final.bed repeat_regions.final.tsv intra_repeat_regions.final.bed \
              intra_repeat_regions.final.tsv inter_repeat_regions.final.bed inter_repeat_regions.final.tsv; do
    if [[ -e "$REPEAT_DIR/$legacy" ]]; then
        mkdir -p "$REPEAT_WORK_DIR"
        mv -f "$REPEAT_DIR/$legacy" "$REPEAT_WORK_DIR/$legacy"
    fi
done
[[ -e "$REPEAT_DIR/dataset_repeat_breadth.final.tsv" && ! -e "$REPEAT_DATASET_SUM" ]] \
    && mv -f "$REPEAT_DIR/dataset_repeat_breadth.final.tsv" "$REPEAT_DATASET_SUM"

mkdir -p "$OUT/00_log" "$REPEAT_DIR" "$REPEAT_WORK_DIR" \
    "$(dirname "$REPEAT_SUM")" "$(dirname "$REPEAT_MAG_SUM")" "$(dirname "$REPEAT_DATASET_SUM")"
exec > >(tee -a "$LOG") 2>&1

echo "[06_repeats] Started: $(date)"
echo "SCRIPT_ID             : $SCRIPT_ID"
echo "REPEAT_FASTA               : $REPEAT_FASTA"
echo "MAG_CONTIG_MAP             : $REPEAT_MAG_MAP"
echo "REPEAT_MIN_ALIGN_PAPER     : >$REPEAT_MIN_ALIGN_PAPER bp"
echo "REPEAT_PIDENT_PAPER        : >=$REPEAT_PIDENT_PAPER%"
echo "REPEAT_CONTIG_HIGH         : >=$REPEAT_CONTIG_HIGH of contig length"
echo "REPEAT_GENOME_HIGH         : >=$REPEAT_GENOME_HIGH of genome length"
echo "INTRA analysis             : same-contig non-diagonal BLASTN"
echo "INTER analysis             : different contigs within the same MAG only"
echo "Contig repeat breadth      : INTRA only"
echo "INTER use                  : MAG total-repeat calculation only"
echo "MAG total repeat breadth   : union(INTRA, INTER); overlap counted once"
echo "Final contig call          : EXCESSIVE_REPEAT_CONTIG=Y when INTRA repeat fraction >= contig threshold"
echo "MAG excessive-contig flag  : HAS_EXCESSIVE_INTRA_REPEAT_CONTIG=Y when >=1 member contig is INTRA-excessive"
echo "Final MAG call             : EXCESSIVE_REPEAT_MAG=Y when MAG-wide total repeat fraction >= genome threshold"
echo "Parallel BLAST jobs        : $REPEAT_JOBS (1 thread each)"
echo "FORCE_RERUN_REPEAT         : $FORCE_RERUN_REPEAT"
echo

for tool in samtools blastn bedtools python3 sort xargs; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: required tool not found in PATH: $tool"
        exit 1
    }
done

[[ -s "$REPEAT_FASTA" ]] || { echo "ERROR: REPEAT_FASTA not found or empty: $REPEAT_FASTA"; exit 1; }
[[ -s "$REPEAT_MAG_MAP" ]] || { echo "ERROR: MAG ownership map not found or empty: $REPEAT_MAG_MAP"; exit 1; }

(( REPEAT_MIN_ALIGN_PAPER >= 1 )) || { echo "ERROR: REPEAT_MIN_ALIGN_PAPER must be >=1"; exit 1; }
(( REPEAT_JOBS >= 1 )) || { echo "ERROR: REPEAT_JOBS must be >=1"; exit 1; }
(( REPEAT_BLAST_THREADS_PER_JOB >= 1 )) || { echo "ERROR: REPEAT_BLAST_THREADS_PER_JOB must be >=1"; exit 1; }
(( REPEAT_SORT_THREADS >= 1 && REPEAT_SORT_THREADS <= THREADS )) || { echo "ERROR: REPEAT_SORT_THREADS must be in [1, THREADS]"; exit 1; }
[[ "$REPEAT_SORT_MEMORY" =~ ^([1-9][0-9]*)(%|[KMGTP])$ ]] || { echo "ERROR: REPEAT_SORT_MEMORY must be a GNU sort size such as 50% or 100G"; exit 1; }
(( REPEAT_CHUNK_BP >= 1 )) || { echo "ERROR: REPEAT_CHUNK_BP must be >=1"; exit 1; }
(( REPEAT_MIN_CONTIG_LEN >= 0 )) || { echo "ERROR: REPEAT_MIN_CONTIG_LEN must be >=0"; exit 1; }
awk -v x="$REPEAT_PIDENT_PAPER" 'BEGIN{exit !(x>=0 && x<=100)}' || { echo "ERROR: REPEAT_PIDENT_PAPER must be 0 to 100"; exit 1; }
awk -v x="$REPEAT_CONTIG_HIGH" 'BEGIN{exit !(x>=0 && x<=1)}' || { echo "ERROR: REPEAT_CONTIG_HIGH must be 0 to 1"; exit 1; }
awk -v x="$REPEAT_GENOME_HIGH" 'BEGIN{exit !(x>=0 && x<=1)}' || { echo "ERROR: REPEAT_GENOME_HIGH must be 0 to 1"; exit 1; }

python3 - "$REPEAT_FASTA" <<'PY_VALIDATE_FASTA'
import sys
path = sys.argv[1]
seen = set()
n = 0
with open(path) as handle:
    for lineno, line in enumerate(handle, 1):
        if not line.startswith('>'):
            continue
        name = line[1:].strip().split()[0] if line[1:].strip() else ''
        if not name:
            raise SystemExit(f'ERROR: blank FASTA identifier at line {lineno}')
        if name in seen:
            raise SystemExit(f'ERROR: duplicate FASTA identifier: {name}')
        seen.add(name)
        n += 1
if n == 0:
    raise SystemExit('ERROR: no FASTA records detected')
print(f'Validated unique FASTA identifiers: {n}')
PY_VALIDATE_FASTA

if [[ ! -s "${REPEAT_FASTA}.fai" || "$REPEAT_FASTA" -nt "${REPEAT_FASTA}.fai" ]]; then
    echo "Creating/updating FASTA index"
    samtools faidx "$REPEAT_FASTA"
else
    echo "Using current FASTA index: ${REPEAT_FASTA}.fai"
fi

python3 - "$REPEAT_MAG_MAP" "${REPEAT_FASTA}.fai" <<'PY_VALIDATE_MAP'
import csv
import sys
from collections import defaultdict

map_path, fai_path = sys.argv[1:]
lengths = {}
with open(fai_path) as handle:
    for line in handle:
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2:
            lengths[f[0]] = int(f[1])

with open(map_path, newline='') as handle:
    rows = list(csv.reader(handle, delimiter='\t'))
if not rows:
    raise SystemExit('ERROR: empty MAG-contig map')
if len(rows[0]) >= 2 and rows[0][0].strip().lower() in {'mag_id','mag','id'}:
    rows = rows[1:]

seen = set()
mags = defaultdict(list)
for row in rows:
    if len(row) < 2:
        continue
    mag, contig = row[0].strip(), row[1].strip()
    if not mag or not contig:
        continue
    if contig in seen:
        raise SystemExit(f'ERROR: contig occurs more than once in MAG map: {contig}')
    if contig not in lengths:
        raise SystemExit(f'ERROR: MAG-map contig absent from assembly FASTA: {contig}')
    seen.add(contig)
    mags[mag].append(contig)
if not mags:
    raise SystemExit('ERROR: no MAG-contig assignments found')
print(f'Validated MAGs: {len(mags)}')
print(f'Validated MAG contigs: {len(seen)}')
print(f'Multi-contig MAGs for INTER analysis: {sum(len(v) > 1 for v in mags.values())}')
PY_VALIDATE_MAP

REPEAT_FINGERPRINT="$(python3 - \
    "$REPEAT_FASTA" "$REPEAT_MAG_MAP" "$REPEAT_MIN_ALIGN_PAPER" \
    "$REPEAT_PIDENT_PAPER" "$REPEAT_CONTIG_HIGH" "$REPEAT_GENOME_HIGH" "$REPEAT_MIN_CONTIG_LEN" \
    "$REPEAT_CHUNK_BP" "$REPEAT_EVIDENCE_CACHE_VERSION" <<'PY_FP'
import hashlib, json, os, subprocess, sys
fasta, mag_map, min_align, pident, contig_high, genome_high, min_contig, chunk_bp, version = sys.argv[1:]
def stat(path):
    s = os.stat(path)
    return [os.path.realpath(path), s.st_size, s.st_mtime_ns]
payload = {
    'fasta': stat(fasta),
    'mag_map': stat(mag_map),
    'min_align': min_align,
    'pident': pident,
    'contig_high': contig_high,
    'genome_high': genome_high,
    'min_contig': min_contig,
    'chunk_bp': chunk_bp,
    'blastn': subprocess.check_output(['blastn','-version'], text=True).splitlines()[0],
    'version': version,
    'method': 'intra_same_contig_plus_inter_sameMAG_total_union',
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_FP
)"

NEED_REPEAT_ANALYSIS=1
CACHED_REPEAT_FINGERPRINT=""
if [[ -s "$REPEAT_CACHE_META" ]]; then
    IFS=$'\t' read -r CACHED_REPEAT_FINGERPRINT _ < "$REPEAT_CACHE_META" || true
fi

if [[ "$FORCE_RERUN_REPEAT" != "Y" \
      && "$CACHED_REPEAT_FINGERPRINT" == "$REPEAT_FINGERPRINT" \
      && -s "$REPEAT_SUM" && -s "$REPEAT_MAG_SUM" \
      && ( ( -s "$REPEAT_HSP" && -s "$REPEAT_REGIONS" ) \
           || ( -s "$SELF_BLAST" && -s "$INTER_BLAST" && -s "$INTRA_REGION_TSV" \
                && -s "$INTER_REGION_TSV" && -s "$REPEAT_REGION_TSV" ) ) ]]; then
    echo "Existing INTRA+INTER repeat evidence matches current inputs; reusing evidence."
    NEED_REPEAT_ANALYSIS=0
else
    echo "Repeat evidence is absent/stale or method changed; rebuilding INTRA+INTER analysis."
fi

merge_bed_file() {
    local input_bed="$1"
    local output_bed="$2"
    if [[ -s "$input_bed" ]]; then
        LC_ALL=C sort --parallel="$REPEAT_SORT_THREADS" -S "$REPEAT_SORT_MEMORY" \
            -k1,1 -k2,2n -k3,3n "$input_bed" \
            | bedtools merge -i - > "$output_bed"
    else
        : > "$output_bed"
    fi
}

if (( NEED_REPEAT_ANALYSIS == 1 )); then
    TMP_ROOT="$(mktemp -d "$REPEAT_DIR/.repeat_run.XXXXXX")"
    trap 'rm -rf "$TMP_ROOT"' EXIT

    TMP_INTRA_HSP="$TMP_ROOT/intra.filtered.tsv"
    TMP_INTER_HSP="$TMP_ROOT/inter.filtered.tsv"
    TMP_INTRA_RAW_BED="$TMP_ROOT/intra.raw.bed"
    TMP_INTER_RAW_BED="$TMP_ROOT/inter.raw.bed"
    TMP_TOTAL_RAW_BED="$TMP_ROOT/total.raw.bed"
    TMP_INTRA_MERGED="$TMP_ROOT/intra.merged.bed"
    TMP_INTER_MERGED="$TMP_ROOT/inter.merged.bed"
    TMP_TOTAL_MERGED="$TMP_ROOT/total.merged.bed"

    printf '%s\n' $'qseqid\tsseqid\tpident\tlength\tqstart\tqend\tsstart\tsend\tqlen\tslen' > "$TMP_INTRA_HSP"
    printf '%s\n' $'MAG_ID\tqseqid\tsseqid\tpident\tlength\tqstart\tqend\tsstart\tsend\tqlen\tslen\torientation' > "$TMP_INTER_HSP"
    : > "$TMP_INTRA_RAW_BED"
    : > "$TMP_INTER_RAW_BED"

    echo "[60A] Grouping contigs for INTRA repeat BLAST"
    INTRA_CHUNK_DIR="$TMP_ROOT/intra_chunks"
    INTRA_DB_DIR="$TMP_ROOT/intra_dbs"
    INTRA_OUT_DIR="$TMP_ROOT/intra_results"
    INTRA_MANIFEST="$TMP_ROOT/intra_chunks.list"
    mkdir -p "$INTRA_CHUNK_DIR" "$INTRA_DB_DIR" "$INTRA_OUT_DIR"

    EFFECTIVE_MIN_CONTIG_LEN="$REPEAT_MIN_CONTIG_LEN"
    if (( EFFECTIVE_MIN_CONTIG_LEN <= REPEAT_MIN_ALIGN_PAPER )); then
        EFFECTIVE_MIN_CONTIG_LEN=$(( REPEAT_MIN_ALIGN_PAPER + 1 ))
    fi

    python3 - "${REPEAT_FASTA}.fai" "$INTRA_CHUNK_DIR" "$INTRA_MANIFEST" \
        "$REPEAT_CHUNK_BP" "$EFFECTIVE_MIN_CONTIG_LEN" "$REPEAT_JOBS" <<'PY_SPLIT'
import sys
from pathlib import Path
fai, out_dir, manifest, max_bp, min_len, jobs = sys.argv[1:]
out_dir = Path(out_dir)
max_bp, min_len, jobs = int(max_bp), int(min_len), int(jobs)
contigs = []
with open(fai) as handle:
    for line in handle:
        name, length = line.split('\t')[:2]
        if int(length) >= min_len:
            contigs.append((int(length), name))
contigs.sort(key=lambda item: (-item[0], item[1]))
total = sum(length for length, _ in contigs)
target = max(1_000_000, min(max_bp, total // max(1, jobs * 4)))
chunks, current, current_bp = [], [], 0
for length, name in contigs:
    if current and current_bp + length > target:
        chunks.append(current)
        current, current_bp = [], 0
    current.append(name)
    current_bp += length
if current:
    chunks.append(current)
with open(manifest, 'w') as handle:
    for idx, names in enumerate(chunks, 1):
        path = out_dir / f'chunk_{idx:06d}.list'
        path.write_text('\n'.join(names) + '\n')
        handle.write(str(path) + '\n')
print(f'Contigs at or above the minimum length: {len(contigs)} in {len(chunks)} chunks (longest first).')
PY_SPLIT

    INTRA_FILTER="$TMP_ROOT/filter_intra.py"
    cat > "$INTRA_FILTER" <<'PY_INTRA_FILTER'
import sys
hsp_out, bed_out, min_align, min_pident = sys.argv[1:]
min_align, min_pident = int(min_align), float(min_pident)
with open(hsp_out, 'w') as hsp, open(bed_out, 'w') as bed:
    for raw in sys.stdin:
        if not raw.strip():
            continue
        f = raw.rstrip('\n').split('\t')
        if len(f) != 10:
            raise SystemExit(f'ERROR: malformed BLAST row ({len(f)} fields)')
        q, s = f[0], f[1]
        if q != s:
            continue
        pident, alen = float(f[2]), int(float(f[3]))
        qs, qe, ss, se = map(int, f[4:8])
        qlen, slen = map(int, f[8:10])
        if qlen != slen:
            raise SystemExit(f'ERROR: same-contig HSP has unequal qlen/slen: {q}')
        if qs == ss and qe == se:
            continue
        if alen <= min_align or pident < min_pident:
            continue
        hsp.write(raw)
        q0, q1 = min(qs, qe) - 1, max(qs, qe)
        s0, s1 = min(ss, se) - 1, max(ss, se)
        bed.write(f'{q}\t{q0}\t{q1}\n')
        bed.write(f'{q}\t{s0}\t{s1}\n')
PY_INTRA_FILTER

    INTRA_WORKER="$TMP_ROOT/run_intra.sh"
    cat > "$INTRA_WORKER" <<'SH_INTRA'
#!/usr/bin/env bash
set -euo pipefail
fasta="$1"; filter_py="$2"; work_root="$3"; out_root="$4"; min_align="$5"; min_pident="$6"; chunk="$7"
base="$(basename "$chunk" .list)"
work="$work_root/$base"
mkdir -p "$work"
while IFS= read -r contig; do
    [[ -n "$contig" ]] || continue
    samtools faidx "$fasta" "$contig" > "$work/contig.fa"
    blastn -task blastn -query "$work/contig.fa" -subject "$work/contig.fa" \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
        2>> "$out_root/$base.blastn.stderr.log"
done < "$chunk" \
    | python3 "$filter_py" "$out_root/$base.hsp.tsv" "$out_root/$base.bed" "$min_align" "$min_pident"
rm -rf "$work"
: > "$out_root/$base.done"
SH_INTRA
    chmod +x "$INTRA_WORKER"

    echo "[60B] Running INTRA repeat BLAST (each contig against itself): ${REPEAT_JOBS} parallel jobs"
    if [[ -s "$INTRA_MANIFEST" ]]; then
        qc_progress_watch "06_repeats][INTRA BLAST chunks" "$(wc -l < "$INTRA_MANIFEST")" \
            "$INTRA_OUT_DIR" '*.done' &
        progress_pid=$!
        xargs -r -n 1 -P "$REPEAT_JOBS" \
            "$INTRA_WORKER" "$REPEAT_FASTA" "$INTRA_FILTER" "$INTRA_DB_DIR" "$INTRA_OUT_DIR" \
            "$REPEAT_MIN_ALIGN_PAPER" "$REPEAT_PIDENT_PAPER" < "$INTRA_MANIFEST"
        qc_progress_stop "$progress_pid"
        while IFS= read -r f; do cat "$f" >> "$TMP_INTRA_HSP"; done \
            < <(find "$INTRA_OUT_DIR" -name '*.hsp.tsv' -type f | LC_ALL=C sort)
        while IFS= read -r f; do cat "$f" >> "$TMP_INTRA_RAW_BED"; done \
            < <(find "$INTRA_OUT_DIR" -name '*.bed' -type f | LC_ALL=C sort)
    fi

    echo "[60C] Merging INTRA repeat intervals"
    merge_bed_file "$TMP_INTRA_RAW_BED" "$TMP_INTRA_MERGED"

    echo "[60D] Preparing multi-contig MAGs for same-MAG INTER repeat BLAST"
    INTER_BUNDLE_ROOT="$TMP_ROOT/inter_MAGs"
    INTER_OUT_DIR="$TMP_ROOT/inter_results"
    INTER_MANIFEST="$TMP_ROOT/inter_bundles.list"
    mkdir -p "$INTER_BUNDLE_ROOT" "$INTER_OUT_DIR"

    python3 - "$REPEAT_MAG_MAP" "$INTER_BUNDLE_ROOT" "$INTER_MANIFEST" "${REPEAT_FASTA}.fai" <<'PY_MAG_BUNDLES'
import csv
import sys
from collections import defaultdict
from pathlib import Path
map_path, root, manifest, fai = sys.argv[1:]
root = Path(root)
lengths = {}
with open(fai) as handle:
    for line in handle:
        name, length = line.split('\t')[:2]
        lengths[name] = int(length)
with open(map_path, newline='') as handle:
    rows = list(csv.reader(handle, delimiter='\t'))
if rows and len(rows[0]) >= 2 and rows[0][0].strip().lower() in {'mag_id','mag','id'}:
    rows = rows[1:]
mags = defaultdict(list)
for row in rows:
    if len(row) >= 2 and row[0].strip() and row[1].strip():
        mags[row[0].strip()].append(row[1].strip())
bundles = []
for i, mag in enumerate(sorted(mags), 1):
    members = mags[mag]
    if len(members) <= 1:
        continue
    bundle = root / f'MAG_{i:06d}'
    bundle.mkdir(parents=True, exist_ok=True)
    (bundle / 'MAG_ID.txt').write_text(mag + '\n')
    (bundle / 'contigs.txt').write_text('\n'.join(members) + '\n')
    bundles.append((sum(lengths.get(c, 0) for c in members) * len(members), bundle))
paths = [bundle for _, bundle in sorted(bundles, key=lambda item: -item[0])]
with open(manifest, 'w') as out:
    for p in paths:
        out.write(str(p) + '\n')
print(f'Multi-contig MAG bundles: {len(paths)}')
PY_MAG_BUNDLES

    INTER_FILTER="$TMP_ROOT/filter_inter.py"
    cat > "$INTER_FILTER" <<'PY_INTER_FILTER'
import sys
mag_id, hsp_out, bed_out, min_align, min_pident = sys.argv[1:]
min_align, min_pident = int(min_align), float(min_pident)
best = {}
for raw in sys.stdin:
    if not raw.strip():
        continue
    f = raw.rstrip('\n').split('\t')
    if len(f) != 10:
        raise SystemExit(f'ERROR: malformed BLAST row ({len(f)} fields)')
    q, s = f[0], f[1]
    if q == s:
        continue
    pident, alen = float(f[2]), int(float(f[3]))
    if alen <= min_align or pident < min_pident:
        continue
    qs, qe, ss, se = map(int, f[4:8])
    qlen, slen = map(int, f[8:10])
    qlo, qhi = sorted((qs, qe)); slo, shi = sorted((ss, se))
    qstrand = 1 if qe >= qs else -1
    sstrand = 1 if se >= ss else -1
    a = (q, qlo, qhi, qlen, qstrand)
    b = (s, slo, shi, slen, sstrand)
    # Canonicalize the two aligned copies so reciprocal q<->s BLAST rows collapse.
    if (a[0], a[1], a[2]) <= (b[0], b[1], b[2]):
        left, right = a, b
    else:
        left, right = b, a
    orientation = 'direct' if left[4] == right[4] else 'inverted'
    key = (left[0], left[1], left[2], right[0], right[1], right[2], orientation)
    candidate = (pident, alen, left, right, orientation)
    old = best.get(key)
    if old is None or (pident, alen) > (old[0], old[1]):
        best[key] = candidate

with open(hsp_out, 'w') as hsp, open(bed_out, 'w') as bed:
    for key in sorted(best):
        pident, alen, left, right, orientation = best[key]
        q, qlo, qhi, qlen, qstrand = left
        s, slo, shi, slen, sstrand = right
        # Coordinates are output 1-based inclusive in ascending order; orientation
        # is retained in a separate explicit column.
        hsp.write(
            f'{mag_id}\t{q}\t{s}\t{pident:.6f}\t{alen}\t{qlo}\t{qhi}\t'
            f'{slo}\t{shi}\t{qlen}\t{slen}\t{orientation}\n'
        )
        bed.write(f'{q}\t{qlo-1}\t{qhi}\n')
        bed.write(f'{s}\t{slo-1}\t{shi}\n')
PY_INTER_FILTER

    INTER_WORKER="$TMP_ROOT/run_inter.sh"
    cat > "$INTER_WORKER" <<'SH_INTER'
#!/usr/bin/env bash
set -euo pipefail
fasta="$1"; filter_py="$2"; out_root="$3"; min_align="$4"; min_pident="$5"; bundle="$6"
base="$(basename "$bundle")"
mag_id="$(<"$bundle/MAG_ID.txt")"
mapfile -t contigs < "$bundle/contigs.txt"
(( ${#contigs[@]} > 1 )) || { : > "$out_root/$base.done"; exit 0; }
mag_fa="$bundle/MAG.fa"
samtools faidx "$fasta" "${contigs[@]}" > "$mag_fa"
samtools faidx "$mag_fa"
for i in "${!contigs[@]}"; do
    others=("${contigs[@]:0:i}" "${contigs[@]:i+1}")
    samtools faidx "$mag_fa" "${contigs[i]}" > "$bundle/query.fa"
    samtools faidx "$mag_fa" "${others[@]}" > "$bundle/others.fa"
    blastn -task blastn -query "$bundle/query.fa" -subject "$bundle/others.fa" \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
        2>> "$out_root/$base.blastn.stderr.log"
done | python3 "$filter_py" "$mag_id" "$out_root/$base.hsp.tsv" "$out_root/$base.bed" "$min_align" "$min_pident"
rm -f "$bundle/query.fa" "$bundle/others.fa" "$mag_fa" "$mag_fa.fai"
: > "$out_root/$base.done"
SH_INTER
    chmod +x "$INTER_WORKER"

    echo "[60E] Running same-MAG INTER repeat BLAST (each contig against the other contigs of its MAG): ${REPEAT_JOBS} parallel jobs"
    if [[ -s "$INTER_MANIFEST" ]]; then
        qc_progress_watch "06_repeats][INTER BLAST MAGs" "$(wc -l < "$INTER_MANIFEST")" \
            "$INTER_OUT_DIR" '*.done' &
        progress_pid=$!
        xargs -r -n 1 -P "$REPEAT_JOBS" \
            "$INTER_WORKER" "$REPEAT_FASTA" "$INTER_FILTER" "$INTER_OUT_DIR" \
            "$REPEAT_MIN_ALIGN_PAPER" "$REPEAT_PIDENT_PAPER" < "$INTER_MANIFEST"
        qc_progress_stop "$progress_pid"
        while IFS= read -r f; do cat "$f" >> "$TMP_INTER_HSP"; done \
            < <(find "$INTER_OUT_DIR" -name '*.hsp.tsv' -type f | LC_ALL=C sort)
        while IFS= read -r f; do cat "$f" >> "$TMP_INTER_RAW_BED"; done \
            < <(find "$INTER_OUT_DIR" -name '*.bed' -type f | LC_ALL=C sort)
    fi

    echo "[60F] Merging INTER intervals and INTRA+INTER TOTAL union"
    merge_bed_file "$TMP_INTER_RAW_BED" "$TMP_INTER_MERGED"
    cat "$TMP_INTRA_RAW_BED" "$TMP_INTER_RAW_BED" > "$TMP_TOTAL_RAW_BED"
    merge_bed_file "$TMP_TOTAL_RAW_BED" "$TMP_TOTAL_MERGED"

    echo "[60G] Writing INTRA/INTER evidence summaries for final normalization"
    python3 "$SCRIPT_DIR/../scripts/repeats_tables.py" evidence \
        "${REPEAT_FASTA}.fai" "$REPEAT_MAG_MAP" \
        "$TMP_INTRA_MERGED" "$TMP_INTER_MERGED" "$TMP_TOTAL_MERGED" \
        "$TMP_INTRA_HSP" "$TMP_INTER_HSP" \
        "$REPEAT_SUM" "$REPEAT_MAG_SUM" \
        "$INTRA_REGION_BED" "$INTRA_REGION_TSV" \
        "$INTER_REGION_BED" "$INTER_REGION_TSV" \
        "$REPEAT_REGION_BED" "$REPEAT_REGION_TSV" \
        "$REPEAT_DATASET_SUM" "$REPEAT_CONTIG_HIGH" "$REPEAT_GENOME_HIGH"

    mv -f "$TMP_INTRA_HSP" "$SELF_BLAST"
    mv -f "$TMP_INTER_HSP" "$INTER_BLAST"
    qc_write_result_cache_atomic "$REPEAT_CACHE_META" "$REPEAT_FINGERPRINT" "$REPEAT_SUM"

    rm -rf "$TMP_ROOT"
    trap - EXIT
fi

echo "[60H] Normalizing final repeat interpretation: contig=INTRA only; MAG=INTRA+INTER union"
python3 "$SCRIPT_DIR/../scripts/repeats_tables.py" normalize \
    "$REPEAT_SUM" "$REPEAT_MAG_SUM" "$REPEAT_DATASET_SUM" \
    "$REPEAT_CONTIG_HIGH" "$REPEAT_GENOME_HIGH"


read -r -a REPEAT_ASSEMBLER_FILES <<< "$ASSEMBLER_FASTAS"
((${#REPEAT_ASSEMBLER_FILES[@]} > 0)) || {
    echo "ERROR: ASSEMBLER_FASTAS is empty"
    exit 1
}

REPEAT_SOURCE_FASTAS=()
for fasta_name in "${REPEAT_ASSEMBLER_FILES[@]}"; do
    [[ "$fasta_name" != */* ]] || {
        echo "ERROR: ASSEMBLER_FASTAS must contain basenames only: $fasta_name"
        exit 1
    }
    source_fasta="$ASSEMBLER_ASSEMBLY_DIR/$fasta_name"
    [[ -s "$source_fasta" ]] || {
        echo "ERROR: assembler FASTA missing or empty: $source_fasta"
        exit 1
    }
    if [[ ! -s "${source_fasta}.fai" || "$source_fasta" -nt "${source_fasta}.fai" ]]; then
        samtools faidx "$source_fasta"
    fi
    REPEAT_SOURCE_FASTAS+=("$source_fasta")
done

python3 "$SCRIPT_DIR/../scripts/repeats_tables.py" assembler \
    "${REPEAT_FASTA}.fai" "$REPEAT_SUM" "$REPEAT_MAG_SUM" \
    "$REPEAT_DATASET_SUM" "${REPEAT_SOURCE_FASTAS[@]}"


rm -f "$REPEAT_DIR/dataset_intracontig_repeat_breadth.final.tsv"

qc_write_result_cache_atomic "$REPEAT_CACHE_META" "$REPEAT_FINGERPRINT" "$REPEAT_SUM"

if [[ -s "$SELF_BLAST" && -s "$INTRA_REGION_TSV" ]]; then
    python3 "$SCRIPT_DIR/../scripts/repeats_tables.py" consolidate "$REPEAT_MAG_MAP" \
        "$SELF_BLAST" "$INTER_BLAST" "$INTRA_REGION_TSV" "$INTER_REGION_TSV" "$REPEAT_REGION_TSV" \
        "$REPEAT_HSP" "$REPEAT_REGIONS"
fi
rm -rf "$REPEAT_WORK_DIR"

for f in "$REPEAT_SUM" "$REPEAT_MAG_SUM" "$REPEAT_HSP" "$REPEAT_REGIONS" \
         "$REPEAT_DATASET_SUM"; do
    [[ -s "$f" ]] || { echo "ERROR: expected repeat output missing/empty: $f"; exit 1; }
done

echo
echo "[06_repeats] Final outputs"
echo "Contig summary (INTRA)            : $REPEAT_SUM"
echo "MAG summary (INTRA+INTER union)   : $REPEAT_MAG_SUM"
echo "Dataset summary (INTRA)           : $REPEAT_DATASET_SUM"
echo "Repeat regions (INTRA/INTER/TOTAL): $REPEAT_REGIONS"
echo "BLAST HSPs (INTRA/INTER)          : $REPEAT_HSP"
echo "Assembler information             : Assembler + assembler-specific INTRA excessive-repeat contig counts"
echo "Contig interpretation             : EXCESSIVE_REPEAT_CONTIG=Y means INTRA repeat fraction >= contig threshold"
echo "MAG excessive-contig annotation : HAS_EXCESSIVE_INTRA_REPEAT_CONTIG=Y means >=1 member contig has INTRA repeat fraction >= contig threshold"
echo "MAG interpretation                : EXCESSIVE_REPEAT_MAG=Y means MAG-wide INTRA+INTER union fraction >= genome threshold"
echo "[06_repeats] Finished: $(date)"
