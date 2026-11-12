#!/usr/bin/env bash
set -euo pipefail

# 06_repeats.sh
#
# Repeat screen for long-read assembly contigs and final MAGs.
#
# Intra-contig repeat
#   - BLASTN self-hit on the same contig (qseqid == sseqid)
#   - exact diagonal self-hit excluded
#   - alignment length >200 bp and nucleotide identity >=80%
#
# Inter-contig repeat
#   - analyzed only between different contigs belonging to the same final MAG
#   - qseqid != sseqid
#   - alignment length >200 bp and nucleotide identity >=80%
#   - reciprocal BLAST HSPs are collapsed to one HSP
#
# Final breadth and calls
#   - Contig-level repeat statistics and EXCESSIVE_REPEAT_CONTIG are based ONLY
#     on INTRA-contig repeats.
#   - EXCESSIVE_REPEAT_CONTIG=Y when the INTRA-contig repeat fraction >=70%.
#   - INTER-contig repeats are calculated only between different contigs within
#     the same final multi-contig MAG.
#   - For MAG-level total repeat content, INTRA and INTER intervals are merged
#     and overlapping intervals are counted once.
#   - HAS_EXCESSIVE_INTRA_REPEAT_CONTIG=Y when at least one member contig has
#     an INTRA-contig repeat fraction >=70%; this is an annotation only.
#   - EXCESSIVE_REPEAT_MAG=Y when the MAG-wide INTRA+INTER union fraction >=70%.
#   - source assembler is reported for contigs/MAGs and excessive-repeat contigs
#     are counted separately for each configured assembler in the dataset summary
#
# Compatibility
#   - self_blast.paper_same_contig.filtered.tsv is retained unchanged in format
#     for 08_false_duplication.sh.

SCRIPT_ID="06_repeats"
# The evidence method ID is part of the repeat cache key.
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
REPEAT_CACHE_META="${REPEAT_CACHE_META:-${REPEAT_SUM}.inputs.sha256}"

# Paper thresholds.
REPEAT_MIN_ALIGN_PAPER="${REPEAT_MIN_ALIGN_PAPER:-200}"
REPEAT_PIDENT_PAPER="${REPEAT_PIDENT_PAPER:-80}"
REPEAT_CONTIG_HIGH="${REPEAT_CONTIG_HIGH:-${REPEAT_HIGH:-0.70}}"
REPEAT_GENOME_HIGH="${REPEAT_GENOME_HIGH:-${REPEAT_HIGH:-0.70}}"

# Parallelism.
REPEAT_CHUNK_BP="${REPEAT_CHUNK_BP:-10000000}"
REPEAT_BLAST_THREADS_PER_JOB="${REPEAT_BLAST_THREADS_PER_JOB:-2}"
REPEAT_JOBS="${REPEAT_JOBS:-$(( THREADS / REPEAT_BLAST_THREADS_PER_JOB ))}"
(( REPEAT_JOBS >= 1 )) || REPEAT_JOBS=1
REPEAT_SORT_THREADS="${REPEAT_SORT_THREADS:-$THREADS}"
REPEAT_SORT_MEMORY="${REPEAT_SORT_MEMORY:-50%}"
REPEAT_MIN_CONTIG_LEN="${REPEAT_MIN_CONTIG_LEN:-0}"

if (( REPEAT_JOBS * REPEAT_BLAST_THREADS_PER_JOB > THREADS )); then
    echo "ERROR: REPEAT_JOBS * REPEAT_BLAST_THREADS_PER_JOB exceeds THREADS"
    exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Main outputs. INTRA evidence is used for contig-level reporting/calls.
# INTER/TOTAL evidence is retained for MAG-level repeat breadth calculation.
SELF_BLAST="${REPEAT_SELF_BLAST:-$REPEAT_DIR/self_blast.paper_same_contig.filtered.tsv}"
INTER_BLAST="${REPEAT_INTER_BLAST:-$REPEAT_DIR/inter_blast.paper_same_MAG_different_contig.filtered.tsv}"
REPEAT_MAG_MAP="${REPEAT_MAG_MAP:-$MAG_CONTIG_MAP}"
REPEAT_MAG_SUM="${REPEAT_MAG_SUM:-$REPEAT_DIR/MAG_repeat_breadth.final.tsv}"
REPEAT_DATASET_SUM="${REPEAT_DATASET_SUM:-$REPEAT_DIR/dataset_repeat_breadth.final.tsv}"
REPEAT_DATASET_SUM_INTUITIVE="${REPEAT_DATASET_SUM_INTUITIVE:-$REPEAT_DIR/dataset_repeat_breadth.final.tsv}"

# Existing repeat_regions.final.* now means TOTAL union and remains the simple
# generic final-region path. Explicit INTRA and INTER files are also written.
REPEAT_REGION_BED="${REPEAT_REGION_BED:-$REPEAT_DIR/repeat_regions.final.bed}"
REPEAT_REGION_TSV="${REPEAT_REGION_TSV:-$REPEAT_DIR/repeat_regions.final.tsv}"
INTRA_REGION_BED="${INTRA_REGION_BED:-$REPEAT_DIR/intra_repeat_regions.final.bed}"
INTRA_REGION_TSV="${INTRA_REGION_TSV:-$REPEAT_DIR/intra_repeat_regions.final.tsv}"
INTER_REGION_BED="${INTER_REGION_BED:-$REPEAT_DIR/inter_repeat_regions.final.bed}"
INTER_REGION_TSV="${INTER_REGION_TSV:-$REPEAT_DIR/inter_repeat_regions.final.tsv}"

mkdir -p "$OUT/00_log" "$REPEAT_DIR" \
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
echo "REPEAT_JOBS                : $REPEAT_JOBS"
echo "REPEAT_BLAST_THREADS/JOB   : $REPEAT_BLAST_THREADS_PER_JOB"
echo "FORCE_RERUN_REPEAT         : $FORCE_RERUN_REPEAT"
echo

for tool in samtools makeblastdb blastn bedtools python3 sort xargs; do
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

# Validate FASTA and MAG map.

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

# Cache fingerprint. INTER analysis depends on both assembly and MAG membership.

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

# Cache validity is tied to the expensive repeat-evidence generation, not to the
if [[ "$FORCE_RERUN_REPEAT" != "Y" \
      && "$CACHED_REPEAT_FINGERPRINT" == "$REPEAT_FINGERPRINT" \
      && -s "$REPEAT_SUM" && -s "$REPEAT_MAG_SUM" \
      && -s "$SELF_BLAST" && -s "$INTER_BLAST" \
      && -e "$INTRA_REGION_BED" && -s "$INTRA_REGION_TSV" \
      && -e "$INTER_REGION_BED" && -s "$INTER_REGION_TSV" \
      && -e "$REPEAT_REGION_BED" && -s "$REPEAT_REGION_TSV" ]]; then
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

    ###########################################################################
    # 60A-C. INTRA-CONTIG repeat BLAST.
    ###########################################################################

    echo "[60A] Preparing intact-contig chunks for INTRA repeat BLAST"
    INTRA_CHUNK_DIR="$TMP_ROOT/intra_chunks"
    INTRA_DB_DIR="$TMP_ROOT/intra_dbs"
    INTRA_OUT_DIR="$TMP_ROOT/intra_results"
    INTRA_MANIFEST="$TMP_ROOT/intra_chunks.list"
    mkdir -p "$INTRA_CHUNK_DIR" "$INTRA_DB_DIR" "$INTRA_OUT_DIR"

    EFFECTIVE_MIN_CONTIG_LEN="$REPEAT_MIN_CONTIG_LEN"
    if (( EFFECTIVE_MIN_CONTIG_LEN <= REPEAT_MIN_ALIGN_PAPER )); then
        EFFECTIVE_MIN_CONTIG_LEN=$(( REPEAT_MIN_ALIGN_PAPER + 1 ))
    fi

    python3 - "$REPEAT_FASTA" "$INTRA_CHUNK_DIR" "$INTRA_MANIFEST" \
        "$REPEAT_CHUNK_BP" "$EFFECTIVE_MIN_CONTIG_LEN" <<'PY_SPLIT'
import sys
from pathlib import Path
fasta, out_dir, manifest, target_bp, min_len = sys.argv[1:]
out_dir = Path(out_dir)
target_bp, min_len = int(target_bp), int(min_len)
idx, chunk_bp, out = 0, 0, None
paths = []
def open_chunk():
    global idx, chunk_bp, out
    idx += 1
    path = out_dir / f'chunk_{idx:06d}.fa'
    out = path.open('w')
    paths.append(path)
    chunk_bp = 0
def emit(header, seq):
    global chunk_bp, out
    if len(seq) < min_len:
        return
    if out is None or (chunk_bp > 0 and chunk_bp + len(seq) > target_bp):
        if out is not None:
            out.close()
        open_chunk()
    out.write(header if header.endswith('\n') else header + '\n')
    for i in range(0, len(seq), 80):
        out.write(seq[i:i+80] + '\n')
    chunk_bp += len(seq)
header, parts = None, []
with open(fasta) as handle:
    for line in handle:
        if line.startswith('>'):
            if header is not None:
                emit(header, ''.join(parts))
            header, parts = line, []
        else:
            parts.append(line.strip())
if header is not None:
    emit(header, ''.join(parts))
if out is not None:
    out.close()
with open(manifest, 'w') as handle:
    for path in paths:
        handle.write(str(path) + '\n')
print(f'Eligible contigs written to {len(paths)} chunks.')
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
filter_py="$1"; db_root="$2"; out_root="$3"; threads="$4"; min_align="$5"; min_pident="$6"; chunk="$7"
base="$(basename "$chunk" .fa)"
db="$db_root/$base/db"
mkdir -p "$(dirname "$db")"
makeblastdb -in "$chunk" -dbtype nucl -out "$db" > "$out_root/$base.makeblastdb.log" 2>&1
blastn -task blastn -query "$chunk" -db "$db" \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
    -num_threads "$threads" 2> "$out_root/$base.blastn.stderr.log" \
    | python3 "$filter_py" "$out_root/$base.hsp.tsv" "$out_root/$base.bed" "$min_align" "$min_pident"
SH_INTRA
    chmod +x "$INTRA_WORKER"

    echo "[60B] Running INTRA repeat BLAST: ${REPEAT_JOBS} jobs x ${REPEAT_BLAST_THREADS_PER_JOB} threads"
    if [[ -s "$INTRA_MANIFEST" ]]; then
        xargs -r -n 1 -P "$REPEAT_JOBS" \
            "$INTRA_WORKER" "$INTRA_FILTER" "$INTRA_DB_DIR" "$INTRA_OUT_DIR" \
            "$REPEAT_BLAST_THREADS_PER_JOB" "$REPEAT_MIN_ALIGN_PAPER" \
            "$REPEAT_PIDENT_PAPER" < "$INTRA_MANIFEST"
        while IFS= read -r f; do cat "$f" >> "$TMP_INTRA_HSP"; done \
            < <(find "$INTRA_OUT_DIR" -name '*.hsp.tsv' -type f | LC_ALL=C sort)
        while IFS= read -r f; do cat "$f" >> "$TMP_INTRA_RAW_BED"; done \
            < <(find "$INTRA_OUT_DIR" -name '*.bed' -type f | LC_ALL=C sort)
    fi

    echo "[60C] Merging INTRA repeat intervals"
    merge_bed_file "$TMP_INTRA_RAW_BED" "$TMP_INTRA_MERGED"

    ###########################################################################
    # 60D-F. INTER-CONTIG repeat BLAST within each multi-contig MAG.
    ###########################################################################

    echo "[60D] Preparing multi-contig MAGs for same-MAG INTER repeat BLAST"
    INTER_BUNDLE_ROOT="$TMP_ROOT/inter_MAGs"
    INTER_OUT_DIR="$TMP_ROOT/inter_results"
    INTER_MANIFEST="$TMP_ROOT/inter_bundles.list"
    mkdir -p "$INTER_BUNDLE_ROOT" "$INTER_OUT_DIR"

    python3 - "$REPEAT_MAG_MAP" "$INTER_BUNDLE_ROOT" "$INTER_MANIFEST" <<'PY_MAG_BUNDLES'
import csv
import sys
from collections import defaultdict
from pathlib import Path
map_path, root, manifest = sys.argv[1:]
root = Path(root)
with open(map_path, newline='') as handle:
    rows = list(csv.reader(handle, delimiter='\t'))
if rows and len(rows[0]) >= 2 and rows[0][0].strip().lower() in {'mag_id','mag','id'}:
    rows = rows[1:]
mags = defaultdict(list)
for row in rows:
    if len(row) >= 2 and row[0].strip() and row[1].strip():
        mags[row[0].strip()].append(row[1].strip())
paths = []
for i, mag in enumerate(sorted(mags), 1):
    members = mags[mag]
    if len(members) <= 1:
        continue
    bundle = root / f'MAG_{i:06d}'
    bundle.mkdir(parents=True, exist_ok=True)
    (bundle / 'MAG_ID.txt').write_text(mag + '\n')
    (bundle / 'contigs.txt').write_text('\n'.join(members) + '\n')
    paths.append(bundle)
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
fasta="$1"; filter_py="$2"; out_root="$3"; threads="$4"; min_align="$5"; min_pident="$6"; bundle="$7"
base="$(basename "$bundle")"
mag_id="$(<"$bundle/MAG_ID.txt")"
mapfile -t contigs < "$bundle/contigs.txt"
(( ${#contigs[@]} > 1 )) || exit 0
mag_fa="$bundle/MAG.fa"
db="$bundle/db/MAG"
mkdir -p "$(dirname "$db")"
samtools faidx "$fasta" "${contigs[@]}" > "$mag_fa"
makeblastdb -in "$mag_fa" -dbtype nucl -out "$db" > "$out_root/$base.makeblastdb.log" 2>&1
blastn -task blastn -query "$mag_fa" -db "$db" \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
    -num_threads "$threads" 2> "$out_root/$base.blastn.stderr.log" \
    | python3 "$filter_py" "$mag_id" "$out_root/$base.hsp.tsv" "$out_root/$base.bed" "$min_align" "$min_pident"
SH_INTER
    chmod +x "$INTER_WORKER"

    echo "[60E] Running same-MAG INTER repeat BLAST: ${REPEAT_JOBS} jobs x ${REPEAT_BLAST_THREADS_PER_JOB} threads"
    if [[ -s "$INTER_MANIFEST" ]]; then
        xargs -r -n 1 -P "$REPEAT_JOBS" \
            "$INTER_WORKER" "$REPEAT_FASTA" "$INTER_FILTER" "$INTER_OUT_DIR" \
            "$REPEAT_BLAST_THREADS_PER_JOB" "$REPEAT_MIN_ALIGN_PAPER" \
            "$REPEAT_PIDENT_PAPER" < "$INTER_MANIFEST"
        while IFS= read -r f; do cat "$f" >> "$TMP_INTER_HSP"; done \
            < <(find "$INTER_OUT_DIR" -name '*.hsp.tsv' -type f | LC_ALL=C sort)
        while IFS= read -r f; do cat "$f" >> "$TMP_INTER_RAW_BED"; done \
            < <(find "$INTER_OUT_DIR" -name '*.bed' -type f | LC_ALL=C sort)
    fi

    echo "[60F] Merging INTER intervals and INTRA+INTER TOTAL union"
    merge_bed_file "$TMP_INTER_RAW_BED" "$TMP_INTER_MERGED"
    cat "$TMP_INTRA_RAW_BED" "$TMP_INTER_RAW_BED" > "$TMP_TOTAL_RAW_BED"
    merge_bed_file "$TMP_TOTAL_RAW_BED" "$TMP_TOTAL_MERGED"

    ###########################################################################
    # 60G. Comprehensive contig/MAG/region summaries.
    ###########################################################################

    echo "[60G] Writing INTRA/INTER evidence summaries for final normalization"
    python3 "$SCRIPT_DIR/../scripts/repeats_summarize_repeat_evidence.py" \
        "${REPEAT_FASTA}.fai" "$REPEAT_MAG_MAP" \
        "$TMP_INTRA_MERGED" "$TMP_INTER_MERGED" "$TMP_TOTAL_MERGED" \
        "$TMP_INTRA_HSP" "$TMP_INTER_HSP" \
        "$REPEAT_SUM" "$REPEAT_MAG_SUM" \
        "$INTRA_REGION_BED" "$INTRA_REGION_TSV" \
        "$INTER_REGION_BED" "$INTER_REGION_TSV" \
        "$REPEAT_REGION_BED" "$REPEAT_REGION_TSV" \
        "$REPEAT_DATASET_SUM_INTUITIVE" "$REPEAT_CONTIG_HIGH" "$REPEAT_GENOME_HIGH"

    # Keep the old configured dataset path as a compatibility copy if it differs
    # from the new intuitive filename.
    if [[ "$REPEAT_DATASET_SUM" != "$REPEAT_DATASET_SUM_INTUITIVE" ]]; then
        cp -f "$REPEAT_DATASET_SUM_INTUITIVE" "$REPEAT_DATASET_SUM"
    fi

    mv -f "$TMP_INTRA_HSP" "$SELF_BLAST"
    mv -f "$TMP_INTER_HSP" "$INTER_BLAST"
    qc_write_result_cache_atomic "$REPEAT_CACHE_META" "$REPEAT_FINGERPRINT" "$REPEAT_SUM"

    rm -rf "$TMP_ROOT"
    trap - EXIT
fi

# Normalize final interpretation/output schema.
#
# IMPORTANT
#   * Contig-level repeat statistics/calls are INTRA-only.
#   * INTER is not classified at the contig level.
#   * INTER is retained only in the MAG table to calculate MAG-wide TOTAL repeat
#     breadth as union(INTRA, INTER), with overlap already removed upstream.
#
# can be corrected without rerunning BLASTN.

echo "[60H] Normalizing final repeat interpretation: contig=INTRA only; MAG=INTRA+INTER union"
python3 "$SCRIPT_DIR/../scripts/repeats_normalize_repeat_tables.py" \
    "$REPEAT_SUM" "$REPEAT_MAG_SUM" "$REPEAT_DATASET_SUM_INTUITIVE" \
    "$REPEAT_CONTIG_HIGH" "$REPEAT_GENOME_HIGH"

if [[ "$REPEAT_DATASET_SUM" != "$REPEAT_DATASET_SUM_INTUITIVE" ]]; then
    cp -f "$REPEAT_DATASET_SUM_INTUITIVE" "$REPEAT_DATASET_SUM"
fi

# Add source-assembler labels and assembler-specific excessive-repeat counts.
#
# This is output-only annotation. It runs after the repeat cache decision so an
# existing valid BLAST/repeat result does not need to be recomputed merely to
# add assembler provenance columns.

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
        echo "ERROR: source assembler FASTA missing or empty: $source_fasta"
        exit 1
    }
    if [[ ! -s "${source_fasta}.fai" || "$source_fasta" -nt "${source_fasta}.fai" ]]; then
        samtools faidx "$source_fasta"
    fi
    REPEAT_SOURCE_FASTAS+=("$source_fasta")
done

python3 "$SCRIPT_DIR/../scripts/repeats_summarize_repeat_by_assembler.py" \
    "${REPEAT_FASTA}.fai" "$REPEAT_SUM" "$REPEAT_MAG_SUM" \
    "$REPEAT_DATASET_SUM_INTUITIVE" "${REPEAT_SOURCE_FASTAS[@]}"

# Keep the legacy configured dataset filename synchronized with the intuitive
# dataset summary after assembler-specific columns are added.
if [[ "$REPEAT_DATASET_SUM" != "$REPEAT_DATASET_SUM_INTUITIVE" ]]; then
    cp -f "$REPEAT_DATASET_SUM_INTUITIVE" "$REPEAT_DATASET_SUM"
fi

# Remove the obsolete compatibility filename from older repeat revisions so a
# stale file is not mistaken for the current dataset summary.
LEGACY_REPEAT_DATASET_SUM="$REPEAT_DIR/dataset_intracontig_repeat_breadth.final.tsv"
if [[ "$LEGACY_REPEAT_DATASET_SUM" != "$REPEAT_DATASET_SUM_INTUITIVE" ]]; then
    rm -f "$LEGACY_REPEAT_DATASET_SUM"
fi


# Refresh cache metadata after final schema normalization/assembler annotation.
# The first field is the evidence fingerprint and is used to decide whether
# expensive BLAST evidence can be reused on the next run.
qc_write_result_cache_atomic "$REPEAT_CACHE_META" "$REPEAT_FINGERPRINT" "$REPEAT_SUM"

# Final sanity checks and concise output map.

for f in "$REPEAT_SUM" "$REPEAT_MAG_SUM" "$SELF_BLAST" "$INTER_BLAST" \
         "$INTRA_REGION_TSV" "$INTER_REGION_TSV" "$REPEAT_REGION_TSV" \
         "$REPEAT_DATASET_SUM_INTUITIVE"; do
    [[ -s "$f" ]] || { echo "ERROR: expected repeat output missing/empty: $f"; exit 1; }
done

echo
echo "[06_repeats] Final outputs"
echo "Contig INTRA-only summary         : $REPEAT_SUM"
echo "MAG INTRA+INTER union summary     : $REPEAT_MAG_SUM"
echo "INTRA HSPs (used by 65)           : $SELF_BLAST"
echo "INTER same-MAG HSPs               : $INTER_BLAST"
echo "INTRA repeat regions              : $INTRA_REGION_TSV"
echo "INTER repeat regions              : $INTER_REGION_TSV"
echo "TOTAL union repeat regions        : $REPEAT_REGION_TSV"
echo "Dataset INTRA-only summary        : $REPEAT_DATASET_SUM_INTUITIVE"
echo "Assembler information             : source_assembler + assembler-specific INTRA excessive-repeat contig counts"
echo "Contig interpretation             : EXCESSIVE_REPEAT_CONTIG=Y means INTRA repeat fraction >= contig threshold"
echo "MAG excessive-contig annotation : HAS_EXCESSIVE_INTRA_REPEAT_CONTIG=Y means >=1 member contig has INTRA repeat fraction >= contig threshold"
echo "MAG interpretation                : EXCESSIVE_REPEAT_MAG=Y means MAG-wide INTRA+INTER union fraction >= genome threshold"
echo "[06_repeats] Finished: $(date)"
