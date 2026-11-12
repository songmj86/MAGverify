#!/usr/bin/env bash
set -Eeuo pipefail

# Profile assembly alleles across source co-assemblies and MAGs
# using all mapped HiFi alignments.

SCRIPT_ID="09_assembly_variants"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

: "${OUT:?ERROR: OUT is not set}"
: "${READS:?ERROR: READS is not set}"
: "${ANVIO_ENV_PREFIX:?ERROR: ANVIO_ENV_PREFIX is not set}"
: "${ANVIO_THREADS:?ERROR: ANVIO_THREADS is not set}"

MAG_CONTIG_MAP="${MAG_CONTIG_MAP:-NA}"
PER_ASSEMBLER_MANIFEST="${PER_ASSEMBLER_MANIFEST:-$ASSEMBLER_MAPPING_MANIFEST}"
ASSEMBLY_VARIANT_SAMPLE_ID="${ASSEMBLY_VARIANT_SAMPLE_ID:-HiFi_reads}"
FORCE_RERUN_ASSEMBLY_VARIANTS="${FORCE_RERUN_ASSEMBLY_VARIANTS:-N}"
# Variant profiling assesses local read support and therefore uses all mapped
# HiFi alignments (primary, secondary and supplementary). The mapped BAM is
# prepared once by 01_mapping.sh with unmapped records excluded; no HiFi
# remapping is performed here.
ASSEMBLY_VARIANT_HIFI_BAM_COLUMN="HiFi_all_BAM"
ASSEMBLY_VARIANT_ALIGNMENT_POPULATION="all_mapped"
ASSEMBLY_VARIANT_PROFILE_ALIGNMENT_CLASSES="primary_secondary_supplementary"

# HiFi allele profiling is performed once per complete assembler-specific
# co-assembly. Assemblers are processed serially; there is no contig batching
# and no concurrent anvi-profile job. anvi-profile starts with 100 threads.
# The 900-GB runtime memory guard remains active. Only if that guard is reached
# is the failed assembler retried with fewer threads, without changing any
# analysis criterion.
ASSEMBLY_VARIANT_MEMORY_BUDGET_GB="${ASSEMBLY_VARIANT_MEMORY_BUDGET_GB:-900}"
ASSEMBLY_VARIANT_MEMORY_RESERVE_GB="${ASSEMBLY_VARIANT_MEMORY_RESERVE_GB:-64}"
ALLELE_ANVIO_MEMORY_POLL_SECONDS="${ALLELE_ANVIO_MEMORY_POLL_SECONDS:-5}"
# lib_config.sh previously used 0 for memory-derived auto mode. Interpret 0 as
# the new fixed default of 100 threads. A positive environment/config override
# is still accepted explicitly.
ALLELE_ANVIO_FULL_PROFILE_THREADS="${ALLELE_ANVIO_FULL_PROFILE_THREADS:-100}"
[[ "$ALLELE_ANVIO_FULL_PROFILE_THREADS" == "0" ]] && ALLELE_ANVIO_FULL_PROFILE_THREADS=100
ALLELE_ANVIO_PREPROCESS_THREADS="${ALLELE_ANVIO_PREPROCESS_THREADS:-70}"
ASSEMBLY_VARIANT_KEEP_FULL_PROFILE="${ASSEMBLY_VARIANT_KEEP_FULL_PROFILE:-N}"
ASSEMBLY_VARIANT_DROP_FILE_CACHE="${ASSEMBLY_VARIANT_DROP_FILE_CACHE:-Y}"
ASSEMBLY_VARIANT_SKIP_GENE_CALLING="${ASSEMBLY_VARIANT_SKIP_GENE_CALLING:-N}"
ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES="${ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES:-Y}"

# Final HiFi profiling cache. All-mapped profiling uses a distinct cache root
# so earlier primary+supplementary or secondary-excluded PROFILE.db files are
# never reused for this analysis.
ASSEMBLY_VARIANT_HIFI_PROFILE_ROOT="${VARIANT_DIR}/hifi_all_mapped"
CANONICAL_HIFI_ANVIO_DIR="$ASSEMBLY_VARIANT_HIFI_PROFILE_ROOT/anvio"
if [[ -z "${ALLELE_ANVIO_DIR+x}" ]]; then
    ALLELE_ANVIO_DIR="$CANONICAL_HIFI_ANVIO_DIR"
fi
ALLELE_BY_ASSEMBLER_DIR="${ALLELE_BY_ASSEMBLER_DIR:-$ALLELE_ANVIO_DIR/by_assembler}"
ALLELE_PROFILE_MANIFEST="${ALLELE_PROFILE_MANIFEST:-${VARIANT_DIR}/anvio_profiles.by_assembler.tsv}"
ALLELE_SOURCE_ROUTING="${ALLELE_SOURCE_ROUTING:-$ALLELE_PROFILE_MANIFEST}"

# Assembly-wide outputs.
ALLELE_VARIANTS="${ALLELE_VARIANTS:-${VARIANT_DIR}/assembly_allele_disagreements.tsv}"
ALLELE_CONTIG_ALL="${ALLELE_CONTIG_ALL:-${VARIANT_DIR}/contig_assembly_allele_support.all_contigs.tsv}"
ALLELE_SUMMARY_ALL="${ALLELE_SUMMARY_ALL:-${VARIANT_DIR}/assembly_allele_support.all_contigs.summary.tsv}"
ALLELE_ASSEMBLER_SUMMARY="${ALLELE_ASSEMBLER_SUMMARY:-${VARIANT_DIR}/assembly_allele_support.by_assembler.tsv}"

# Final-MAG subset outputs. Existing names are retained for downstream compatibility.
ALLELE_CONTIG="${ALLELE_CONTIG:-${VARIANT_DIR}/contig_assembly_allele_support.final.tsv}"
ALLELE_SUMMARY="${ALLELE_SUMMARY:-${VARIANT_DIR}/assembly_allele_support.summary.tsv}"
ALLELE_VARIANTS_MAG="${ALLELE_VARIANTS_MAG:-${VARIANT_DIR}/assembly_allele_disagreements.by_MAG.tsv}"
ALLELE_MAG="${ALLELE_MAG:-${VARIANT_DIR}/assembly_allele_support.by_MAG.tsv}"
ALLELE_CACHE_META="${ALLELE_CACHE_META:-${VARIANT_DIR}/assembly_allele_support.inputs.sha256}"

ANVIO_REFORMAT_FASTA="$ANVIO_ENV_PREFIX/bin/anvi-script-reformat-fasta"
ANVIO_REFORMAT_BAM="$ANVIO_ENV_PREFIX/bin/anvi-script-reformat-bam"
ANVIO_GEN_CONTIGS_DB="$ANVIO_ENV_PREFIX/bin/anvi-gen-contigs-database"
ANVIO_PROFILE="$ANVIO_ENV_PREFIX/bin/anvi-profile"
ANVIO_EXPORT_TABLE="$ANVIO_ENV_PREFIX/bin/anvi-export-table"
ANVIO_PYTHON="$ANVIO_ENV_PREFIX/bin/python"
PARSER="$SCRIPT_DIR/../scripts/variants_summarize_assembly_variants.py"
MERGER="$SCRIPT_DIR/../scripts/variants_merge_assembly_variant_assembler_results.py"
MAG_SUMMARIZER="$SCRIPT_DIR/../scripts/variants_summarize_assembly_variants_by_mag.py"

mkdir -p "$OUT/00_log" "${VARIANT_DIR}" "$ALLELE_BY_ASSEMBLER_DIR"
LOG="$OUT/00_log/09_assembly_variants.log"
exec > >(tee -a "$LOG") 2>&1

for path in "$ANVIO_REFORMAT_FASTA" "$ANVIO_REFORMAT_BAM" "$ANVIO_GEN_CONTIGS_DB" \
            "$ANVIO_PROFILE" "$ANVIO_EXPORT_TABLE" "$ANVIO_PYTHON"; do
    [[ -x "$path" ]] || { echo "ERROR: required executable is missing or not executable: $path"; exit 1; }
done
for path in "$PARSER" "$MERGER" "$MAG_SUMMARIZER"; do
    [[ -s "$path" ]] || { echo "ERROR: required Python script is missing or empty: $path"; exit 1; }
done
for tool in samtools python3 awk sha256sum readlink cmp; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing tool: $tool"; exit 1; }
done
[[ -s "$READS" ]] || { echo "ERROR: missing pooled HiFi reads: $READS"; exit 1; }
[[ -s "$MAG_CONTIG_MAP" ]] || { echo "ERROR: MAG_CONTIG_MAP is required: $MAG_CONTIG_MAP"; exit 1; }
[[ -s "$PER_ASSEMBLER_MANIFEST" ]] || { echo "ERROR: missing assembler mapping manifest; rerun 01_mapping.sh: $PER_ASSEMBLER_MANIFEST"; exit 1; }
[[ "$FORCE_RERUN_ASSEMBLY_VARIANTS" =~ ^[YN]$ ]] || { echo "ERROR: FORCE_RERUN_ASSEMBLY_VARIANTS must be Y or N"; exit 1; }
[[ "$ASSEMBLY_VARIANT_KEEP_FULL_PROFILE" =~ ^[YN]$ ]] || { echo "ERROR: ASSEMBLY_VARIANT_KEEP_FULL_PROFILE must be Y or N"; exit 1; }
[[ "$ASSEMBLY_VARIANT_DROP_FILE_CACHE" =~ ^[YN]$ ]] || { echo "ERROR: ASSEMBLY_VARIANT_DROP_FILE_CACHE must be Y or N"; exit 1; }
[[ "$ASSEMBLY_VARIANT_SKIP_GENE_CALLING" =~ ^[YN]$ ]] || { echo "ERROR: ASSEMBLY_VARIANT_SKIP_GENE_CALLING must be Y or N"; exit 1; }
[[ "$ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES" =~ ^[YN]$ ]] || { echo "ERROR: ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES must be Y or N"; exit 1; }
[[ "$ASSEMBLY_VARIANT_SAMPLE_ID" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "ERROR: invalid ASSEMBLY_VARIANT_SAMPLE_ID: $ASSEMBLY_VARIANT_SAMPLE_ID"; exit 1; }
for numeric_var in ANVIO_THREADS ASSEMBLY_VARIANT_MEMORY_BUDGET_GB ASSEMBLY_VARIANT_MEMORY_RESERVE_GB ALLELE_ANVIO_MEMORY_POLL_SECONDS ALLELE_ANVIO_FULL_PROFILE_THREADS ALLELE_ANVIO_PREPROCESS_THREADS; do
    numeric_value="${!numeric_var}"
    [[ "$numeric_value" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: $numeric_var must be a positive integer (got: $numeric_value)"; exit 1; }
done
if (( ALLELE_ANVIO_FULL_PROFILE_THREADS > ANVIO_THREADS )); then
    echo "ERROR: requested anvi-profile threads (${ALLELE_ANVIO_FULL_PROFILE_THREADS}) exceed the detected CPU allocation (${ANVIO_THREADS})."
    echo "       Allocate at least ${ALLELE_ANVIO_FULL_PROFILE_THREADS} CPUs or explicitly lower ALLELE_ANVIO_FULL_PROFILE_THREADS."
    exit 1
fi
(( ALLELE_ANVIO_PREPROCESS_THREADS > ANVIO_THREADS )) && ALLELE_ANVIO_PREPROCESS_THREADS="$ANVIO_THREADS"

MEM_AVAILABLE_GB="$(awk '/^MemAvailable:/ {print int($2/1024/1024); exit}' /proc/meminfo)"
[[ "$MEM_AVAILABLE_GB" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: could not detect MemAvailable from /proc/meminfo"; exit 1; }
if (( MEM_AVAILABLE_GB <= ASSEMBLY_VARIANT_MEMORY_RESERVE_GB )); then
    echo "ERROR: MemAvailable=${MEM_AVAILABLE_GB} GB is not above the reserved ${ASSEMBLY_VARIANT_MEMORY_RESERVE_GB} GB."
    exit 1
fi
PROFILE_MEMORY_EFFECTIVE_GB="$ASSEMBLY_VARIANT_MEMORY_BUDGET_GB"
max_by_available=$(( MEM_AVAILABLE_GB - ASSEMBLY_VARIANT_MEMORY_RESERVE_GB ))
(( PROFILE_MEMORY_EFFECTIVE_GB > max_by_available )) && PROFILE_MEMORY_EFFECTIVE_GB="$max_by_available"

ANVIO_VERSION="$($ANVIO_PYTHON - <<'PY_ANVIO_VERSION'
import anvio
print(anvio.__version__)
PY_ANVIO_VERSION
)"
PARSER_SHA256="$(sha256sum "$PARSER" | awk '{print $1}')"
MERGER_SHA256="$(sha256sum "$MERGER" | awk '{print $1}')"
MAG_SUMMARIZER_SHA256="$(sha256sum "$MAG_SUMMARIZER" | awk '{print $1}')"
MM2_INDEX_CHUNK="${MINIMAP2_INDEX_CHUNK:-256G}"
if [[ "${FORCE_SINGLE_PART_MINIMAP2:-Y}" == "Y" ]]; then
    MM2_INDEX_CHUNK="${MINIMAP2_INDEX_CHUNK_SINGLE:-900G}"
fi
MINIMAP2_MAX_SECONDARY="${MINIMAP2_MAX_SECONDARY:-5}"

run_anvio_reformat_bam() {
    local input_bam="$1"
    local report_file="$2"
    local output_bam="$3"
    local help_text
    help_text="$("$ANVIO_REFORMAT_BAM" --help 2>&1 || true)"

    # anvi'o documents --list/-l for the rename report. Prefer that option
    # explicitly; some help text also mentions --report-file while describing
    # anvi-script-reformat-fasta, which is not sufficient evidence that -r is
    # the rename-list option for anvi-script-reformat-bam.
    if grep -q -- '--list' <<< "$help_text"; then
        "$ANVIO_REFORMAT_BAM" "$input_bam" --list "$report_file" --output "$output_bam"
    elif grep -qE '(^|[[:space:],])-l([=[:space:],]|$)' <<< "$help_text"; then
        "$ANVIO_REFORMAT_BAM" "$input_bam" -l "$report_file" -o "$output_bam"
    elif grep -q -- '--report-file' <<< "$help_text"; then
        "$ANVIO_REFORMAT_BAM" -r "$report_file" -o "$output_bam" "$input_bam"
    else
        echo "ERROR: could not determine anvi-script-reformat-bam rename-list option" >&2
        return 1
    fi
}


validate_anvio_renamed_inputs() {
    local original_fai="$1"
    local renamed_fai="$2"
    local name_map_file="$3"
    local renamed_bam="$4"

    # This is the decisive downstream compatibility test: every BAM reference
    # name and length must match the FASTA that will be used to build CONTIGS.db.
    qc_bam_matches_fai "$renamed_bam" "$renamed_fai" || return 1

    # anvi'o rename reports encountered across versions/workflows can be written
    # in either original->renamed or renamed->original orientation. Original
    # entries may also be full FASTA deflines, while .fai/BAM references use the
    # first whitespace-delimited token. Detect the orientation from the source
    # and renamed FASTA indices instead of assuming a fixed column order.
    python3 - "$original_fai" "$renamed_fai" "$name_map_file" <<'PY_VALIDATE_ANVIO_PREP'
import sys

original_fai, renamed_fai, name_map_path = sys.argv[1:]

def load_fai(path):
    rows = []
    seen = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            name, length = fields[0], int(fields[1])
            if name in seen:
                raise SystemExit(f"duplicate FASTA sequence ID in {path}: {name}")
            seen.add(name)
            rows.append((name, length))
    if not rows:
        raise SystemExit(f"empty FASTA index: {path}")
    return rows

def normalize_original_defline(value):
    value = value.strip()
    if value.startswith(">"):
        value = value[1:]
    return value.split(None, 1)[0] if value else value

original_rows = load_fai(original_fai)
renamed_rows = load_fai(renamed_fai)
if len(original_rows) != len(renamed_rows):
    raise SystemExit(
        f"renamed FASTA sequence count differs from source: "
        f"{len(original_rows)} != {len(renamed_rows)}"
    )

for i, ((old_name, old_len), (new_name, new_len)) in enumerate(
    zip(original_rows, renamed_rows), start=1
):
    if old_len != new_len:
        raise SystemExit(
            f"sequence length changed during anvi'o FASTA rename at row {i}: "
            f"{old_name}={old_len}, {new_name}={new_len}"
        )

pairs = []
with open(name_map_path, encoding="utf-8-sig") as handle:
    for line_no, line in enumerate(handle, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            raise SystemExit(f"invalid anvi'o rename-report row {line_no}")
        first, second = fields[0].strip(), fields[1].strip()
        if not first or not second:
            raise SystemExit(f"blank anvi'o rename-report ID at row {line_no}")
        pairs.append((line_no, first, second))

if len(pairs) != len(renamed_rows):
    raise SystemExit(
        f"anvi'o rename-report row count differs from renamed FASTA: "
        f"{len(pairs)} != {len(renamed_rows)}"
    )

original = {name: length for name, length in original_rows}
renamed = {name: length for name, length in renamed_rows}

def validate_orientation(mode):
    normalized = []
    original_seen = set()
    renamed_seen = set()
    for line_no, first, second in pairs:
        if mode == "original_renamed":
            old = normalize_original_defline(first)
            new = second
        else:
            new = first
            old = normalize_original_defline(second)

        if old not in original or new not in renamed:
            return None
        if old in original_seen or new in renamed_seen:
            return None
        if original[old] != renamed[new]:
            return None
        original_seen.add(old)
        renamed_seen.add(new)
        normalized.append((old, new))

    if original_seen != set(original) or renamed_seen != set(renamed):
        return None
    return normalized

old_new = validate_orientation("original_renamed")
new_old = validate_orientation("renamed_original")

if old_new is None and new_old is None:
    line_no, first, second = pairs[0]
    raise SystemExit(
        "anvi'o rename-report does not match the source and renamed FASTA indices; "
        f"first record at line {line_no}: {first!r} <-> {second!r}"
    )
if old_new is not None and new_old is not None and old_new != new_old:
    raise SystemExit("anvi'o rename-report orientation is ambiguous")

orientation = "original->renamed" if old_new is not None else "renamed->original"
print(f"[40-LR] anvi'o rename-report orientation: {orientation}", file=sys.stderr)
PY_VALIDATE_ANVIO_PREP
}


reheader_bam_to_renamed_fasta() {
    local input_bam="$1"
    local original_fai="$2"
    local renamed_fai="$3"
    local output_bam="$4"
    local header_tmp="${output_bam}.header.tmp.$$"

    # Fallback used only when anvi-script-reformat-bam does not produce a BAM
    # whose dictionary matches assembly.anvio.fa. Because --simplify-names does
    # not remove/reorder sequences here, the original and renamed FAI rows give
    # an unambiguous old->new reference-name map. Reheadering changes only BAM
    # reference names; alignments are not remapped or otherwise modified.
    python3 - "$input_bam" "$original_fai" "$renamed_fai" "$header_tmp" <<'PY_REHEADER_ANVIO_BAM'
import subprocess
import sys

bam, original_fai, renamed_fai, output_header = sys.argv[1:]

def load_fai(path):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                rows.append((fields[0], int(fields[1])))
    if not rows:
        raise SystemExit(f"empty FASTA index: {path}")
    return rows

old_rows = load_fai(original_fai)
new_rows = load_fai(renamed_fai)
if len(old_rows) != len(new_rows):
    raise SystemExit("cannot reheader BAM: source and renamed FASTA counts differ")

rename = {}
old_lengths = {}
new_seen = set()
for i, ((old_name, old_len), (new_name, new_len)) in enumerate(zip(old_rows, new_rows), start=1):
    if old_len != new_len:
        raise SystemExit(
            f"cannot reheader BAM: sequence length differs at row {i}: "
            f"{old_name}={old_len}, {new_name}={new_len}"
        )
    if old_name in rename:
        raise SystemExit(f"cannot reheader BAM: duplicate source ID {old_name}")
    if new_name in new_seen:
        raise SystemExit(f"cannot reheader BAM: duplicate renamed ID {new_name}")
    rename[old_name] = new_name
    old_lengths[old_name] = old_len
    new_seen.add(new_name)

header = subprocess.check_output(["samtools", "view", "-H", bam], text=True)
seen = set()
out = []
for line in header.splitlines():
    if not line.startswith("@SQ\t"):
        out.append(line)
        continue
    fields = line.split("\t")
    sn_i = None
    ln = None
    for i, field in enumerate(fields[1:], start=1):
        if field.startswith("SN:"):
            sn_i = i
        elif field.startswith("LN:"):
            ln = int(field[3:])
    if sn_i is None or ln is None:
        raise SystemExit("cannot reheader BAM: malformed @SQ line")
    old_name = fields[sn_i][3:]
    if old_name not in rename:
        raise SystemExit(f"cannot reheader BAM: BAM reference not in source FASTA: {old_name}")
    old_len = old_lengths[old_name]
    if ln != old_len:
        raise SystemExit(
            f"cannot reheader BAM: BAM/source length mismatch for {old_name}: {ln} != {old_len}"
        )
    fields[sn_i] = "SN:" + rename[old_name]
    seen.add(old_name)
    out.append("\t".join(fields))

if seen != set(rename):
    missing = sorted(set(rename) - seen)
    raise SystemExit(
        f"cannot reheader BAM: {len(missing)} source references are absent from BAM header"
    )

with open(output_header, "w", encoding="utf-8") as handle:
    handle.write("\n".join(out) + "\n")
PY_REHEADER_ANVIO_BAM

    samtools reheader "$header_tmp" "$input_bam" > "$output_bam"
    rm -f "$header_tmp"
}

write_anvio_split_map() {
    local db="$1"
    local output="$2"
    python3 - "$db" "$output" <<'PY_SPLIT_MAP'
import sqlite3
import sys

db_path, out_path = sys.argv[1:]
con = sqlite3.connect(db_path)
try:
    rows = con.execute(
        "SELECT split, parent FROM splits_basic_info ORDER BY parent, order_in_parent, split"
    ).fetchall()
finally:
    con.close()
if not rows:
    raise SystemExit(f"ERROR: no splits found in contigs DB: {db_path}")
with open(out_path, "w", encoding="utf-8", newline="") as out:
    out.write("split_name\tparent_contig\n")
    for split, parent in rows:
        out.write(f"{split}\t{parent}\n")
print(f"split_records={len(rows)}")
PY_SPLIT_MAP
}

count_anvio_splits() {
    local db="$1"
    python3 - "$db" <<'PY_COUNT_SPLITS'
import sqlite3, sys
con = sqlite3.connect(sys.argv[1])
try:
    print(con.execute("SELECT COUNT(*) FROM splits_basic_info").fetchone()[0])
finally:
    con.close()
PY_COUNT_SPLITS
}

process_tree_rss_kb() {
    local root_pid="$1"
    python3 - "$root_pid" <<'PY_TREE_RSS'
import os, sys
root = int(sys.argv[1])
children = {}
for entry in os.scandir('/proc'):
    if not entry.name.isdigit():
        continue
    pid = int(entry.name)
    try:
        with open(f'/proc/{pid}/stat', encoding='utf-8') as h:
            fields = h.read().split()
        ppid = int(fields[3])
    except Exception:
        continue
    children.setdefault(ppid, []).append(pid)
stack=[root]
seen=set()
total=0
while stack:
    pid=stack.pop()
    if pid in seen:
        continue
    seen.add(pid)
    stack.extend(children.get(pid, []))
    try:
        with open(f'/proc/{pid}/status', encoding='utf-8') as h:
            for line in h:
                if line.startswith('VmRSS:'):
                    total += int(line.split()[1])
                    break
    except Exception:
        pass
print(total)
PY_TREE_RSS
}

terminate_process_tree() {
    local root_pid="$1"
    python3 - "$root_pid" <<'PY_KILL_TREE'
import os, signal, sys, time
root = int(sys.argv[1])
def descendants(root):
    children={}
    for entry in os.scandir('/proc'):
        if not entry.name.isdigit():
            continue
        pid=int(entry.name)
        try:
            with open(f'/proc/{pid}/stat', encoding='utf-8') as h:
                fields=h.read().split()
            children.setdefault(int(fields[3]), []).append(pid)
        except Exception:
            pass
    out=[]; stack=[root]
    while stack:
        pid=stack.pop()
        if pid in out:
            continue
        out.append(pid); stack.extend(children.get(pid, []))
    return out
pids=descendants(root)
for pid in reversed(pids):
    try: os.kill(pid, signal.SIGTERM)
    except ProcessLookupError: pass
time.sleep(2)
for pid in reversed(pids):
    try: os.kill(pid, signal.SIGKILL)
    except ProcessLookupError: pass
PY_KILL_TREE
}

install_if_changed() {
    local src="$1" dst="$2"
    if [[ -s "$dst" ]] && cmp -s -- "$src" "$dst"; then
        rm -f -- "$src"
    else
        mv -f -- "$src" "$dst"
    fi
}

drop_file_cache_best_effort() {
    [[ "$ASSEMBLY_VARIANT_DROP_FILE_CACHE" == "Y" ]] || return 0
    python3 - "$@" <<'PY_DROP_CACHE' || true
import os, sys
if not hasattr(os, 'posix_fadvise') or not hasattr(os, 'POSIX_FADV_DONTNEED'):
    raise SystemExit(0)
for path in sys.argv[1:]:
    try:
        fd=os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)
    except Exception:
        pass
PY_DROP_CACHE
}

run_full_anvi_profile_guarded() {
    local input_bam="$1" contigs_db="$2" output_dir="$3" assembler="$4" threads="$5"
    local budget_kb=$(( PROFILE_MEMORY_EFFECTIVE_GB * 1024 * 1024 ))
    local reserve_kb=$(( ASSEMBLY_VARIANT_MEMORY_RESERVE_GB * 1024 * 1024 ))
    local pid rc rss_kb avail_kb peak_kb=0 guard_hit=0

    rm -rf "$output_dir"
    "$ANVIO_PROFILE" -i "$input_bam" -o "$output_dir" -T "$threads" \
        -S "$ASSEMBLY_VARIANT_SAMPLE_ID" -c "$contigs_db" \
        --min-contig-length 1000 \
        --skip-hierarchical-clustering &
    pid=$!

    while kill -0 "$pid" 2>/dev/null; do
        sleep "$ALLELE_ANVIO_MEMORY_POLL_SECONDS"
        kill -0 "$pid" 2>/dev/null || break
        rss_kb="$(process_tree_rss_kb "$pid" 2>/dev/null || echo 0)"
        [[ "$rss_kb" =~ ^[0-9]+$ ]] || rss_kb=0
        (( rss_kb > peak_kb )) && peak_kb="$rss_kb"
        avail_kb="$(awk '/^MemAvailable:/ {print $2; exit}' /proc/meminfo)"
        [[ "$avail_kb" =~ ^[0-9]+$ ]] || avail_kb=0
        if (( rss_kb >= budget_kb || avail_kb <= reserve_kb )); then
            echo "[40-LR][$assembler] Memory guard reached: process_tree_RSS=$((rss_kb/1024/1024)) GB, MemAvailable=$((avail_kb/1024/1024)) GB; stopping this attempt."
            guard_hit=1
            terminate_process_tree "$pid"
            break
        fi
    done

    set +e
    wait "$pid"
    rc=$?
    set -e
    echo "[40-LR][$assembler] anvi-profile attempt peak monitored RSS: $((peak_kb/1024/1024)) GB (threads=$threads)."
    (( guard_hit == 0 )) || return 99
    return "$rc"
}

# Read and validate the selected full-co-assembly HiFi alignment population.
MANIFEST_RECORDS="$(
python3 - "$PER_ASSEMBLER_MANIFEST" "$ASSEMBLY_VARIANT_HIFI_BAM_COLUMN" <<'PY_READ_LR_MANIFEST'
import csv, re, sys
path, bam_column = sys.argv[1:]
required = {"Assembler", "Assembly_FASTA", "HiFi_all_BAM", bam_column}
rows, seen = [], set()
with open(path, encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    missing = sorted(required - set(reader.fieldnames or []))
    if missing:
        raise SystemExit("ERROR: assembler mapping manifest lacks required LR columns: " + ", ".join(missing))
    for line_no, raw in enumerate(reader, 2):
        assembler = (raw.get("Assembler") or "").strip()
        assembly = (raw.get("Assembly_FASTA") or "").strip()
        hifi_all = (raw.get("HiFi_all_BAM") or "").strip()
        hifi = (raw.get(bam_column) or "").strip()
        if not assembler or not re.fullmatch(r"[A-Za-z0-9._-]+", assembler):
            raise SystemExit(f"ERROR: invalid assembler name at manifest line {line_no}: {assembler!r}")
        if assembler in seen:
            raise SystemExit(f"ERROR: duplicate assembler in mapping manifest: {assembler}")
        seen.add(assembler)
        if not assembly or assembly in {"NA", "None"}:
            raise SystemExit(f"ERROR: Assembly_FASTA is unavailable for {assembler}")
        if not hifi_all or hifi_all in {"NA", "None"}:
            raise SystemExit(f"ERROR: HiFi_all_BAM is unavailable for {assembler}; rerun 01_mapping.sh")
        if not hifi or hifi in {"NA", "None"}:
            raise SystemExit(f"ERROR: {bam_column} is unavailable for {assembler}; rerun 01_mapping.sh")
        rows.append((assembler, assembly, hifi_all, hifi))
if not rows:
    raise SystemExit("ERROR: assembler mapping manifest contains no data rows")
for row in rows:
    print("\t".join(row))
PY_READ_LR_MANIFEST
)"

if [[ "$FORCE_RERUN_ASSEMBLY_VARIANTS" == "Y" ]]; then
    echo "[40-LR] FORCE_RERUN_ASSEMBLY_VARIANTS=Y; removing cached assembler-specific anvi'o results."
    rm -rf "$ALLELE_BY_ASSEMBLER_DIR"
    mkdir -p "$ALLELE_BY_ASSEMBLER_DIR"
    rm -f "$ALLELE_PROFILE_MANIFEST" "$ALLELE_VARIANTS" "$ALLELE_CONTIG_ALL" \
          "$ALLELE_SUMMARY_ALL" "$ALLELE_ASSEMBLER_SUMMARY" "$ALLELE_CONTIG" \
          "$ALLELE_SUMMARY" "$ALLELE_VARIANTS_MAG" "$ALLELE_MAG" "$ALLELE_CACHE_META"
fi

PROFILE_MANIFEST_TMP="${ALLELE_PROFILE_MANIFEST}.tmp.$$"
printf 'Assembler\tAssembly_FASTA\tHiFi_all_BAM\tHiFi_profile_BAM\tAlignment_population\tSNV_raw_TSV\tINDEL_raw_TSV\tDetail_TSV\tContig_TSV\tSummary_TSV\n' > "$PROFILE_MANIFEST_TMP"

GLOBAL_FP_PARTS=()

echo "[09_assembly_variants] Started: $(date)"
echo "SCRIPT_ID: $SCRIPT_ID"
echo "anvi'o version: $ANVIO_VERSION"
echo "HiFi scope: every >=1,000-bp contig in each corresponding full co-assembly; assembly-allele classification is site-level from anvi'o-reported variability"
echo "HiFi reference population: each assembler is profiled separately"
echo "HiFi alignment population: $ASSEMBLY_VARIANT_ALIGNMENT_POPULATION ($ASSEMBLY_VARIANT_HIFI_BAM_COLUMN)"
echo "MAG outputs: subset from the assembly-wide classifications; no second HiFi mapping/profile"
echo "anvi-profile execution: one complete assembler-specific co-assembly at a time; no contig batching; assemblers serial"
echo "Memory planning: configured_budget=${ASSEMBLY_VARIANT_MEMORY_BUDGET_GB} GB; MemAvailable_at_start=${MEM_AVAILABLE_GB} GB; host_reserve=${ASSEMBLY_VARIANT_MEMORY_RESERVE_GB} GB; effective_budget=${PROFILE_MEMORY_EFFECTIVE_GB} GB"
echo "anvi-profile threads: selected=${ALLELE_ANVIO_FULL_PROFILE_THREADS}; CPU_allocation=${ANVIO_THREADS}; fixed default=100"
echo "Runtime memory guard: poll=${ALLELE_ANVIO_MEMORY_POLL_SECONDS}s; retry with fewer threads if the effective budget is reached"
echo "anvi'o preprocessing/index threads: $ALLELE_ANVIO_PREPROCESS_THREADS"
echo "Full PROFILE.db retention: $ASSEMBLY_VARIANT_KEEP_FULL_PROFILE"

while IFS=$'\t' read -r assembler assembly hifi_all_bam hifi_bam; do
    echo "[40-LR][$assembler] Validating corresponding full co-assembly inputs."
    [[ -s "$assembly" ]] || { echo "ERROR: missing corresponding assembly for $assembler: $assembly"; exit 1; }
    [[ -s "$hifi_bam" ]] || { echo "ERROR: missing corresponding HiFi BAM for $assembler: $hifi_bam"; exit 1; }
    samtools faidx "$assembly"
    samtools quickcheck -v "$hifi_bam" || { echo "ERROR: HiFi BAM failed quickcheck: $hifi_bam"; exit 1; }
    qc_bam_matches_fai "$hifi_bam" "${assembly}.fai" || { echo "ERROR: HiFi BAM/reference mismatch for $assembler"; exit 1; }
    if [[ ! -s "${hifi_bam}.bai" && ! -s "${hifi_bam}.csi" ]]; then
        samtools index -@ "$ALLELE_ANVIO_PREPROCESS_THREADS" -c "$hifi_bam"
    fi

    hifi_reference_fp="$(qc_compute_mapping_fingerprint "$assembly" "$READS" "$PRESET" "$MM2_INDEX_CHUNK")"
    hifi_all_fp="$(qc_hash_text \
        "$hifi_reference_fp" \
        "assembler=$assembler" \
        "minimap2=-ax_${PRESET}_-p_1_-N_${MINIMAP2_MAX_SECONDARY}_--secondary-seq" \
        "retained=mapped_records_only" \
        "sort=coordinate")"
    expected_hifi_fp="$hifi_all_fp"
    qc_cache_matches "${hifi_bam}.inputs.sha256" "$expected_hifi_fp" || {
        echo "ERROR: stale or missing $ASSEMBLY_VARIANT_HIFI_BAM_COLUMN cache for $assembler; rerun 01_mapping.sh"
        exit 1
    }

    n_unmapped="$(samtools view -c -f 4 "$hifi_bam")"
    (( n_unmapped == 0 )) || { echo "ERROR: unmapped records remain in $hifi_bam"; exit 1; }

    adir="$ALLELE_BY_ASSEMBLER_DIR/$assembler"
    mkdir -p "$adir"
    reformatted_fa="$adir/assembly.anvio.fa"
    name_map="$adir/assembly.anvio.name_map.tsv"
    reformatted_bam="$adir/HiFi_all.anvio.sorted.bam"
    contigs_db="$adir/CONTIGS.db"
    contigs_db_cache="$adir/CONTIGS.inputs.sha256"
    profile_dir="$adir/PROFILE"
    profile_db="$profile_dir/PROFILE.db"
    profile_batches="$adir/PROFILE_BATCHES"  # legacy batch-local directory; removed after successful full profiling
    split_map="$adir/anvio_splits_to_contigs.tsv"
    snv_raw="$adir/variable_nucleotides.anvio.tsv"
    indel_raw="$adir/indels.anvio.tsv"
    detail_tsv="$adir/assembly_allele_disagreements.tsv"
    contig_tsv="$adir/contig_assembly_allele_support.tsv"
    summary_tsv="$adir/assembly_allele_support.summary.tsv"
    preprocess_cache="$adir/preprocess.inputs.sha256"
    profile_cache="$adir/profile.inputs.sha256"
    classification_cache="$adir/classification.inputs.sha256"

    preprocess_fp="$(python3 - "$assembly" "$hifi_bam" "$expected_hifi_fp" "$ANVIO_VERSION" <<'PY_PREP_FP'
import hashlib, json, os, sys
assembly, bam, mapping_fp, anvio_version = sys.argv[1:]
def st(p):
    s=os.stat(p); return [os.path.realpath(p), s.st_size, s.st_mtime_ns]
payload={
    "workflow":"assembly_variant_anvio_preprocess_name_map",
    "assembly":st(assembly), "bam":st(bam),
    "mapping_fingerprint":mapping_fp, "anvio_version":anvio_version,
    "reformat_fasta":"simplify_names", "reformat_bam":"rename_dictionary_no_remap",
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_PREP_FP
)"

    if [[ -s "$preprocess_cache" ]] && ! qc_cache_matches "$preprocess_cache" "$preprocess_fp"; then
        echo "[40-LR][$assembler] Source assembly/BAM changed; rebuilding renamed anvi'o FASTA/BAM and full-reference profile."
        rm -rf "$profile_dir" "$profile_batches"
        rm -f "$reformatted_fa" "${reformatted_fa}.fai" "$name_map" \
              "$reformatted_bam" "${reformatted_bam}.bai" "${reformatted_bam}.csi" \
              "$contigs_db" "$contigs_db_cache" "$split_map" "$snv_raw" "$indel_raw" \
              "$detail_tsv" "$contig_tsv" "$summary_tsv" \
              "$preprocess_cache" "$profile_cache" "$classification_cache"
    fi

    # Reuse expensive renamed FASTA/BAM intermediates from earlier or failed runs.
    if [[ ! -s "$reformatted_fa" || ! -s "$name_map" ]]; then
        echo "[40-LR][$assembler] Reformatting the corresponding full co-assembly FASTA for anvi'o."
        "$ANVIO_REFORMAT_FASTA" "$assembly" --simplify-names --report-file "$name_map.tmp.$$" -o "$reformatted_fa.tmp.$$"
        mv -f "$reformatted_fa.tmp.$$" "$reformatted_fa"
        mv -f "$name_map.tmp.$$" "$name_map"
    fi
    samtools faidx "$reformatted_fa"

    # A previous interrupted/failed run may have left a non-empty but invalid
    # renamed BAM. Validate it before reuse and rebuild automatically if needed.
    if [[ -s "$reformatted_bam" ]] && ! qc_bam_matches_fai "$reformatted_bam" "${reformatted_fa}.fai"; then
        echo "[40-LR][$assembler] Existing renamed anvi'o BAM is inconsistent; rebuilding it automatically."
        rm -f "$reformatted_bam" "${reformatted_bam}.bai" "${reformatted_bam}.csi"
    fi

    if [[ ! -s "$reformatted_bam" ]]; then
        renamed_tmp="$reformatted_bam.tmp.$$.bam"
        rm -f "$renamed_tmp" "${renamed_tmp}.bai" "${renamed_tmp}.csi"
        echo "[40-LR][$assembler] Renaming the existing full-co-assembly HiFi BAM dictionary; reads are not remapped."

        if run_anvio_reformat_bam "$hifi_bam" "$name_map" "$renamed_tmp" \
           && qc_bam_matches_fai "$renamed_tmp" "${reformatted_fa}.fai"; then
            mv -f "$renamed_tmp" "$reformatted_bam"
        else
            echo "[40-LR][$assembler][WARN] anvi-script-reformat-bam output did not match the renamed FASTA; using a reference-order samtools reheader fallback."
            rm -f "$renamed_tmp"
            reheader_bam_to_renamed_fasta \
                "$hifi_bam" "${assembly}.fai" "${reformatted_fa}.fai" "$renamed_tmp"
            qc_bam_matches_fai "$renamed_tmp" "${reformatted_fa}.fai" || {
                rm -f "$renamed_tmp"
                echo "ERROR: both anvi'o BAM rename and samtools reheader fallback failed for $assembler"
                exit 1
            }
            mv -f "$renamed_tmp" "$reformatted_bam"
        fi
        samtools index -@ "$ALLELE_ANVIO_PREPROCESS_THREADS" -c "$reformatted_bam"
    elif [[ ! -s "${reformatted_bam}.bai" && ! -s "${reformatted_bam}.csi" ]]; then
        samtools index -@ "$ALLELE_ANVIO_PREPROCESS_THREADS" -c "$reformatted_bam"
    fi

    validate_anvio_renamed_inputs "${assembly}.fai" "${reformatted_fa}.fai" "$name_map" "$reformatted_bam" || {
        echo "ERROR: renamed anvi'o FASTA/BAM are genuinely inconsistent for $assembler after automatic repair"
        exit 1
    }
    qc_write_cache_atomic "$preprocess_cache" "$preprocess_fp"

    contigs_db_fp="$(python3 - "$preprocess_fp" "$ANVIO_VERSION" "$ASSEMBLY_VARIANT_SKIP_GENE_CALLING" <<'PY_CONTIGS_DB_FP'
import hashlib, json, sys
prep_fp, anvio_version, skip_gene_calling = sys.argv[1:]
payload={
    "workflow":"assembly_variant_full_contigs_db_v1",
    "preprocess":prep_fp,
    "anvio_version":anvio_version,
    "skip_gene_calling":skip_gene_calling,
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_CONTIGS_DB_FP
)"
    if [[ -s "$contigs_db_cache" ]] && ! qc_cache_matches "$contigs_db_cache" "$contigs_db_fp"; then
        rm -f "$contigs_db" "$contigs_db_cache" "$split_map"
        rm -rf "$profile_dir"
    fi
    if [[ ! -s "$contigs_db" ]]; then
        echo "[40-LR][$assembler] Building one full co-assembly CONTIGS.db."
        rm -f "$contigs_db.tmp.$$.db"
        contigs_db_args=( -f "$reformatted_fa" -o "$contigs_db.tmp.$$.db" -T "$ALLELE_ANVIO_PREPROCESS_THREADS" )
        if [[ "$ASSEMBLY_VARIANT_SKIP_GENE_CALLING" == "Y" ]]; then
            contigs_db_args+=( --skip-gene-calling )
        fi
        "$ANVIO_GEN_CONTIGS_DB" "${contigs_db_args[@]}"
        [[ -s "$contigs_db.tmp.$$.db" ]] || { echo "ERROR: full CONTIGS.db was not created for $assembler"; exit 1; }
        mv -f "$contigs_db.tmp.$$.db" "$contigs_db"
        qc_write_cache_atomic "$contigs_db_cache" "$contigs_db_fp"
    elif [[ ! -s "$contigs_db_cache" ]]; then
        # Existing full DB from an earlier compatible workflow is adopted only
        # after validating that it contains splits.
        split_probe="$(count_anvio_splits "$contigs_db")"
        [[ "$split_probe" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: existing CONTIGS.db is invalid for $assembler"; exit 1; }
        qc_write_cache_atomic "$contigs_db_cache" "$contigs_db_fp"
    fi

    if [[ ! -s "$split_map" ]]; then
        write_anvio_split_map "$contigs_db" "$split_map.tmp.$$"
        mv -f "$split_map.tmp.$$" "$split_map"
    fi
    split_count="$(count_anvio_splits "$contigs_db")"
    echo "[40-LR][$assembler] Full-reference CONTIGS.db splits: $split_count."
    if (( split_count > 200000 && ALLELE_ANVIO_FULL_PROFILE_THREADS > 20 )); then
        echo "[40-LR][$assembler] NOTE: >200,000 splits with >20 threads can be memory-intensive in anvi-profile; runtime memory guard is enabled."
    fi

    profile_fp="$(python3 - "$preprocess_fp" "$contigs_db_fp" "$ANVIO_VERSION" "$ASSEMBLY_VARIANT_SAMPLE_ID" "$ASSEMBLY_VARIANT_PROFILE_ALIGNMENT_CLASSES" <<'PY_PROFILE_FP'
import hashlib, json, sys
prep_fp, db_fp, anvio_version, sample_id, alignment_classes = sys.argv[1:]
payload={
    "workflow":"assembly_variant_anvi_profile_assembler_reference",
    "preprocess":prep_fp,
    "contigs_db":db_fp,
    "anvio_version":anvio_version,
    "sample_id":sample_id,
    "min_contig_length":1000,
    "minimum_coverage_to_calculate_snvs":10,
    "profile_snvs":True,
    "profile_indels":True,
    "full_variability":False,
    "percent_identity_filter":None,
    "fetch_filter":None,
    "alignment_classes":alignment_classes,
    "mapping_reference":"existing_corresponding_full_coassembly_competitive_bam",
    "profiling_reference":"complete_corresponding_assembler_coassembly_single_profile",
    "contig_batching":False,
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_PROFILE_FP
)"

    rm -rf "$adir"/PROFILE.tmp.*
    if ! qc_cache_matches "$profile_cache" "$profile_fp"; then
        rm -f "$snv_raw" "$indel_raw" "$detail_tsv" "$contig_tsv" "$summary_tsv" \
              "$profile_cache" "$classification_cache"
        # A legacy completed full-reference PROFILE.db can still be reused.  It
        # is never confused with PROFILE_BATCHES because the locations differ.
        if [[ -s "$profile_db" ]]; then
            echo "[40-LR][$assembler] Reusing existing full-reference PROFILE.db and exporting raw tables."
        else
            echo "[40-LR][$assembler] Running one full-reference anvi-profile; starting threads=$ALLELE_ANVIO_FULL_PROFILE_THREADS, effective memory budget=${PROFILE_MEMORY_EFFECTIVE_GB} GB."
            attempt_threads="$ALLELE_ANVIO_FULL_PROFILE_THREADS"
            while :; do
                tmp_profile="$adir/PROFILE.tmp.$$"
                if run_full_anvi_profile_guarded "$reformatted_bam" "$contigs_db" "$tmp_profile" "$assembler" "$attempt_threads"; then
                    [[ -s "$tmp_profile/PROFILE.db" ]] || { echo "ERROR: anvi-profile finished without PROFILE.db for $assembler"; exit 1; }
                    rm -rf "$profile_dir"
                    mv "$tmp_profile" "$profile_dir"
                    break
                else
                    rc=$?
                    rm -rf "$tmp_profile"
                    if (( rc == 99 || rc == 137 )); then
                        next_threads=$(( attempt_threads * 3 / 4 ))
                        (( next_threads < 1 )) && next_threads=1
                        if (( next_threads >= attempt_threads )); then
                            next_threads=$(( attempt_threads - 1 ))
                        fi
                        if (( next_threads < 1 )); then
                            echo "ERROR: anvi-profile exceeded the memory budget even at the minimum thread count for $assembler."
                            exit 1
                        fi
                        echo "[40-LR][$assembler] Retrying full-reference anvi-profile with $next_threads threads after memory/OOM stop."
                        attempt_threads="$next_threads"
                    else
                        echo "ERROR: anvi-profile failed for $assembler with exit code $rc."
                        exit "$rc"
                    fi
                fi
            done
        fi

        "$ANVIO_EXPORT_TABLE" --table variable_nucleotides "$profile_db" -o "$snv_raw.tmp.$$"
        "$ANVIO_EXPORT_TABLE" --table indels "$profile_db" -o "$indel_raw.tmp.$$"
        mv -f "$snv_raw.tmp.$$" "$snv_raw"
        mv -f "$indel_raw.tmp.$$" "$indel_raw"
        [[ -s "$snv_raw" && -s "$indel_raw" && -s "$split_map" ]] || {
            echo "ERROR: full-reference anvi'o outputs are missing for $assembler"
            exit 1
        }
        qc_write_cache_atomic "$profile_cache" "$profile_fp"

        if [[ "$ASSEMBLY_VARIANT_KEEP_FULL_PROFILE" == "N" ]]; then
            rm -rf "$profile_dir"
        fi
        if [[ "$ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES" == "Y" && -d "$profile_batches" ]]; then
            echo "[40-LR][$assembler] Removing obsolete PROFILE_BATCHES after successful full-reference profiling."
            rm -rf "$profile_batches"
        fi
        drop_file_cache_best_effort "$assembly" "$hifi_bam" "$reformatted_fa" "$reformatted_bam" "$contigs_db"
    fi

    classification_fp="$(python3 - "$profile_fp" "$PARSER_SHA256" "$MAG_CONTIG_MAP" "$split_map" <<'PY_CLASS_FP'
import hashlib, json, os, sys
profile_fp, parser_sha, mag_map, split_map = sys.argv[1:]
def st(path):
    s=os.stat(path)
    return [os.path.realpath(path), s.st_size, s.st_mtime_ns]
payload={"profile":profile_fp,"parser":parser_sha,"mag_map":st(mag_map),"split_map":st(split_map),
         "criteria":"Trigodet_SNV_dominant_alt_fraction_gt_0.5;_INDEL_event_fraction_gt_0.5;_unsupported_equals_all_coverage"}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_CLASS_FP
)"
    if ! qc_cache_matches "$classification_cache" "$classification_fp"; then
        rm -f "$detail_tsv" "$contig_tsv" "$summary_tsv" "$classification_cache"
    fi
    if [[ ! -s "$detail_tsv" || ! -s "$contig_tsv" || ! -s "$summary_tsv" ]]; then
        echo "[40-LR][$assembler] Classifying minority and unsupported assembly alleles at anvi'o-reported variable sites."
        python3 "$PARSER" \
            --fai "${assembly}.fai" \
            --split-map "$split_map" \
            --snv-table "$snv_raw" \
            --indel-table "$indel_raw" \
            --mag-map "$MAG_CONTIG_MAP" \
            --name-map "$name_map" \
            --renamed-fai "${reformatted_fa}.fai" \
            --detail-out "$detail_tsv.tmp.$$" \
            --contig-out "$contig_tsv.tmp.$$" \
            --summary-out "$summary_tsv.tmp.$$"
        mv -f "$detail_tsv.tmp.$$" "$detail_tsv"
        mv -f "$contig_tsv.tmp.$$" "$contig_tsv"
        mv -f "$summary_tsv.tmp.$$" "$summary_tsv"
        qc_write_cache_atomic "$classification_cache" "$classification_fp"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$assembler" "$assembly" "$hifi_all_bam" "$hifi_bam" "$ASSEMBLY_VARIANT_ALIGNMENT_POPULATION" "$snv_raw" "$indel_raw" \
        "$detail_tsv" "$contig_tsv" "$summary_tsv" >> "$PROFILE_MANIFEST_TMP"
    GLOBAL_FP_PARTS+=("$assembler|$classification_fp")
done <<< "$MANIFEST_RECORDS"

install_if_changed "$PROFILE_MANIFEST_TMP" "$ALLELE_PROFILE_MANIFEST"

echo "HiFi BAM validation: mapped-only BAM confirmed for every assembler; secondary alignments are retained when reported."

echo "[40-LR] Merging assembler-specific classifications."
python3 "$MERGER" \
    --manifest "$ALLELE_PROFILE_MANIFEST" \
    --detail-out "$ALLELE_VARIANTS.tmp.$$" \
    --contig-out "$ALLELE_CONTIG_ALL.tmp.$$" \
    --summary-out "$ALLELE_SUMMARY_ALL.tmp.$$" \
    --assembler-summary-out "$ALLELE_ASSEMBLER_SUMMARY.tmp.$$"
install_if_changed "$ALLELE_VARIANTS.tmp.$$" "$ALLELE_VARIANTS"
mv -f "$ALLELE_CONTIG_ALL.tmp.$$" "$ALLELE_CONTIG_ALL"
mv -f "$ALLELE_SUMMARY_ALL.tmp.$$" "$ALLELE_SUMMARY_ALL"
mv -f "$ALLELE_ASSEMBLER_SUMMARY.tmp.$$" "$ALLELE_ASSEMBLER_SUMMARY"

echo "[40-LR] Subsetting assembly-wide results to final long-read MAG contigs."
python3 "$MAG_SUMMARIZER" \
    --contig-summary "$ALLELE_CONTIG_ALL" \
    --detail "$ALLELE_VARIANTS" \
    --mag-contig-map "$MAG_CONTIG_MAP" \
    --site-output "$ALLELE_VARIANTS_MAG.tmp.$$" \
    --summary-output "$ALLELE_MAG.tmp.$$" \
    --contig-output "$ALLELE_CONTIG.tmp.$$" \
    --overall-summary-output "$ALLELE_SUMMARY.tmp.$$"
install_if_changed "$ALLELE_VARIANTS_MAG.tmp.$$" "$ALLELE_VARIANTS_MAG"
mv -f "$ALLELE_MAG.tmp.$$" "$ALLELE_MAG"
mv -f "$ALLELE_CONTIG.tmp.$$" "$ALLELE_CONTIG"
mv -f "$ALLELE_SUMMARY.tmp.$$" "$ALLELE_SUMMARY"

ALLELE_FINGERPRINT="$(qc_hash_text "$SCRIPT_ID" "$PARSER_SHA256" "$MERGER_SHA256" "$MAG_SUMMARIZER_SHA256" "${GLOBAL_FP_PARTS[@]}")"
qc_write_cache_atomic "$ALLELE_CACHE_META" "$ALLELE_FINGERPRINT"

echo "[40-LR] Outputs:"
echo "  Per-assembler profile manifest: $ALLELE_PROFILE_MANIFEST"
echo "  Assembly-wide disagreement loci: $ALLELE_VARIANTS"
echo "  Assembly-wide by contig:          $ALLELE_CONTIG_ALL"
echo "  Assembly-wide summary:            $ALLELE_SUMMARY_ALL"
echo "  By assembler:                     $ALLELE_ASSEMBLER_SUMMARY"
echo "  Final-MAG disagreement loci:      $ALLELE_VARIANTS_MAG"
echo "  Final-MAG by contig:              $ALLELE_CONTIG"
echo "  Final-MAG by MAG:                 $ALLELE_MAG"
echo "  Final-MAG summary:                $ALLELE_SUMMARY"
echo "[40-LR] HiFi assembly-variant stage complete: $(date)"

echo "[09_assembly_variants] Finished: $(date)"
