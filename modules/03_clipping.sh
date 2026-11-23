#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_ID="03_clipping"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found on PATH" >&2; exit 1; }
python3 "$CLIP_VALIDATE_SCRIPT" check \
    --assembly "$ASSEMBLY_FA" \
    --mag-dir "$MAG_DIR" \
    --manifest "$ASSEMBLER_MAPPING_MANIFEST"

CLIP_CACHE_DIR="${CLIPPING_DIR}/.cache"
CLIP_ASSEMBLER_DIR="${CLIP_CACHE_DIR}/by_assembler"
CLIP_MANIFEST_COPY="${CLIP_CACHE_DIR}/mapping_manifest.tsv"

if [[ -d "${CLIPPING_DIR}/by_assembler" && ! -e "$CLIP_ASSEMBLER_DIR" ]]; then
    mkdir -p "$CLIP_CACHE_DIR"
    mv "${CLIPPING_DIR}/by_assembler" "$CLIP_ASSEMBLER_DIR"
fi
LEGACY_CLIPPING_FILES=(
    clipping.all_regions.default.tsv clipping.internal.default.tsv clipping.terminal.default.tsv
    clipping.all_regions.stringent.tsv clipping.internal.stringent.tsv clipping.terminal.stringent.tsv
    clipping.site.all_regions.tsv clipping.internal.filtered.tsv clipping.terminal.filtered.tsv
    clipping.all_regions.Y.tsv clipping.internal.Y.tsv clipping.terminal.Y.tsv
    clipping.region_density.tsv clipping.region_density_ratio.tsv clipping.mapping_manifest.tsv
    clipping.default.tsv clipping.stringent.tsv clipping.Y.tsv clipping.criteria.tsv
    MAG_clipping_sites.final.tsv MAG_clipping_sites.Y.tsv MAG_clipping.Y.tsv
    anvio_find_misassemblies-clipping.txt anvio_find_misassemblies-zero_cov.txt
    clipping.filtered.tsv clipping.all_regions.filtered.tsv
    .cache/clipping.site.tsv.inputs.sha256
)

mkdir -p "$OUT/00_log" "$CLIP_ASSEMBLER_DIR"

LOG="$OUT/00_log/03_clipping.log"
exec > >(tee -a "$LOG") 2>&1

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

ensure_bam_index() {
    local bam="$1"
    if [[ ! -s "$bam.bai" && ! -s "$bam.csi" ]]; then
        samtools index -@ "${SAMTOOLS_THREADS:-$THREADS}" -c "$bam" \
            || samtools index -@ "${SAMTOOLS_THREADS:-$THREADS}" "$bam"
    fi
}

copy_atomic() {
    local src="$1" dst="$2" tmp="${2}.tmp.$$"
    cp -f -- "$src" "$tmp"
    mv -f -- "$tmp" "$dst"
}


write_mag_outputs() {
    python3 "$CLIP_SUMMARIZE_SCRIPT" mag \
        --clipping "$CLIP_ALL_FINAL" \
        --mag-contig-map "$MAG_CONTIG_MAP" \
        --fai "${ASSEMBLY_FA}.fai" \
        --site-output "$CLIP_ALL_FINAL" \
        --summary-output "$CLIP_MAG_SUMMARY" \
        --contig-summary-output "$CLIP_CONTIG_SUMMARY"
}

print_outputs() {
    echo "Outputs:"
    for path in "$CLIP_ALL_FINAL" "$CLIP_CONTIG_SUMMARY" "$CLIP_MAG_SUMMARY" "$CLIP_ASSEMBLER_SUMMARY"; do
        echo "  $path"
    done
}

for mode in FORCE_RERUN_CLIPPING CLIPPING_POSTPROCESS_ONLY CLIPPING_MAG_SUMMARY_ONLY; do
    value="${!mode}"
    [[ "$value" == "Y" || "$value" == "N" ]] || fail "$mode must be Y or N"
done
[[ "$ANVIO_DEFAULT_CLIPPING_RATIO" == "0.80" || "$ANVIO_DEFAULT_CLIPPING_RATIO" == "0.8" ]] \
    || fail "ANVIO_DEFAULT_CLIPPING_RATIO must remain the anvi'o default 0.80"
[[ "$MIN_DIST_TO_END" == "0" ]] \
    || fail "MIN_DIST_TO_END must be 0 for complete-contig terminal clipping analysis"
[[ "$TERMINAL_CLIPPING_WINDOW_BP" =~ ^[1-9][0-9]*$ ]] \
    || fail "TERMINAL_CLIPPING_WINDOW_BP must be a positive integer"
[[ "$STRINGENT_CLIP_MIN_DEPTH" == "10" ]] \
    || fail "STRINGENT_CLIP_MIN_DEPTH must remain 10"
[[ "$STRINGENT_CLIP_RATIO" == "1.0" || "$STRINGENT_CLIP_RATIO" == "1" ]] \
    || fail "STRINGENT_CLIP_RATIO must remain 1.0"

for tool in samtools python3 sha256sum; do
    command -v "$tool" >/dev/null 2>&1 || fail "missing tool: $tool"
done

[[ -s "$ASSEMBLY_FA" ]] || fail "required input missing or empty: $ASSEMBLY_FA"
[[ -s "$MAG_MAP_BUILDER" ]] || fail "MAG-map builder missing: $MAG_MAP_BUILDER"

[[ -s "${ASSEMBLY_FA}.fai" && "${ASSEMBLY_FA}.fai" -nt "$ASSEMBLY_FA" ]] || samtools faidx "$ASSEMBLY_FA"
mkdir -p "$(dirname "$MAG_CONTIG_MAP")"
python3 "$MAG_MAP_BUILDER" \
    --mag-dir "$MAG_DIR" \
    --extensions "$MAG_EXTENSIONS" \
    --assembly-fai "${ASSEMBLY_FA}.fai" \
    --output "${MAG_CONTIG_MAP}.new.$$"
if [[ -s "$MAG_CONTIG_MAP" ]] && cmp -s "${MAG_CONTIG_MAP}.new.$$" "$MAG_CONTIG_MAP"; then
    rm -f "${MAG_CONTIG_MAP}.new.$$"
else
    mv -f "${MAG_CONTIG_MAP}.new.$$" "$MAG_CONTIG_MAP"
fi

if [[ ! -s "$ASSEMBLER_MAPPING_MANIFEST" ]]; then
    fail "assembler mapping manifest is missing: $ASSEMBLER_MAPPING_MANIFEST. Run 01_mapping.sh before 03_clipping.sh."
fi

for required in \
    "$MAG_CONTIG_MAP" "$ASSEMBLER_MAPPING_MANIFEST" \
    "$CLIP_VALIDATE_SCRIPT" "$CLIP_SUMMARIZE_SCRIPT"; do
    [[ -s "$required" ]] || fail "required input missing or empty: $required"
done
[[ -x "$ANVIO_FIND_MISASSEMBLIES" ]] \
    || fail "anvi-script-find-misassemblies is not executable: $ANVIO_FIND_MISASSEMBLIES"

if [[ "$CLIPPING_MAG_SUMMARY_ONLY" == "Y" ]]; then
    qc_validate_current_clipping_result || exit 1
    [[ -s "$CLIP_ALL_FINAL" ]] || fail "all-region clipping result is required for MAG summary: $CLIP_ALL_FINAL"
    echo "[03_clipping] Regenerating intuitive ALL-region MAG clipping summary."
    write_mag_outputs
    print_outputs
    exit 0
fi

python3 - "$ASSEMBLER_MAPPING_MANIFEST" <<'PY_MANIFEST_SCHEMA'
import csv
import sys
path = sys.argv[1]
with open(path, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'Assembler', 'Assembly_FASTA', 'HiFi_all_BAM'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'ERROR: assembler mapping manifest lacks columns: {sorted(missing)}')
    rows = list(reader)
if not rows:
    raise SystemExit('ERROR: assembler mapping manifest contains no rows')
for i, row in enumerate(rows, start=2):
    assembler = (row.get('Assembler') or '').strip()
    if not assembler or not all(c.isalnum() or c in '._-' for c in assembler):
        raise SystemExit(f'ERROR: invalid assembler label at row {i}: {assembler!r}')
PY_MANIFEST_SCHEMA

copy_atomic "$ASSEMBLER_MAPPING_MANIFEST" "$CLIP_MANIFEST_COPY"

while IFS=$'\t' read -r assembler assembly hifi_all; do
    [[ -s "$assembly" ]] || fail "assembly missing: $assembly"
    samtools faidx "$assembly"
    [[ -s "$hifi_all" && -s "$(qc_bam_cache "$hifi_all")" ]] \
        || fail "HiFi mapping missing for $assembler; rerun 01_mapping.sh"
    samtools quickcheck -v "$hifi_all" \
        || fail "samtools quickcheck failed: $hifi_all"
    qc_bam_matches_fai "$hifi_all" "${assembly}.fai" \
        || fail "HiFi BAM/assembly mismatch: $assembler"
    ensure_bam_index "$hifi_all"
done < <(qc_manifest_columns "$ASSEMBLER_MAPPING_MANIFEST" Assembler Assembly_FASTA HiFi_all_BAM)

python3 - "$ASSEMBLER_MAPPING_MANIFEST" "${ASSEMBLY_FA}.fai" <<'PY_SOURCE_UNION'
import csv
import sys
manifest, assemblies_fai = sys.argv[1:]

def fai(path):
    out = {}
    with open(path, encoding='utf-8') as handle:
        for raw in handle:
            fields = raw.rstrip('\n\r').split('\t')
            if len(fields) < 2:
                raise SystemExit(f'ERROR: malformed FAI row: {path}')
            name, length = fields[0], int(fields[1])
            if name in out:
                raise SystemExit(f'ERROR: duplicate contig ID in FAI: {name}')
            out[name] = length
    return out

combined = fai(assemblies_fai)
union = {}
with open(manifest, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        assembler = row['Assembler'].strip()
        source = fai(row['Assembly_FASTA'].strip() + '.fai')
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
for contig, length in combined.items():
    if union[contig][1] != length:
        raise SystemExit(f'ERROR: assembler/all-assemblies contig length mismatch: {contig}')
print(f'Validated assembler contig union: {len(combined)} contigs')
PY_SOURCE_UNION

if [[ -n "$ANVIO_ENV_PREFIX" && -x "$ANVIO_ENV_PREFIX/bin/python" ]]; then
    ANVIO_PYTHON="$ANVIO_ENV_PREFIX/bin/python"
elif command -v python3 >/dev/null 2>&1; then
    ANVIO_PYTHON="$(command -v python3)"
else
    fail "anvi'o Python not found; install anvi'o in and activate the MAGverify Conda environment"
fi
ANVIO_VERSION="$("$ANVIO_PYTHON" - <<'PY_ANVIO_VERSION'
import anvio
print(anvio.__version__)
PY_ANVIO_VERSION
)"

FINAL_INPUTS=()

clip_one_assembler() {
    local assembler="$1" assembly="$2" hifi_all="$3"
    local source_fai asm_dir asm_prefix asm_raw asm_zero asm_raw_meta asm_default
    local asm_stringent asm_final asm_final_meta raw_base_fp source_fai_sha raw_fp
    local asm_final_fp output

    echo
    echo "[03_clipping][$assembler] Assembly: $assembly"
    [[ -s "$assembly" ]] || fail "assembly missing: $assembly"
    samtools faidx "$assembly"
    source_fai="${assembly}.fai"

    [[ -s "$hifi_all" ]] || fail "HiFi_all_BAM missing: $hifi_all"
    [[ -s "$(qc_bam_cache "$hifi_all")" ]] \
        || fail "HiFi BAM cache missing: $(qc_bam_cache "$hifi_all"); rerun map"
    samtools quickcheck -v "$hifi_all" || fail "samtools quickcheck failed: $hifi_all"
    qc_bam_matches_fai "$hifi_all" "$source_fai" \
        || fail "HiFi BAM/assembly mismatch: $assembler"
    ensure_bam_index "$hifi_all"

    asm_dir="$CLIP_ASSEMBLER_DIR/$assembler"
    mkdir -p "$asm_dir"
    asm_prefix="$asm_dir/anvio_find_misassemblies"
    asm_raw="${asm_prefix}-clipping.txt"
    asm_zero="${asm_prefix}-zero_cov.txt"
    asm_raw_meta="${asm_raw}.inputs.sha256"
    asm_default="$asm_dir/clipping.default.tsv"
    asm_stringent="$asm_dir/clipping.stringent.tsv"
    asm_final="$asm_dir/clipping.site.tsv"
    asm_final_meta="${asm_final}.inputs.sha256"

    raw_base_fp="$(qc_compute_anvio_clipping_raw_fingerprint \
        "$(qc_bam_cache "$hifi_all")" "$ANVIO_DEFAULT_CLIPPING_RATIO" \
        "$MIN_DIST_TO_END" "$ANVIO_CLIPPING_ALGORITHM_ID" \
        "$ANVIO_FIND_MISASSEMBLIES" "$ANVIO_VERSION")"
    source_fai_sha="$(qc_file_sha256 "$source_fai")"
    raw_fp="$(qc_hash_text \
        "$raw_base_fp" "assembler=$assembler" \
        "source_fai_sha256=$source_fai_sha")"

    if [[ "$CLIPPING_POSTPROCESS_ONLY" == "Y" ]]; then
        [[ -s "$asm_raw" && -s "$asm_zero" ]] \
            || fail "postprocess-only mode requires existing raw outputs for $assembler"
        qc_write_result_cache_atomic "$asm_raw_meta" "$raw_fp" "$asm_raw"
        echo "[03_clipping][$assembler] Reusing raw anvi'o output by explicit postprocess-only request."
    elif [[ "$FORCE_RERUN_CLIPPING" == "Y" ]] \
      || ! qc_result_cache_matches "$asm_raw_meta" "$raw_fp" "$asm_raw" \
      || [[ ! -s "$asm_zero" ]]; then
        echo "[03_clipping][$assembler] Running anvi-script-find-misassemblies with --min-dist-to-end 0."
        rm -f "$asm_raw" "$asm_zero" "$asm_raw_meta"
        "$ANVIO_FIND_MISASSEMBLIES" \
            -b "$hifi_all" \
            -o "$asm_prefix" \
            --clipping-ratio "$ANVIO_DEFAULT_CLIPPING_RATIO" \
            --min-dist-to-end "$MIN_DIST_TO_END" \
            -T "$ANVIO_THREADS" \
            --just-do-it
        [[ -s "$asm_raw" ]] || fail "anvi'o clipping output missing: $asm_raw"
        [[ -s "$asm_zero" ]] || fail "anvi'o zero-coverage output missing: $asm_zero"
        qc_write_result_cache_atomic "$asm_raw_meta" "$raw_fp" "$asm_raw"
    else
        echo "[03_clipping][$assembler] Reusing fingerprinted raw anvi'o output."
    fi

    asm_final_fp="$(qc_compute_anvio_clipping_result_fingerprint \
        "$asm_raw_meta" "$asm_raw" \
        "$STRINGENT_CLIP_MIN_DEPTH" "$STRINGENT_CLIP_RATIO" \
        "$CLIPPING_NORMALIZATION_ALGORITHM_ID" "$CLIPPING_VALIDATION_ALGORITHM_ID" \
        "$CLIP_VALIDATE_SCRIPT" "$CLIP_VALIDATE_SCRIPT")"
    asm_final_fp="$(qc_hash_text "$asm_final_fp" "assembler=$assembler" "reference_boundary_policy=unresolved")"

    if [[ "$FORCE_RERUN_CLIPPING" != "Y" && "$CLIPPING_POSTPROCESS_ONLY" != "Y" \
          && -s "$asm_default" && -s "$asm_stringent" ]] \
       && qc_result_cache_matches "$asm_final_meta" "$asm_final_fp" "$asm_final"; then
        echo "[03_clipping][$assembler] Reusing read-level clipping result."
    else
        echo "[03_clipping][$assembler] Normalizing and validating exact HiFi clipping positions."
        rm -f "$asm_default" "$asm_stringent" "$asm_final" "$asm_final_meta"
        python3 "$CLIP_VALIDATE_SCRIPT" normalize \
            --anvio-clipping "$asm_raw" \
            --fai "$source_fai" \
            --default-output "$asm_default" \
            --stringent-output "$asm_stringent" \
            --stringent-min-depth "$STRINGENT_CLIP_MIN_DEPTH" \
            --stringent-ratio "$STRINGENT_CLIP_RATIO"
        python3 "$CLIP_VALIDATE_SCRIPT" validate \
            --candidates "$asm_default" \
            --bam "$hifi_all" \
            --output "$asm_final" \
            --threads "$CLIPPING_THREADS" \
            --high-confidence-min-coverage "$STRINGENT_CLIP_MIN_DEPTH" \
            --high-confidence-ratio "$STRINGENT_CLIP_RATIO"
        qc_write_result_cache_atomic "$asm_final_meta" "$asm_final_fp" "$asm_final"
    fi

    for output in "$asm_default" "$asm_stringent" "$asm_final"; do
        [[ -s "$output" ]] || fail "assembler clipping output missing: $output"
    done

    : > "$asm_dir/.clipping_done"
}

CLIP_ROWS=()
while IFS=$'\t' read -r assembler assembly hifi_all _rest; do
    [[ "$assembler" == "Assembler" || -z "$assembler" ]] && continue
    CLIP_ROWS+=("$assembler"$'\t'"$assembly"$'\t'"$hifi_all")
done < "$ASSEMBLER_MAPPING_MANIFEST"
CLIP_JOBS="$(qc_parallel_jobs "${#CLIP_ROWS[@]}")"
ANVIO_THREADS="$(qc_split_threads "$ANVIO_THREADS" "$CLIP_JOBS")"
CLIPPING_THREADS="$(qc_split_threads "$CLIPPING_THREADS" "$CLIP_JOBS")"
echo "[03_clipping] Parallel assembler jobs: $CLIP_JOBS (anvi'o -T $ANVIO_THREADS, validation threads $CLIPPING_THREADS per job)"
for row in "${CLIP_ROWS[@]}"; do
    IFS=$'\t' read -r assembler assembly hifi_all <<< "$row"
    mkdir -p "$CLIP_ASSEMBLER_DIR/$assembler"
    rm -f "$CLIP_ASSEMBLER_DIR/$assembler/.clipping_done"
    ( clip_one_assembler "$assembler" "$assembly" "$hifi_all" ) &
    qc_wait_for_slot "$CLIP_JOBS"
done
qc_wait_all
for row in "${CLIP_ROWS[@]}"; do
    IFS=$'\t' read -r assembler assembly hifi_all <<< "$row"
    asm_dir="$CLIP_ASSEMBLER_DIR/$assembler"
    [[ -e "$asm_dir/.clipping_done" ]] || fail "clipping failed for $assembler"
    FINAL_INPUTS+=("$assembler=$asm_dir/clipping.site.tsv")
done

(( ${#FINAL_INPUTS[@]} > 0 )) || fail "no assembler rows were processed"

CMD=(
    python3 "$CLIP_SUMMARIZE_SCRIPT" sites
    --mapping-manifest "$ASSEMBLER_MAPPING_MANIFEST"
    --assemblies-fai "${ASSEMBLY_FA}.fai"
    --terminal-window "$TERMINAL_CLIPPING_WINDOW_BP"
)
for input in "${FINAL_INPUTS[@]}"; do CMD+=(--input "$input"); done
CMD+=(--all-output "$CLIP_ALL_FINAL" --internal-output "$CLIP_FINAL")
"${CMD[@]}"

python3 "$CLIP_SUMMARIZE_SCRIPT" assembler \
    --mapping-manifest "$ASSEMBLER_MAPPING_MANIFEST" \
    --terminal-window "$TERMINAL_CLIPPING_WINDOW_BP" \
    --sites "$CLIP_ALL_FINAL" \
    --output "$CLIP_ASSEMBLER_SUMMARY"

GLOBAL_FINAL_FP="$(qc_compute_source_specific_clipping_fingerprint \
    "$ASSEMBLER_MAPPING_MANIFEST" "$ANVIO_DEFAULT_CLIPPING_RATIO" \
    "$MIN_DIST_TO_END" "$TERMINAL_CLIPPING_WINDOW_BP" \
    "$STRINGENT_CLIP_MIN_DEPTH" "$STRINGENT_CLIP_RATIO" \
    "$ANVIO_CLIPPING_ALGORITHM_ID" "$CLIPPING_NORMALIZATION_ALGORITHM_ID" \
    "$CLIPPING_VALIDATION_ALGORITHM_ID" "$CLIP_VALIDATE_SCRIPT" \
    "$CLIP_VALIDATE_SCRIPT" "$CLIP_SUMMARIZE_SCRIPT" \
    "$ANVIO_FIND_MISASSEMBLIES" "$ANVIO_VERSION")"
qc_write_result_cache_atomic "$CLIP_FINAL_CACHE_META" "$GLOBAL_FINAL_FP" "$CLIP_FINAL"

write_mag_outputs

echo
echo "[03_clipping] Per-assembler clipping analysis complete: $(date)"
echo "SCRIPT_ID              : $SCRIPT_ID"
echo "Candidate ratio             : anvi'o documented default >$ANVIO_DEFAULT_CLIPPING_RATIO"
echo "Terminal exclusion          : disabled (--min-dist-to-end $MIN_DIST_TO_END)"
echo "Terminal classification     : first/last $TERMINAL_CLIPPING_WINDOW_BP bp"
echo "HiFi mapping context        : all contigs of each assembler, all mapped reads"
echo "Permissive clipping audit    : Clipping_ge1=Y requires >=1 distinct clipped HiFi read + zero spanning reads"
echo "Final clipping criterion     : Clipping=Y requires >=2 distinct clipped HiFi reads + zero spanning reads"
echo "I/D/N policy                 : audit only; spanning requires mapped reference sequence on both sides regardless of I/D/N pattern"
echo "Stringent audit tier         : Clipping=Y + coverage >=$STRINGENT_CLIP_MIN_DEPTH + clipping ratio =$STRINGENT_CLIP_RATIO"
echo "Exact reference boundaries  : retained for terminal statistics; BOUNDARY_UNRESOLVED; never Clipping=Y"

python3 - "$CLIP_ALL_FINAL" "$CLIP_ASSEMBLER_SUMMARY" <<'PY_SUMMARY'
import csv
import sys
final_path, assembler_path = sys.argv[1:]
with open(final_path, encoding='utf-8-sig', newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
print(f"All-region clipping candidates      : {len(rows)}")
print(f"Terminal candidates                 : {sum(r.get('region_type') == 'TERMINAL' for r in rows)}")
print(f"Internal candidates                 : {sum(r.get('region_type') == 'INTERNAL' for r in rows)}")
print(f"Clipping_ge1=Y audit                : {sum(r.get('Clipping_ge1') == 'Y' for r in rows)}")
print(f"Clipping=Y (>=2 reads)              : {sum(r.get('Clipping') == 'Y' for r in rows)}")
print(f"Stringent_clipping=Y                : {sum(r.get('Stringent_clipping') == 'Y' for r in rows)}")
print(f"Terminal Clipping=Y                 : {sum(r.get('Clipping') == 'Y' and r.get('region_type') == 'TERMINAL' for r in rows)}")
print(f"Internal Clipping=Y                 : {sum(r.get('Clipping') == 'Y' and r.get('region_type') == 'INTERNAL' for r in rows)}")
print(f"Reference-boundary raw candidates   : {sum(r.get('reference_boundary') == 'Y' for r in rows)}")
with open(assembler_path, encoding='utf-8-sig', newline='') as handle:
    summary = list(csv.DictReader(handle, delimiter='\t'))
all_row = next((r for r in summary if r.get('Assembler') == 'ALL'), None)
if all_row:
    print("Terminal/internal Clipping=Y density ratio: " + all_row['terminal_to_internal_Clipping_Y_sites_density_ratio'])
PY_SUMMARY

for legacy in "${LEGACY_CLIPPING_FILES[@]}"; do
    rm -f "${CLIPPING_DIR}/$legacy"
done

print_outputs
echo "[03_clipping] Finished: $(date)"
