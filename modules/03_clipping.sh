#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_ID="03_clipping"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found on PATH" >&2; exit 1; }
python3 "$SCRIPT_DIR/../scripts/clipping_check_clipping_inputs.py" \
    --assembly "$ASSEMBLY_FA" \
    --mag-dir "$MAG_DIR" \
    --manifest "$ASSEMBLER_MAPPING_MANIFEST"

CLIP_ALL_DEFAULT="${CLIPPING_DIR}/clipping.all_regions.default.tsv"
CLIP_INTERNAL_DEFAULT="${CLIPPING_DIR}/clipping.internal.default.tsv"
CLIP_TERMINAL_DEFAULT="${CLIPPING_DIR}/clipping.terminal.default.tsv"
CLIP_ALL_STRINGENT="${CLIPPING_DIR}/clipping.all_regions.stringent.tsv"
CLIP_INTERNAL_STRINGENT="${CLIPPING_DIR}/clipping.internal.stringent.tsv"
CLIP_TERMINAL_STRINGENT="${CLIPPING_DIR}/clipping.terminal.stringent.tsv"
CLIP_ALL_FINAL="${CLIPPING_DIR}/clipping.all_regions.filtered.tsv"
CLIP_INTERNAL_FINAL="${CLIPPING_DIR}/clipping.internal.filtered.tsv"
CLIP_TERMINAL_FINAL="${CLIPPING_DIR}/clipping.terminal.filtered.tsv"
CLIP_ALL_POSITIVE="${CLIPPING_DIR}/clipping.all_regions.Y.tsv"
CLIP_INTERNAL_POSITIVE="${CLIPPING_DIR}/clipping.internal.Y.tsv"
CLIP_TERMINAL_POSITIVE="${CLIPPING_DIR}/clipping.terminal.Y.tsv"
CLIP_REGION_SUMMARY="${CLIPPING_DIR}/clipping.region_density.tsv"
CLIP_REGION_RATIO="${CLIPPING_DIR}/clipping.region_density_ratio.tsv"
CLIP_SOURCE_MANIFEST="${CLIPPING_DIR}/clipping.source_mapping_manifest.tsv"

mkdir -p "$OUT/00_log" "${CLIPPING_DIR}" "${CLIPPING_DIR}/by_assembler"

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

combine_tsvs() {
    local output="$1"
    shift
    python3 - "$output" "$@" <<'PY_COMBINE_TSV'
import csv
import os
import sys
from pathlib import Path

out, *paths = sys.argv[1:]
if not paths:
    raise SystemExit("ERROR: no TSV inputs supplied")
header = None
rows = []
for path in paths:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"ERROR: TSV input missing or empty: {path}")
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            fields = next(reader)
        except StopIteration:
            raise SystemExit(f"ERROR: TSV has no header: {path}")
        if header is None:
            header = fields
        elif fields != header:
            raise SystemExit(f"ERROR: TSV schemas differ: {path}")
        rows.extend(reader)
Path(out).parent.mkdir(parents=True, exist_ok=True)
tmp = out + f".tmp.{os.getpid()}"
with open(tmp, "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
os.replace(tmp, out)
PY_COMBINE_TSV
}

write_mag_outputs() {
    local clipping_input="$CLIP_ALL_FINAL"
    python3 "$CLIP_MAG_SUMMARIZER_SCRIPT" \
        --clipping "$clipping_input" \
        --mag-contig-map "$MAG_CONTIG_MAP" \
        --fai "${ASSEMBLY_FA}.fai" \
        --positive-site-output "$CLIP_POSITIVE" \
        --mag-site-output "$CLIP_MAG_SITE" \
        --mag-positive-site-output "$CLIP_MAG_POSITIVE_SITE" \
        --summary-output "$CLIP_MAG_SUMMARY" \
        --positive-summary-output "$CLIP_MAG_POSITIVE" \
        --contig-summary-output "$CLIP_CONTIG_SUMMARY"
}

print_outputs() {
    echo "Outputs:"
    for path in \
        "$CLIP_SOURCE_MANIFEST" \
        "$ANVIO_CLIP_RAW" "$ANVIO_ZERO_RAW" \
        "$CLIP_ALL_DEFAULT" "$CLIP_INTERNAL_DEFAULT" "$CLIP_TERMINAL_DEFAULT" \
        "$CLIP_ALL_STRINGENT" "$CLIP_INTERNAL_STRINGENT" "$CLIP_TERMINAL_STRINGENT" \
        "$CLIP_ALL_FINAL" "$CLIP_INTERNAL_FINAL" "$CLIP_TERMINAL_FINAL" \
        "$CLIP_ALL_POSITIVE" "$CLIP_INTERNAL_POSITIVE" "$CLIP_TERMINAL_POSITIVE" \
        "$CLIP_REGION_SUMMARY" "$CLIP_REGION_RATIO" \
        "$CLIP_DEFAULT" "$CLIP_STRINGENT" "$CLIP_FINAL" "$CLIP_POSITIVE" \
        "$CLIP_MAG_SITE" "$CLIP_MAG_POSITIVE_SITE" \
        "$CLIP_MAG_SUMMARY" "$CLIP_MAG_POSITIVE" "$CLIP_CONTIG_SUMMARY" "$CLIP_CRITERIA"; do
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

samtools faidx "$ASSEMBLY_FA"
if [[ ! -s "$MAG_CONTIG_MAP" ]]; then
    echo "[03_clipping] MAG_CONTIG_MAP is absent; creating shared prep map."
fi
python3 "$MAG_MAP_BUILDER" \
    --mag-dir "$MAG_DIR" \
    --extensions "$MAG_EXTENSIONS" \
    --assembly-fai "${ASSEMBLY_FA}.fai" \
    --output "$MAG_CONTIG_MAP"

if [[ ! -s "$ASSEMBLER_MAPPING_MANIFEST" ]]; then
    fail "assembler mapping manifest is missing: $ASSEMBLER_MAPPING_MANIFEST. Run 01_mapping.sh before 03_clipping.sh."
fi

for required in \
    "$MAG_CONTIG_MAP" "$ASSEMBLER_MAPPING_MANIFEST" \
    "$CLIP_NORMALIZER_SCRIPT" "$CLIP_SPANNING_VALIDATOR" \
    "$CLIP_MAG_SUMMARIZER_SCRIPT" "$CLIP_REGION_CLASSIFIER_SCRIPT" \
    "$CLIP_REGION_SUMMARIZER_SCRIPT"; do
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

copy_atomic "$ASSEMBLER_MAPPING_MANIFEST" "$CLIP_SOURCE_MANIFEST"

while IFS=$'\t' read -r assembler assembly hifi_all; do
    [[ -s "$assembly" ]] || fail "source assembly missing: $assembly"
    samtools faidx "$assembly"
    [[ -s "$hifi_all" && -s "$(qc_bam_cache "$hifi_all")" ]] \
        || fail "source-assembler HiFi mapping missing for $assembler; rerun 01_mapping.sh"
    samtools quickcheck -v "$hifi_all" \
        || fail "samtools quickcheck failed: $hifi_all"
    qc_bam_matches_fai "$hifi_all" "${assembly}.fai" \
        || fail "HiFi BAM/source assembly mismatch: $assembler"
    ensure_bam_index "$hifi_all"
done < <(qc_manifest_columns "$ASSEMBLER_MAPPING_MANIFEST" Assembler Assembly_FASTA HiFi_all_BAM)

python3 - "$ASSEMBLER_MAPPING_MANIFEST" "${ASSEMBLY_FA}.fai" <<'PY_SOURCE_UNION'
import csv
import sys
manifest, combined_fai = sys.argv[1:]

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

combined = fai(combined_fai)
union = {}
with open(manifest, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        assembler = row['Assembler'].strip()
        source = fai(row['Assembly_FASTA'].strip() + '.fai')
        for contig, length in source.items():
            if contig in union:
                raise SystemExit(f'ERROR: contig ID occurs in multiple source assemblers: {contig}')
            union[contig] = (assembler, length)
if set(union) != set(combined):
    missing = sorted(set(combined) - set(union))[:20]
    extra = sorted(set(union) - set(combined))[:20]
    raise SystemExit(
        'ERROR: source-assembler contig union differs from concatenated assembly; '
        f'missing_from_sources={missing}, extra_in_sources={extra}'
    )
for contig, length in combined.items():
    if union[contig][1] != length:
        raise SystemExit(f'ERROR: source/concatenated contig length mismatch: {contig}')
print(f'Validated source-assembler contig union: {len(combined)} contigs')
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

DEFAULT_INPUTS=()
STRINGENT_INPUTS=()
FINAL_INPUTS=()
RAW_CLIP_FILES=()
RAW_ZERO_FILES=()

clip_one_assembler() {
    local assembler="$1" assembly="$2" hifi_all="$3"
    local source_fai asm_dir asm_prefix asm_raw asm_zero asm_raw_meta asm_default
    local asm_stringent asm_final asm_final_meta raw_base_fp source_fai_sha raw_fp
    local asm_final_fp output

    echo
    echo "[03_clipping][$assembler] Source assembly: $assembly"
    [[ -s "$assembly" ]] || fail "source assembly missing: $assembly"
    samtools faidx "$assembly"
    source_fai="${assembly}.fai"

    [[ -s "$hifi_all" ]] || fail "source-assembler HiFi_all_BAM missing: $hifi_all"
    [[ -s "$(qc_bam_cache "$hifi_all")" ]] \
        || fail "source-assembler HiFi BAM cache missing: $(qc_bam_cache "$hifi_all"); rerun map"
    samtools quickcheck -v "$hifi_all" || fail "samtools quickcheck failed: $hifi_all"
    qc_bam_matches_fai "$hifi_all" "$source_fai" \
        || fail "HiFi BAM/source assembly mismatch: $assembler"
    ensure_bam_index "$hifi_all"

    asm_dir="${CLIPPING_DIR}/by_assembler/$assembler"
    mkdir -p "$asm_dir"
    asm_prefix="$asm_dir/anvio_find_misassemblies"
    asm_raw="${asm_prefix}-clipping.txt"
    asm_zero="${asm_prefix}-zero_cov.txt"
    asm_raw_meta="${asm_raw}.inputs.sha256"
    asm_default="$asm_dir/clipping.default.tsv"
    asm_stringent="$asm_dir/clipping.stringent.tsv"
    asm_final="$asm_dir/clipping.filtered.tsv"
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
            || fail "postprocess-only mode requires existing source raw outputs for $assembler"
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
        "$CLIP_NORMALIZER_SCRIPT" "$CLIP_SPANNING_VALIDATOR")"
    asm_final_fp="$(qc_hash_text "$asm_final_fp" "assembler=$assembler" "reference_boundary_policy=unresolved")"

    if [[ "$FORCE_RERUN_CLIPPING" != "Y" && "$CLIPPING_POSTPROCESS_ONLY" != "Y" \
          && -s "$asm_default" && -s "$asm_stringent" ]] \
       && qc_result_cache_matches "$asm_final_meta" "$asm_final_fp" "$asm_final"; then
        echo "[03_clipping][$assembler] Reusing read-level clipping result."
    else
        echo "[03_clipping][$assembler] Normalizing and validating exact HiFi clipping positions."
        rm -f "$asm_default" "$asm_stringent" "$asm_final" "$asm_final_meta"
        python3 "$CLIP_NORMALIZER_SCRIPT" \
            --anvio-clipping "$asm_raw" \
            --fai "$source_fai" \
            --default-output "$asm_default" \
            --stringent-output "$asm_stringent" \
            --stringent-min-depth "$STRINGENT_CLIP_MIN_DEPTH" \
            --stringent-ratio "$STRINGENT_CLIP_RATIO"
        python3 "$CLIP_SPANNING_VALIDATOR" \
            --candidates "$asm_default" \
            --bam "$hifi_all" \
            --output "$asm_final" \
            --threads "$CLIPPING_THREADS" \
            --high-confidence-min-coverage "$STRINGENT_CLIP_MIN_DEPTH" \
            --high-confidence-ratio "$STRINGENT_CLIP_RATIO"
        qc_write_result_cache_atomic "$asm_final_meta" "$asm_final_fp" "$asm_final"
    fi

    for output in "$asm_default" "$asm_stringent" "$asm_final"; do
        [[ -s "$output" ]] || fail "source clipping output missing: $output"
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
    mkdir -p "${CLIPPING_DIR}/by_assembler/$assembler"
    rm -f "${CLIPPING_DIR}/by_assembler/$assembler/.clipping_done"
    ( clip_one_assembler "$assembler" "$assembly" "$hifi_all" ) &
    qc_wait_for_slot "$CLIP_JOBS"
done
wait
for row in "${CLIP_ROWS[@]}"; do
    IFS=$'\t' read -r assembler assembly hifi_all <<< "$row"
    asm_dir="${CLIPPING_DIR}/by_assembler/$assembler"
    [[ -e "$asm_dir/.clipping_done" ]] || fail "clipping failed for $assembler"
    DEFAULT_INPUTS+=("$assembler=$asm_dir/clipping.default.tsv")
    STRINGENT_INPUTS+=("$assembler=$asm_dir/clipping.stringent.tsv")
    FINAL_INPUTS+=("$assembler=$asm_dir/clipping.filtered.tsv")
    RAW_CLIP_FILES+=("$asm_dir/anvio_find_misassemblies-clipping.txt")
    RAW_ZERO_FILES+=("$asm_dir/anvio_find_misassemblies-zero_cov.txt")
done

(( ${#FINAL_INPUTS[@]} > 0 )) || fail "no source assembler rows were processed"

combine_tsvs "$ANVIO_CLIP_RAW" "${RAW_CLIP_FILES[@]}"
combine_tsvs "$ANVIO_ZERO_RAW" "${RAW_ZERO_FILES[@]}"

CLASSIFY_CMD=(
    python3 "$CLIP_REGION_CLASSIFIER_SCRIPT"
    --mapping-manifest "$ASSEMBLER_MAPPING_MANIFEST"
    --combined-fai "${ASSEMBLY_FA}.fai"
    --terminal-window "$TERMINAL_CLIPPING_WINDOW_BP"
)
CMD=("${CLASSIFY_CMD[@]}")
for input in "${DEFAULT_INPUTS[@]}"; do CMD+=(--input "$input"); done
CMD+=(
    --all-output "$CLIP_ALL_DEFAULT"
    --internal-output "$CLIP_INTERNAL_DEFAULT"
    --terminal-output "$CLIP_TERMINAL_DEFAULT"
)
"${CMD[@]}"

CMD=("${CLASSIFY_CMD[@]}")
for input in "${STRINGENT_INPUTS[@]}"; do CMD+=(--input "$input"); done
CMD+=(
    --all-output "$CLIP_ALL_STRINGENT"
    --internal-output "$CLIP_INTERNAL_STRINGENT"
    --terminal-output "$CLIP_TERMINAL_STRINGENT"
)
"${CMD[@]}"

CMD=("${CLASSIFY_CMD[@]}")
for input in "${FINAL_INPUTS[@]}"; do CMD+=(--input "$input"); done
CMD+=(
    --all-output "$CLIP_ALL_FINAL"
    --internal-output "$CLIP_INTERNAL_FINAL"
    --terminal-output "$CLIP_TERMINAL_FINAL"
    --positive-output "$CLIP_ALL_POSITIVE"
    --internal-positive-output "$CLIP_INTERNAL_POSITIVE"
    --terminal-positive-output "$CLIP_TERMINAL_POSITIVE"
)
"${CMD[@]}"

copy_atomic "$CLIP_INTERNAL_DEFAULT" "$CLIP_DEFAULT"
copy_atomic "$CLIP_INTERNAL_STRINGENT" "$CLIP_STRINGENT"
copy_atomic "$CLIP_INTERNAL_FINAL" "$CLIP_FINAL"

python3 "$CLIP_REGION_SUMMARIZER_SCRIPT" \
    --mapping-manifest "$ASSEMBLER_MAPPING_MANIFEST" \
    --terminal-window "$TERMINAL_CLIPPING_WINDOW_BP" \
    --default "$CLIP_ALL_DEFAULT" \
    --stringent "$CLIP_ALL_STRINGENT" \
    --final "$CLIP_ALL_FINAL" \
    --summary-output "$CLIP_REGION_SUMMARY" \
    --ratio-output "$CLIP_REGION_RATIO"

cat > "$CLIP_CRITERIA" <<EOF_CRITERIA
field	definition	alignment_population	counting_unit
source_reference	Full original co-assembly from the assembler that produced the evaluated contig	assembler-specific	contig
clipping_candidate	Position reported by anvi-script-find-misassemblies with its documented default clipping-ratio threshold >0.8; --min-dist-to-end is intentionally set to 0 and no minimum coverage is imposed	HiFi ALL mapped	exact base
region_type	TERMINAL for the first/last ${TERMINAL_CLIPPING_WINDOW_BP} bp of a contig; INTERNAL otherwise; classification occurs after identical clipping detection/validation	NA	exact base
reference_boundary	Y only at exact first/last reference base; retained for raw terminal statistics but two-sided spanning is impossible and Clipping cannot become Y	NA	exact base
exact_position_validation	PASS when a non-boundary candidate has >=2 distinct exact-position clipped HiFi QNAMEs and zero spanning QNAMEs; INSUFFICIENT_CLIPPED_QNAMES when only the >=1-read audit criterion is met; BOUNDARY_UNRESOLVED at exact reference edges; FAIL otherwise	HiFi ALL mapped	unique QNAME
Clipping_ge1	Y for the permissive audit call when a non-boundary site has >=1 distinct clipped HiFi QNAME and zero spanning QNAMEs	HiFi ALL mapped	exact base
exact_position_spanning	Main counter-evidence: the same HiFi alignment contains actual M/= /X aligned reference bases on both sides of the candidate; I/D/N patterns do not alter this call and are retained only as audit columns	HiFi ALL mapped	unique QNAME
exact_spanning_min_flank_bp_min/median/max	Distribution of the shorter M/= /X-aligned flank among spanning HiFi QNAMEs. For each QNAME, the alignment with the largest shorter flank is retained; audit only and never alters Clipping	HiFi ALL mapped	bp per unique QNAME
n_unique_spanning_with_I_qnames	Spanning HiFi QNAMEs with a candidate-adjacent insertion; audit only and never alters Clipping	HiFi ALL mapped	unique QNAME
spanning_I_length_bp_min/median/max	Candidate-adjacent insertion length among spanning HiFi QNAMEs. For each QNAME, the largest candidate-local insertion observed across its spanning alignments is retained; audit only	HiFi ALL mapped	bp per unique QNAME
n_unique_spanning_with_D_qnames	Spanning HiFi QNAMEs with a deletion consuming the candidate base; audit only and never alters Clipping	HiFi ALL mapped	unique QNAME
spanning_D_length_bp_min/median/max	Length of the deletion consuming the candidate base among spanning HiFi QNAMEs. For each QNAME, the largest observed candidate-local deletion is retained; audit only	HiFi ALL mapped	bp per unique QNAME
n_unique_spanning_with_N_qnames	Spanning HiFi QNAMEs with a reference skip consuming the candidate base; audit only and never alters Clipping	HiFi ALL mapped	unique QNAME
spanning_N_length_bp_min/median/max	Length of the reference skip consuming the candidate base among spanning HiFi QNAMEs. For each QNAME, the largest observed candidate-local skip is retained; audit only	HiFi ALL mapped	bp per unique QNAME
n_unique_spanning_with_any_IDN_qnames	Union of spanning HiFi QNAMEs carrying candidate-local I, D or N audit events	HiFi ALL mapped	unique QNAME
stringent_candidate	anvi'o coverage >=${STRINGENT_CLIP_MIN_DEPTH} and clipping ratio =${STRINGENT_CLIP_RATIO}	HiFi ALL mapped	exact base
clipping_confidence	STRINGENT for validated sites meeting the stringent audit threshold; DEFAULT_ONLY for validated sites below that audit threshold; NOT_VALIDATED otherwise	HiFi ALL mapped	site
Clipping	Final call: Y for non-boundary sites with >=2 distinct clipped HiFi QNAMEs and zero spanning QNAMEs	HiFi ALL mapped	exact base
Stringent_clipping	Y only when Clipping=Y plus coverage >=${STRINGENT_CLIP_MIN_DEPTH} and clipping ratio =${STRINGENT_CLIP_RATIO}; retained as a separate audit call	HiFi ALL mapped	exact base
Internal_clipping_ge1	MAG-level audit Y when >=1 Clipping_ge1=Y site occurs in an INTERNAL region; N otherwise	HiFi ALL mapped	MAG
Terminal_clipping_ge1	MAG-level audit Y when >=1 Clipping_ge1=Y site occurs in a TERMINAL region; N otherwise	HiFi ALL mapped	MAG
MAG_Clipping_ge1	MAG-level audit call reported as Clipping_ge1=Y when Internal_clipping_ge1=Y or Terminal_clipping_ge1=Y; N otherwise	HiFi ALL mapped	MAG
Internal_clipping	MAG-level final Y when >=1 final Clipping=Y site occurs in an INTERNAL region; N otherwise	HiFi ALL mapped	MAG
Terminal_clipping	MAG-level final Y when >=1 final Clipping=Y site occurs in a TERMINAL region; N otherwise	HiFi ALL mapped	MAG
MAG_Clipping	MAG-level final call reported as Clipping=Y when Internal_clipping=Y or Terminal_clipping=Y; N otherwise	HiFi ALL mapped	MAG
terminal_internal_density	Site counts are normalized by the total analyzed bp in each region and reported per Mbp; terminal/internal ratios are enrichment ratios, not direct error-rate ratios	NA	region
EOF_CRITERIA

GLOBAL_FINAL_FP="$(qc_compute_source_specific_clipping_fingerprint \
    "$ASSEMBLER_MAPPING_MANIFEST" "$ANVIO_DEFAULT_CLIPPING_RATIO" \
    "$MIN_DIST_TO_END" "$TERMINAL_CLIPPING_WINDOW_BP" \
    "$STRINGENT_CLIP_MIN_DEPTH" "$STRINGENT_CLIP_RATIO" \
    "$ANVIO_CLIPPING_ALGORITHM_ID" "$CLIPPING_NORMALIZATION_ALGORITHM_ID" \
    "$CLIPPING_VALIDATION_ALGORITHM_ID" "$CLIP_NORMALIZER_SCRIPT" \
    "$CLIP_SPANNING_VALIDATOR" "$CLIP_REGION_CLASSIFIER_SCRIPT" \
    "$ANVIO_FIND_MISASSEMBLIES" "$ANVIO_VERSION")"
qc_write_result_cache_atomic "$CLIP_FINAL_CACHE_META" "$GLOBAL_FINAL_FP" "$CLIP_FINAL"

write_mag_outputs

echo
echo "[03_clipping] Source-specific clipping analysis complete: $(date)"
echo "SCRIPT_ID              : $SCRIPT_ID"
echo "Candidate ratio             : anvi'o documented default >$ANVIO_DEFAULT_CLIPPING_RATIO"
echo "Terminal exclusion          : disabled (--min-dist-to-end $MIN_DIST_TO_END)"
echo "Terminal classification     : first/last $TERMINAL_CLIPPING_WINDOW_BP bp"
echo "HiFi mapping context        : source-assembler full assembly, ALL mapped"
echo "Permissive clipping audit    : Clipping_ge1=Y requires >=1 distinct clipped HiFi read + zero spanning reads"
echo "Final clipping criterion     : Clipping=Y requires >=2 distinct clipped HiFi reads + zero spanning reads"
echo "I/D/N policy                 : audit only; spanning requires mapped reference sequence on both sides regardless of I/D/N pattern"
echo "Stringent audit tier         : Clipping=Y + coverage >=$STRINGENT_CLIP_MIN_DEPTH + clipping ratio =$STRINGENT_CLIP_RATIO"
echo "Exact reference boundaries  : retained for terminal statistics; BOUNDARY_UNRESOLVED; never Clipping=Y"

python3 - "$CLIP_ALL_FINAL" "$CLIP_REGION_RATIO" <<'PY_SUMMARY'
import csv
import sys
final_path, ratio_path = sys.argv[1:]
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
with open(ratio_path, encoding='utf-8-sig', newline='') as handle:
    ratios = list(csv.DictReader(handle, delimiter='\t'))
all_row = next((r for r in ratios if r.get('source_assembler') == 'ALL'), None)
if all_row:
    print("Terminal/internal Clipping=Y density ratio: " + all_row['terminal_to_internal_Clipping_Y_density_ratio'])
PY_SUMMARY

print_outputs
echo "[03_clipping] Finished: $(date)"
