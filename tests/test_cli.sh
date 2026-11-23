#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../modules" && pwd)"
TEST_ROOT="$(mktemp -d)"
trap 'rm -rf "$TEST_ROOT"' EXIT

check_config() {
    local stage="$1" expected_input_variable="$2" input="$3" out="$4" expected_output_variable="$5" expected_output="$6"
    QC_STAGE="$stage" bash -c '
      set -euo pipefail
      config="$1"; input="$2"; out="$3"; input_variable="$4"; output_variable="$5"; expected_output="$6"
      set -- "$input" "$out"
      source "$config/lib_config.sh"
      [[ "${!input_variable}" == "$input" ]]
      [[ "${!output_variable}" == "$expected_output" ]]
    ' _ "$HERE" "$input" "$out" "$expected_input_variable" "$expected_output_variable" "$expected_output"
}

stages=(01_mapping 03_clipping 07_chimera 02_circularity 04_false_closure 09_assembly_variants 05_phantom 06_repeats 08_false_duplication 10_summary)
inputs=(READS ASSEMBLY_FA GTDB_GENE_TAX_TSV MAG_DIR ASSEMBLY_FA ASSEMBLY_FA ASSEMBLY_FA ASSEMBLY_FA ASSEMBLY_FA ASSEMBLY_FA)
outputs=(MAPPING_DIR CLIPPING_DIR CHIMERA_ROOT_DIR CIRCULARITY_DIR FALSE_CLOSURE_DIR VARIANT_DIR PHANTOM_DIR REPEAT_DIR FALSE_DUP_DIR SUMMARY_DIR)
for i in "${!stages[@]}"; do
    stage="${stages[$i]}"
    in="$TEST_ROOT/inputs/${stage}.fa"
    out="$TEST_ROOT/analysis/${stage}_custom"
    check_config "$stage" "${inputs[$i]}" "$in" "$out" "${outputs[$i]}" "$out"
done

QC_STAGE=03_clipping bash -c '
  config="$1"; set -- /inputs/combined.fa /outputs/clipping
  source "$config/lib_config.sh"
  [[ "$CLIP_FINAL" == /outputs/clipping/.cache/clipping.site.internal.tsv ]]
  [[ "$CLIP_ALL_FINAL" == /outputs/clipping/.cache/clipping.site.tsv ]]
  [[ "$OUT" == /outputs ]]
' _ "$HERE"
QC_STAGE=01_mapping bash -c '
  config="$1"; set -- /inputs/reads.fastq.gz /outputs/mapping
  source "$config/lib_config.sh"
  [[ "$ASSEMBLER_MAPPING_MANIFEST" == /outputs/mapping/.cache/mapping_manifest.tsv ]]
' _ "$HERE"
QC_STAGE=07_chimera bash -c '
  config="$1"; set -- /inputs/taxonomy.tsv /outputs/chimera
  source "$config/lib_config.sh"
  [[ "$GTDB_CHIMERA_DIR" == /outputs/chimera/.cache ]]
' _ "$HERE"
QC_STAGE=02_circularity bash -c '
  config="$1"; set -- /inputs/mags /outputs/circularity
  source "$config/lib_config.sh"
  [[ "$CIRCULARITY_ROTATED_RESULT_DIR" == /outputs/circularity ]]
' _ "$HERE"
QC_STAGE=10_summary bash -c '
  config="$1"; set -- /inputs/combined.fa /outputs/summary
  source "$config/lib_config.sh"
  [[ "$SUMMARY_MAG" == /outputs/summary/summary.MAG.tsv ]]
' _ "$HERE"

QC_STAGE=03_clipping CLIPPING_INPUT_ASSEMBLY_FA=/other/assembly.fa bash -c '
  config="$1"; set -- /cli/assembly.fa /cli/clipping
  source "$config/lib_config.sh"
  [[ "$ASSEMBLY_FA" == /cli/assembly.fa ]]
  [[ "$CLIPPING_DIR" == /cli/clipping ]]
' _ "$HERE"

QC_STAGE=02_circularity OUT="$TEST_ROOT/flag_out" bash -c '
  config="$1"; set -- --mag-dir /flag/mags --outdir /flag/out
  source "$config/lib_config.sh"
  [[ "$MAG_DIR" == "${MAG_DIR:-}" ]]
  [[ "$CIRCULARITY_DIR" == "$OUT/02_circularity" ]]
' _ "$HERE"

if QC_STAGE=06_repeats bash -c 'config="$1"; unset OUT; source "$config/lib_config.sh"' _ "$HERE" >/dev/null 2>&1; then
    echo 'ERROR: lib_config.sh accepted a missing output directory' >&2
    exit 1
fi
[[ ! -e "$HERE/results" ]] || { echo "ERROR: $HERE/results was created" >&2; exit 1; }

if QC_STAGE=03_clipping bash -c 'config="$1"; set -- /input.fa; source "$config/lib_config.sh"' _ "$HERE" >/dev/null 2>&1; then
    echo 'ERROR: single-argument clipping invocation unexpectedly succeeded' >&2
    exit 1
fi

if bash "$HERE/07_chimera.sh" __taxonomy_internal /inputs/assembly.fa /outputs/taxonomy \
    >"$TEST_ROOT/taxonomy.log" 2>&1; then
    echo 'ERROR: taxonomy unexpectedly accepted missing required settings' >&2
    exit 1
fi
grep -q 'MMSEQS_TAX_DB and MAG_CONTIG_MAP must be set' "$TEST_ROOT/taxonomy.log"

echo 'Positional input/output configuration tests passed.'
