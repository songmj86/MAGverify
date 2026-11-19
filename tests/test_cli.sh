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

# Derived filenames and subdirectories must follow positional output paths.
QC_STAGE=03_clipping bash -c '
  config="$1"; set -- /inputs/combined.fa /outputs/clipping
  source "$config/lib_config.sh"
  [[ "$CLIP_FINAL" == /outputs/clipping/clipping.filtered.tsv ]]
  [[ "$OUT" == /outputs ]]
' _ "$HERE"
QC_STAGE=01_mapping bash -c '
  config="$1"; set -- /inputs/reads.fastq.gz /outputs/mapping
  source "$config/lib_config.sh"
  [[ "$ASSEMBLER_MAPPING_MANIFEST" == /outputs/mapping/mapping_manifest.tsv ]]
' _ "$HERE"
QC_STAGE=07_chimera bash -c '
  config="$1"; set -- /inputs/taxonomy.tsv /outputs/chimera
  source "$config/lib_config.sh"
  [[ "$GTDB_CHIMERA_DIR" == /outputs/chimera/GTDB_R232 ]]
' _ "$HERE"
QC_STAGE=02_circularity bash -c '
  config="$1"; set -- /inputs/mags /outputs/circularity
  source "$config/lib_config.sh"
  [[ "$CIRCULARITY_ROTATED_RESULT_DIR" == /outputs/circularity ]]
' _ "$HERE"
QC_STAGE=10_summary bash -c '
  config="$1"; set -- /inputs/combined.fa /outputs/summary
  source "$config/lib_config.sh"
  [[ "$SUMMARY_MAG" == /outputs/summary/integrated_QC_summary.by_MAG.tsv ]]
' _ "$HERE"

# An explicit CLI input supersedes an environment-supplied per-module input.
QC_STAGE=03_clipping CLIPPING_INPUT_ASSEMBLY_FA=/other/assembly.fa bash -c '
  config="$1"; set -- /cli/assembly.fa /cli/clipping
  source "$config/lib_config.sh"
  [[ "$ASSEMBLY_FA" == /cli/assembly.fa ]]
  [[ "$CLIPPING_DIR" == /cli/clipping ]]
' _ "$HERE"

# Legacy flag-based circularity interface remains unchanged.
QC_STAGE=02_circularity bash -c '
  config="$1"; set -- --mag-dir /flag/mags --outdir /flag/out
  source "$config/lib_config.sh"
  [[ "$MAG_DIR" == "${MAG_DIR:-}" ]]
  [[ "$CIRCULARITY_DIR" == "$OUT/02_circularity" ]]
' _ "$HERE"

# A single positional argument for a two-argument module is a usage error.
if QC_STAGE=03_clipping bash -c 'config="$1"; set -- /input.fa; source "$config/lib_config.sh"' _ "$HERE" >/dev/null 2>&1; then
    echo 'ERROR: single-argument clipping invocation unexpectedly succeeded' >&2
    exit 1
fi

# The taxonomy helper accepts two arguments, but still requires its DB and MAG map.
if bash "$HERE/07_chimera.sh" __taxonomy_internal /inputs/assembly.fa /outputs/taxonomy \
    >"$TEST_ROOT/taxonomy.log" 2>&1; then
    echo 'ERROR: taxonomy unexpectedly accepted missing required settings' >&2
    exit 1
fi
grep -q 'MMSEQS_TAX_DB and MAG_CONTIG_MAP must be set' "$TEST_ROOT/taxonomy.log"

echo 'Positional input/output configuration tests passed.'
