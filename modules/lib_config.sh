#!/usr/bin/env bash
set -euo pipefail
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QC_STAGE="${QC_STAGE:-$(basename "${BASH_SOURCE[1]:-lib_config.sh}" .sh)}"

QC_CLI_INPUT=""
QC_CLI_OUTPUT=""
if [[ "$QC_STAGE" == "01_mapping" ||
      "$QC_STAGE" == "03_clipping" || "$QC_STAGE" == "07_chimera" ||
      "$QC_STAGE" == "02_circularity" || "$QC_STAGE" == "04_false_closure" ||
      "$QC_STAGE" == "09_assembly_variants" || "$QC_STAGE" == "05_phantom" ||
      "$QC_STAGE" == "06_repeats" || "$QC_STAGE" == "08_false_duplication" ||
      "$QC_STAGE" == "10_summary" ]]; then
    if [[ "$QC_STAGE" != "02_circularity" || "${1:-}" != -* ]]; then
        if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
            case "$QC_STAGE" in
                01_mapping) description='Map PacBio HiFi reads to each original assembly'; input='HiFi FASTQ or FASTQ.GZ'; extras='ASSEMBLER_ASSEMBLY_DIR, ASSEMBLER_FASTAS'; output='OUTPUT_DIR/mapping_manifest.tsv' ;;
            03_clipping) description='Identify clipping and validate read support'; input='All-assemblies FASTA'; extras='MAG_DIR, ASSEMBLER_MAPPING_MANIFEST with assembler HiFi BAMs from 01_mapping; MAG map generated automatically'; output='OUTPUT_DIR/clipping.site.tsv' ;;
            07_chimera) description='Assess genome-wide GTDB taxonomic transitions'; input='GTDB R232 gene-taxonomy TSV'; extras='ASSEMBLY_FA, MAG_CONTIG_MAP, mapping manifest, 03_clipping result'; output='OUTPUT_DIR/chimera.contig.tsv' ;;
            02_circularity) description='Validate end-to-start closure of single-contig MAGs'; input='Directory of MAG FASTA files'; extras='READS, ASSEMBLER_ASSEMBLY_DIR, ASSEMBLER_FASTAS'; output='OUTPUT_DIR/circularity.MAG.tsv' ;;
            04_false_closure) description='Evaluate false-closure candidates'; input='All-assemblies FASTA'; extras='READS, MAG_CONTIG_MAP, mapping, clipping, circularity results'; output='OUTPUT_DIR/false_closure.MAG.tsv' ;;
            09_assembly_variants) description='Quantify assembly allele support using anvi-o'; input='All-assemblies FASTA'; extras='READS, MAG_CONTIG_MAP, mapping manifest, active MAGverify Conda environment'; output='OUTPUT_DIR/variants.contig.tsv' ;;
            05_phantom) description='Assess zero-coverage sequences and read k-mer evidence'; input='All-assemblies FASTA'; extras='READS must be four-line FASTQ(.GZ), MAG_CONTIG_MAP, mapping manifest'; output='OUTPUT_DIR/phantom.region.tsv' ;;
            06_repeats) description='Quantify intra- and inter-contig repeats'; input='All-assemblies FASTA'; extras='MAG_CONTIG_MAP, assembler FASTAs'; output='OUTPUT_DIR/repeats.contig.tsv' ;;
            08_false_duplication) description='Assess candidate false duplications'; input='All-assemblies FASTA'; extras='MAG_CONTIG_MAP, mapping manifest, clipping result'; output='OUTPUT_DIR/false_duplication.contig.tsv' ;;
            10_summary) description='Integrate and validate enabled QC results'; input='All-assemblies FASTA'; extras='MAG_CONTIG_MAP and outputs from all enabled analyses'; output='OUTPUT_DIR/summary.MAG.tsv' ;;
            esac
            printf 'Usage: bash %s INPUT OUTPUT_DIR\n\n%s\nInput: %s\nAdditional inputs/settings: %s\nMain output: %s\n\n' \
                "$QC_STAGE.sh" "$description" "$input" "$extras" "$output"
            printf 'Use MAGverify MODULE --help to specify paths directly.\n' 
            exit 0
        fi
        if { (( $# > 2 )) && [[ "$QC_STAGE" != "02_circularity" ]]; } ||
           { (( $# == 1 )) && [[ "$QC_STAGE" != "07_chimera" ]]; }; then
            printf 'ERROR: usage: bash %s [INPUT OUTPUT_DIR]\n' "$QC_STAGE.sh" >&2
            exit 2
        fi
        QC_CLI_INPUT="${1:-}"
        QC_CLI_OUTPUT="${2:-}"
    fi
fi
if [[ -n "$QC_CLI_OUTPUT" ]]; then
    OUT="${OUT:-$(dirname -- "$QC_CLI_OUTPUT")}"
    case "$QC_STAGE" in
        01_mapping) MAPPING_DIR="$QC_CLI_OUTPUT" ;;
        03_clipping) CLIPPING_DIR="$QC_CLI_OUTPUT" ;;
        07_chimera) CHIMERA_ROOT_DIR="$QC_CLI_OUTPUT" ;;
        02_circularity) CIRCULARITY_DIR="$QC_CLI_OUTPUT"; CIRCULARITY_ROTATED_RESULT_DIR="$QC_CLI_OUTPUT" ;;
        04_false_closure) FALSE_CLOSURE_DIR="$QC_CLI_OUTPUT"; FALSE_CLOSURE_RESULT_DIR="$QC_CLI_OUTPUT" ;;
        09_assembly_variants) VARIANT_DIR="$QC_CLI_OUTPUT" ;;
        05_phantom) PHANTOM_DIR="$QC_CLI_OUTPUT" ;;
        06_repeats) REPEAT_DIR="$QC_CLI_OUTPUT" ;;
        08_false_duplication) FALSE_DUP_DIR="$QC_CLI_OUTPUT" ;;
        10_summary) SUMMARY_DIR="$QC_CLI_OUTPUT" ;;
    esac
fi

source "$CONFIG_DIR/defaults.sh"

MINIMAP2_BIN="${MINIMAP2_BIN:-${MINIMAP2_ENV:+$MINIMAP2_ENV/bin/minimap2}}"
MERYL_BIN="${MERYL_BIN:-${MERYL_ENV:+$MERYL_ENV/bin/meryl}}"
if [[ -n "${MINIMAP2_ENV:-}" && -x "$MINIMAP2_BIN" ]]; then
    export PATH="$MINIMAP2_ENV/bin:$PATH"
fi
if [[ -n "${MERYL_ENV:-}" && -x "$MERYL_BIN" ]]; then
    export PATH="$MERYL_ENV/bin:$PATH"
fi
if [[ -n "${ANVIO_ENV_PREFIX:-}" && -d "$ANVIO_ENV_PREFIX/bin" ]]; then
    export PATH="$ANVIO_ENV_PREFIX/bin:$PATH"
fi

detect_max_threads() {
    local value=""
    for value in \
        "${SLURM_CPUS_PER_TASK:-}" \
        "${SLURM_CPUS_ON_NODE:-}" \
        "${PBS_NP:-}" \
        "${NSLOTS:-}"; do
        if [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
            printf '%s\n' "$value"
            return 0
        fi
    done
    if command -v nproc >/dev/null 2>&1; then
        value="$(nproc 2>/dev/null || true)"
    elif command -v getconf >/dev/null 2>&1; then
        value="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
    fi
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        value=1
    fi
    printf '%s\n' "$value"
}

detect_max_memory_gb() {
    local kb="" pages="" size=""
    if [[ -r /proc/meminfo ]]; then
        kb="$(awk '/^MemTotal:/{print $2; exit}' /proc/meminfo 2>/dev/null || true)"
    fi
    if [[ ! "$kb" =~ ^[1-9][0-9]*$ ]] && command -v getconf >/dev/null 2>&1; then
        pages="$(getconf _PHYS_PAGES 2>/dev/null || true)"
        size="$(getconf PAGE_SIZE 2>/dev/null || true)"
        if [[ "$pages" =~ ^[1-9][0-9]*$ && "$size" =~ ^[1-9][0-9]*$ ]]; then
            kb=$(( pages * size / 1024 ))
        fi
    fi
    [[ "$kb" =~ ^[1-9][0-9]*$ ]] || kb=0
    printf '%s\n' "$(( kb / 1024 / 1024 ))"
}

# Without -t or -m, MAGverify leaves a fifth of the machine free.
RESOURCE_DEFAULT_PERCENT="${RESOURCE_DEFAULT_PERCENT:-80}"

MAX_THREADS="$(detect_max_threads)"
DEFAULT_THREADS=$(( MAX_THREADS * RESOURCE_DEFAULT_PERCENT / 100 ))
(( DEFAULT_THREADS >= 1 )) || DEFAULT_THREADS=1
THREADS="${THREADS:-$DEFAULT_THREADS}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: THREADS must be a positive integer: $THREADS" >&2
    return 1 2>/dev/null || exit 1
}
if (( THREADS > MAX_THREADS )); then
    echo "ERROR: THREADS=$THREADS exceeds the detected allocation of $MAX_THREADS" >&2
    return 1 2>/dev/null || exit 1
fi
MAX_MEMORY_GB="$(detect_max_memory_gb)"
DEFAULT_MEMORY_GB=$(( MAX_MEMORY_GB * RESOURCE_DEFAULT_PERCENT / 100 ))
if (( DEFAULT_MEMORY_GB >= 1 )); then
    MEMORY_GB="${MEMORY_GB:-$DEFAULT_MEMORY_GB}"
else
    MEMORY_GB="${MEMORY_GB:-}"
fi
if [[ -n "$MEMORY_GB" ]]; then
    [[ "$MEMORY_GB" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: MEMORY_GB must be a positive integer: $MEMORY_GB" >&2
        return 1 2>/dev/null || exit 1
    }
    if (( MAX_MEMORY_GB > 0 && MEMORY_GB > MAX_MEMORY_GB )); then
        echo "ERROR: MEMORY_GB=$MEMORY_GB exceeds the detected memory of ${MAX_MEMORY_GB} GB" >&2
        return 1 2>/dev/null || exit 1
    fi
fi

PRESET="map-hifi"

RUN_CLIPPING="${RUN_CLIPPING:-Y}"
RUN_TAXONOMIC_CHIMERA="${RUN_TAXONOMIC_CHIMERA:-Y}"
RUN_CIRCULARITY="${RUN_CIRCULARITY:-Y}"
RUN_FALSE_CLOSURE="${RUN_FALSE_CLOSURE:-Y}"
RUN_ASSEMBLY_VARIANTS="${RUN_ASSEMBLY_VARIANTS:-Y}"
RUN_PHANTOM="${RUN_PHANTOM:-Y}"
RUN_REPEAT="${RUN_REPEAT:-Y}"
RUN_FALSE_DUPLICATION="${RUN_FALSE_DUPLICATION:-Y}"

ASSEMBLY_FA="${ASSEMBLY_FA:-}"

READS="${READS:-}"

MAG_DIR="${MAG_DIR:-}"
MAG_EXTENSIONS="fa fasta fna"

ASSEMBLER_ASSEMBLY_DIR="${ASSEMBLER_ASSEMBLY_DIR:-}"
ASSEMBLER_FASTAS="${ASSEMBLER_FASTAS:-hifiasm_meta.fa metaflye.fa metaMDBG.fa}"

if [[ -z "${OUT:-}" ]]; then
    echo "ERROR: output directory is not set; run modules through MAGverify with -o (or set OUT)" >&2
    exit 1
fi

# Stage input/output paths (previously stage_paths.sh).
MAPPING_INPUT_READS="${MAPPING_INPUT_READS:-}"
MAPPING_INPUT_ASSEMBLER_ASSEMBLY_DIR="${MAPPING_INPUT_ASSEMBLER_ASSEMBLY_DIR:-}"
MAPPING_INPUT_ASSEMBLER_FASTAS="${MAPPING_INPUT_ASSEMBLER_FASTAS:-}"

MAPPING_DIR="${MAPPING_DIR:-$OUT/01_mapping}"
ASSEMBLER_MAPPING_MANIFEST="${ASSEMBLER_MAPPING_MANIFEST:-$MAPPING_DIR/.cache/mapping_manifest.tsv}"

CIRCULARITY_INPUT_MAG_DIR="${CIRCULARITY_INPUT_MAG_DIR:-}"
CIRCULARITY_INPUT_READS="${CIRCULARITY_INPUT_READS:-}"
CIRCULARITY_INPUT_ASSEMBLER_ASSEMBLY_DIR="${CIRCULARITY_INPUT_ASSEMBLER_ASSEMBLY_DIR:-}"
CIRCULARITY_INPUT_ASSEMBLER_FASTAS="${CIRCULARITY_INPUT_ASSEMBLER_FASTAS:-}"
CIRCULARITY_INPUT_CIRCULARITY_TSV="${CIRCULARITY_INPUT_CIRCULARITY_TSV:-}"

CIRCULARITY_DIR="${CIRCULARITY_DIR:-$OUT/02_circularity}"

CIRCULARITY_ONLY_FINAL="${CIRCULARITY_ONLY_FINAL:-}"

CLIPPING_INPUT_ASSEMBLY_FA="${CLIPPING_INPUT_ASSEMBLY_FA:-}"
CLIPPING_INPUT_MAG_DIR="${CLIPPING_INPUT_MAG_DIR:-}"
CLIPPING_INPUT_READS="${CLIPPING_INPUT_READS:-}"
CLIPPING_INPUT_MAG_CONTIG_MAP="${CLIPPING_INPUT_MAG_CONTIG_MAP:-}"
CLIPPING_INPUT_ASSEMBLER_MAPPING_MANIFEST="${CLIPPING_INPUT_ASSEMBLER_MAPPING_MANIFEST:-}"

CLIPPING_DIR="${CLIPPING_DIR:-$OUT/03_clipping}"

CLIP_FINAL="${CLIP_FINAL:-}"
CLIP_ALL_FINAL="${CLIP_ALL_FINAL:-}"

FALSE_CLOSURE_INPUT_READS="${FALSE_CLOSURE_INPUT_READS:-}"
FALSE_CLOSURE_INPUT_ASSEMBLY_FA="${FALSE_CLOSURE_INPUT_ASSEMBLY_FA:-}"
FALSE_CLOSURE_INPUT_MAG_CONTIG_MAP="${FALSE_CLOSURE_INPUT_MAG_CONTIG_MAP:-}"
FALSE_CLOSURE_INPUT_ASSEMBLER_MAPPING_MANIFEST="${FALSE_CLOSURE_INPUT_ASSEMBLER_MAPPING_MANIFEST:-}"
FALSE_CLOSURE_INPUT_CLIP_FINAL="${FALSE_CLOSURE_INPUT_CLIP_FINAL:-}"
FALSE_CLOSURE_INPUT_CLIP_ALL_FINAL="${FALSE_CLOSURE_INPUT_CLIP_ALL_FINAL:-}"
FALSE_CLOSURE_INPUT_CIRCULARITY_ONLY_FINAL="${FALSE_CLOSURE_INPUT_CIRCULARITY_ONLY_FINAL:-}"
FALSE_CLOSURE_INPUT_CIRCULARITY_TSV="${FALSE_CLOSURE_INPUT_CIRCULARITY_TSV:-}"
FALSE_CLOSURE_INPUT_SCG_TSV="${FALSE_CLOSURE_INPUT_SCG_TSV:-}"

FALSE_CLOSURE_DIR="${FALSE_CLOSURE_DIR:-$OUT/04_false_closure}"

FALSE_CLOSURE_FINAL="${FALSE_CLOSURE_FINAL:-}"

PHANTOM_INPUT_ASSEMBLY_FA="${PHANTOM_INPUT_ASSEMBLY_FA:-}"
PHANTOM_INPUT_READS="${PHANTOM_INPUT_READS:-}"
PHANTOM_INPUT_MAG_CONTIG_MAP="${PHANTOM_INPUT_MAG_CONTIG_MAP:-}"
PHANTOM_INPUT_ASSEMBLER_MAPPING_MANIFEST="${PHANTOM_INPUT_ASSEMBLER_MAPPING_MANIFEST:-}"

PHANTOM_DIR="${PHANTOM_DIR:-$OUT/05_phantom}"

PHANTOM_FINAL="${PHANTOM_FINAL:-}"

REPEAT_INPUT_ASSEMBLY_FA="${REPEAT_INPUT_ASSEMBLY_FA:-}"
REPEAT_INPUT_MAG_CONTIG_MAP="${REPEAT_INPUT_MAG_CONTIG_MAP:-}"
REPEAT_INPUT_ASSEMBLER_ASSEMBLY_DIR="${REPEAT_INPUT_ASSEMBLER_ASSEMBLY_DIR:-}"
REPEAT_INPUT_ASSEMBLER_FASTAS="${REPEAT_INPUT_ASSEMBLER_FASTAS:-}"

REPEAT_DIR="${REPEAT_DIR:-$OUT/06_repeats}"

REPEAT_SUM="${REPEAT_SUM:-}"

CHIMERA_INPUT_ASSEMBLY_FA="${CHIMERA_INPUT_ASSEMBLY_FA:-}"
CHIMERA_INPUT_READS="${CHIMERA_INPUT_READS:-}"
CHIMERA_INPUT_MAG_CONTIG_MAP="${CHIMERA_INPUT_MAG_CONTIG_MAP:-}"
CHIMERA_INPUT_ASSEMBLER_MAPPING_MANIFEST="${CHIMERA_INPUT_ASSEMBLER_MAPPING_MANIFEST:-}"
CHIMERA_INPUT_GTDB_GENE_TAX_TSV="${CHIMERA_INPUT_GTDB_GENE_TAX_TSV:-}"
CHIMERA_INPUT_MAG_GTDB_TAXONOMY_TSV="${CHIMERA_INPUT_MAG_GTDB_TAXONOMY_TSV:-}"
CHIMERA_INPUT_CLIP_FINAL="${CHIMERA_INPUT_CLIP_FINAL:-}"
CHIMERA_INPUT_SCG_TSV="${CHIMERA_INPUT_SCG_TSV:-}"

CHIMERA_ROOT_DIR="${CHIMERA_ROOT_DIR:-$OUT/07_chimera}"

GTDB_CHIMERA_CONTIG="${GTDB_CHIMERA_CONTIG:-}"
GTDB_CHIMERA_MAG="${GTDB_CHIMERA_MAG:-}"

FALSE_DUP_INPUT_ASSEMBLY_FA="${FALSE_DUP_INPUT_ASSEMBLY_FA:-}"
FALSE_DUP_INPUT_READS="${FALSE_DUP_INPUT_READS:-}"
FALSE_DUP_INPUT_MAG_CONTIG_MAP="${FALSE_DUP_INPUT_MAG_CONTIG_MAP:-}"
FALSE_DUP_INPUT_ASSEMBLER_MAPPING_MANIFEST="${FALSE_DUP_INPUT_ASSEMBLER_MAPPING_MANIFEST:-}"
FALSE_DUP_INPUT_CLIP_FINAL="${FALSE_DUP_INPUT_CLIP_FINAL:-}"

FALSE_DUP_DIR="${FALSE_DUP_DIR:-$OUT/08_false_duplication}"

FALSE_DUP_CONTIG_SUM="${FALSE_DUP_CONTIG_SUM:-}"

VARIANT_INPUT_ASSEMBLY_FA="${VARIANT_INPUT_ASSEMBLY_FA:-}"
VARIANT_INPUT_READS="${VARIANT_INPUT_READS:-}"
VARIANT_INPUT_MAG_CONTIG_MAP="${VARIANT_INPUT_MAG_CONTIG_MAP:-}"
VARIANT_INPUT_ASSEMBLER_MAPPING_MANIFEST="${VARIANT_INPUT_ASSEMBLER_MAPPING_MANIFEST:-}"

VARIANT_DIR="${VARIANT_DIR:-$OUT/09_assembly_variants}"

ALLELE_CONTIG="${ALLELE_CONTIG:-}"
ALLELE_MAG="${ALLELE_MAG:-}"

SUMMARY_INPUT_ASSEMBLY_FA="${SUMMARY_INPUT_ASSEMBLY_FA:-}"
SUMMARY_INPUT_MAG_CONTIG_MAP="${SUMMARY_INPUT_MAG_CONTIG_MAP:-}"
SUMMARY_INPUT_CLIP_FINAL="${SUMMARY_INPUT_CLIP_FINAL:-}"
SUMMARY_INPUT_ALLELE_CONTIG="${SUMMARY_INPUT_ALLELE_CONTIG:-}"
SUMMARY_INPUT_PHANTOM_FINAL="${SUMMARY_INPUT_PHANTOM_FINAL:-}"
SUMMARY_INPUT_REPEAT_SUM="${SUMMARY_INPUT_REPEAT_SUM:-}"
SUMMARY_INPUT_FALSE_DUP_CONTIG_SUM="${SUMMARY_INPUT_FALSE_DUP_CONTIG_SUM:-}"
SUMMARY_INPUT_CIRCULARITY_ONLY_FINAL="${SUMMARY_INPUT_CIRCULARITY_ONLY_FINAL:-}"
SUMMARY_INPUT_FALSE_CLOSURE_FINAL="${SUMMARY_INPUT_FALSE_CLOSURE_FINAL:-}"

SUMMARY_INPUT_REPEAT_MAG_SUM="${SUMMARY_INPUT_REPEAT_MAG_SUM:-}"
SUMMARY_INPUT_FALSE_DUP_CLIPPING_CONTIG_SUM="${SUMMARY_INPUT_FALSE_DUP_CLIPPING_CONTIG_SUM:-}"
SUMMARY_INPUT_FALSE_DUP_CLIPPING_MAG_SUM="${SUMMARY_INPUT_FALSE_DUP_CLIPPING_MAG_SUM:-}"
SUMMARY_INPUT_CIRCULARITY_ONLY_AUDIT="${SUMMARY_INPUT_CIRCULARITY_ONLY_AUDIT:-}"
SUMMARY_INPUT_FALSE_CLOSURE_AUDIT="${SUMMARY_INPUT_FALSE_CLOSURE_AUDIT:-}"
SUMMARY_INPUT_GTDB_CHIMERA_REPORT_CONTIG="${SUMMARY_INPUT_GTDB_CHIMERA_REPORT_CONTIG:-}"

SUMMARY_DIR="${SUMMARY_DIR:-$OUT/.cache/summary}"

SUMMARY_CONTIG="${SUMMARY_CONTIG:-}"
SUMMARY_MAG_CONTIG="${SUMMARY_MAG_CONTIG:-}"
SUMMARY_MAG="${SUMMARY_MAG:-}"
QC_STAGE="${QC_STAGE:-$(basename "${BASH_SOURCE[1]:-lib_config.sh}" .sh)}"
case "$QC_STAGE" in
    01_mapping) QC_STAGE_PREFIX=MAPPING ;;
    03_clipping) QC_STAGE_PREFIX=CLIPPING ;;
    07_chimera) QC_STAGE_PREFIX=CHIMERA ;;
    02_circularity) QC_STAGE_PREFIX=CIRCULARITY ;;
    04_false_closure) QC_STAGE_PREFIX=FALSE_CLOSURE ;;
    09_assembly_variants) QC_STAGE_PREFIX=VARIANT ;;
    05_phantom) QC_STAGE_PREFIX=PHANTOM ;;
    06_repeats) QC_STAGE_PREFIX=REPEAT ;;
    08_false_duplication) QC_STAGE_PREFIX=FALSE_DUP ;;
    10_summary) QC_STAGE_PREFIX=SUMMARY ;;
    *) QC_STAGE_PREFIX="" ;;
esac
if [[ -n "$QC_STAGE_PREFIX" ]]; then
    for input_name in ASSEMBLY_FA READS MAG_DIR MAG_CONTIG_MAP \
                      ASSEMBLER_ASSEMBLY_DIR ASSEMBLER_FASTAS \
                      ASSEMBLER_MAPPING_MANIFEST CIRCULARITY_TSV \
                      GTDB_GENE_TAX_TSV MAG_GTDB_TAXONOMY_TSV \
                      SCG_TSV CLIP_FINAL CLIP_ALL_FINAL \
                      ALLELE_CONTIG PHANTOM_FINAL REPEAT_SUM REPEAT_MAG_SUM \
                      FALSE_DUP_CONTIG_SUM FALSE_DUP_CLIPPING_CONTIG_SUM \
                      FALSE_DUP_CLIPPING_MAG_SUM CIRCULARITY_ONLY_FINAL \
                      CIRCULARITY_ONLY_AUDIT FALSE_CLOSURE_FINAL \
                      FALSE_CLOSURE_AUDIT GTDB_CHIMERA_REPORT_CONTIG; do
        override_name="${QC_STAGE_PREFIX}_INPUT_${input_name}"
        if [[ -n "${!override_name:-}" ]]; then
            printf -v "$input_name" '%s' "${!override_name}"
        fi
    done
fi

if [[ -n "$QC_CLI_INPUT" ]]; then
    case "$QC_STAGE" in
        01_mapping) READS="$QC_CLI_INPUT" ;;
        07_chimera) GTDB_GENE_TAX_TSV="$QC_CLI_INPUT" ;;
        02_circularity) MAG_DIR="$QC_CLI_INPUT" ;;
        *) ASSEMBLY_FA="$QC_CLI_INPUT" ;;
    esac
fi

PREP_DIR="${PREP_DIR:-$OUT/00_prep}"
MAG_CONTIG_MAP="${MAG_CONTIG_MAP:-$PREP_DIR/MAG_contig_map.tsv}"

CIRCULARITY_TSV="${CIRCULARITY_TSV:-}"

GENE_CALLS_DIR="${GENE_CALLS_DIR:-$PREP_DIR}"
SCG_OUT_DIR="${SCG_OUT_DIR:-$PREP_DIR/.cache/scg}"
SCG_TSV="${SCG_TSV:-$PREP_DIR/scg.contig.tsv}"

ANVIO_DEFAULT_CLIPPING_RATIO=0.80
MIN_DIST_TO_END=0
TERMINAL_CLIPPING_WINDOW_BP=100

STRINGENT_CLIP_MIN_DEPTH=10
STRINGENT_CLIP_RATIO=1.0

ANVIO_CLIPPING_ALGORITHM_ID="anvi_script_find_misassemblies_source_assembler_all_mapped_end0_v7"
CLIPPING_NORMALIZATION_ALGORITHM_ID="anvio_exact_default_and_stringent_source_specific_v2"
CLIPPING_VALIDATION_ALGORITHM_ID="qname_source_assembler_all_mapped_boundary_safe_exact_idn_length_audit_v9"

SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-$THREADS}"
CLIPPING_THREADS="${CLIPPING_THREADS:-$THREADS}"
SORT_MEM="2G"

FORCE_SINGLE_PART_MINIMAP2="Y"
MINIMAP2_INDEX_CHUNK_SINGLE="900G"
MINIMAP2_INDEX_CHUNK="256G"
MINIMAP2_MAX_SECONDARY=5

INDEX_CSI="Y"

FORCE_ASSEMBLY_REMAP="${FORCE_ASSEMBLY_REMAP:-N}"

ANVIO_ENV_PREFIX="${ANVIO_ENV_PREFIX:-${CONDA_PREFIX:-}}"
if [[ -z "$ANVIO_ENV_PREFIX" ]] && command -v anvi-profile >/dev/null 2>&1; then
    ANVIO_ENV_PREFIX="$(cd "$(dirname "$(command -v anvi-profile)")/.." && pwd -P)"
fi
ANVIO_THREADS="${ANVIO_THREADS:-$THREADS}"

FORCE_RERUN_CLIPPING="${FORCE_RERUN_CLIPPING:-N}"
CLIPPING_POSTPROCESS_ONLY="${CLIPPING_POSTPROCESS_ONLY:-N}"
CLIPPING_MAG_SUMMARY_ONLY="${CLIPPING_MAG_SUMMARY_ONLY:-N}"

CIRCULARITY_SCRIPT="$CONFIG_DIR/02_circularity.sh"
CIRCULARITY_ROTATED_RESULT_DIR="${CIRCULARITY_ROTATED_RESULT_DIR:-$CIRCULARITY_DIR}"

CIRCULARITY_HIFI_ANCHOR_BP=2000
CIRCULARITY_MIN_HIFI_QNAMES=2

CIRCULARITY_THREADS="${CIRCULARITY_THREADS:-$THREADS}"
CIRCULARITY_SORT_THREADS="${CIRCULARITY_SORT_THREADS:-$(( THREADS < 20 ? THREADS : 20 ))}"
CIRCULARITY_SORT_MEMORY="2G"
CIRCULARITY_MINIMAP2_INDEX_CHUNK="$MINIMAP2_INDEX_CHUNK_SINGLE"
CIRCULARITY_MINIMAP2_MAX_SECONDARY="$MINIMAP2_MAX_SECONDARY"

PREMATURE_MAX_LENGTH_BP=500000
PREMATURE_MIN_RIBOSOMAL_FAMILIES=3
FALSE_CLOSURE_MIN_LINKED_READS_PER_SIDE=2

ASSEMBLY_VARIANT_SAMPLE_ID="${ASSEMBLY_VARIANT_SAMPLE_ID:-HiFi_reads}"
FORCE_RERUN_ASSEMBLY_VARIANTS="${FORCE_RERUN_ASSEMBLY_VARIANTS:-N}"

ASSEMBLY_VARIANT_MEMORY_BUDGET_GB="${ASSEMBLY_VARIANT_MEMORY_BUDGET_GB:-${MEMORY_GB:-900}}"
ASSEMBLY_VARIANT_MEMORY_RESERVE_GB="${ASSEMBLY_VARIANT_MEMORY_RESERVE_GB:-64}"
ALLELE_ANVIO_MEMORY_PER_THREAD_GB="${ALLELE_ANVIO_MEMORY_PER_THREAD_GB:-16}"
ALLELE_ANVIO_MEMORY_POLL_SECONDS="${ALLELE_ANVIO_MEMORY_POLL_SECONDS:-5}"
ALLELE_ANVIO_FULL_PROFILE_THREADS="${ALLELE_ANVIO_FULL_PROFILE_THREADS:-$THREADS}"
ALLELE_ANVIO_PREPROCESS_THREADS="${ALLELE_ANVIO_PREPROCESS_THREADS:-$THREADS}"
ASSEMBLY_VARIANT_KEEP_FULL_PROFILE="${ASSEMBLY_VARIANT_KEEP_FULL_PROFILE:-N}"
ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES="${ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES:-Y}"
ASSEMBLY_VARIANT_DROP_FILE_CACHE="${ASSEMBLY_VARIANT_DROP_FILE_CACHE:-Y}"
ASSEMBLY_VARIANT_SKIP_GENE_CALLING="${ASSEMBLY_VARIANT_SKIP_GENE_CALLING:-N}"

ZERO_MIN_BP=501
PHANTOM_BLAST_THREADS="$THREADS"
PHANTOM_KMER_K=21
FORCE_RERUN_PHANTOM="N"

REPEAT_MIN_ALIGN_PAPER=200
REPEAT_PIDENT_PAPER=80
REPEAT_CONTIG_HIGH=0.70
REPEAT_GENOME_HIGH=0.70
REPEAT_HIGH="$REPEAT_CONTIG_HIGH"
REPEAT_CHUNK_BP=10000000
REPEAT_BLAST_THREADS_PER_JOB=1
REPEAT_JOBS="${REPEAT_JOBS:-$THREADS}"
(( REPEAT_JOBS >= 1 )) || REPEAT_JOBS=1
REPEAT_SORT_THREADS="$THREADS"
REPEAT_SORT_MEMORY="50%"
REPEAT_MIN_CONTIG_LEN=0
FORCE_RERUN_REPEAT="N"

FORCE_RERUN_FALSE_DUPLICATION="N"
FALSE_DUP_MAPPING_MANIFEST="${FALSE_DUP_MAPPING_MANIFEST:-${ASSEMBLER_MAPPING_MANIFEST:-$MAPPING_DIR/.cache/mapping_manifest.tsv}}"
FALSE_DUP_MIN_COPY_BP="${FALSE_DUP_MIN_COPY_BP:-1000}"
FALSE_DUP_MIN_PIDENT="${FALSE_DUP_MIN_PIDENT:-95}"
FALSE_DUP_HIFI_ANCHOR_BP=500
FALSE_DUP_HIFI_MIN_SPANNING_READS="${FALSE_DUP_HIFI_MIN_SPANNING_READS:-2}"
FALSE_DUP_HIFI_MIN_REACHING_READS="${FALSE_DUP_HIFI_MIN_REACHING_READS:-2}"

if [[ -z "${ANVIO_FIND_MISASSEMBLIES:-}" ]]; then
    if [[ -n "$ANVIO_ENV_PREFIX" && -x "$ANVIO_ENV_PREFIX/bin/anvi-script-find-misassemblies" ]]; then
        ANVIO_FIND_MISASSEMBLIES="$ANVIO_ENV_PREFIX/bin/anvi-script-find-misassemblies"
    else
        ANVIO_FIND_MISASSEMBLIES="$(command -v anvi-script-find-misassemblies || true)"
    fi
fi
CLIP_ALL_FINAL="${CLIP_ALL_FINAL:-${CLIPPING_DIR}/.cache/clipping.site.tsv}"
CLIP_CONTIG_SUMMARY="${CLIP_CONTIG_SUMMARY:-${CLIPPING_DIR}/.cache/clipping.contig.tsv}"
CLIP_MAG_SUMMARY="${CLIP_MAG_SUMMARY:-${CLIPPING_DIR}/.cache/clipping.MAG.tsv}"
CLIP_ASSEMBLER_SUMMARY="${CLIP_ASSEMBLER_SUMMARY:-${CLIPPING_DIR}/.cache/clipping.assembler.tsv}"
CLIP_FINAL="${CLIP_FINAL:-${CLIPPING_DIR}/.cache/clipping.site.internal.tsv}"
CLIP_FINAL_CACHE_META="${CLIP_FINAL_CACHE_META:-${CLIPPING_DIR}/.cache/clipping.site.internal.tsv.inputs.sha256}"

CLIP_VALIDATE_SCRIPT="${CLIP_VALIDATE_SCRIPT:-$CONFIG_DIR/../scripts/clipping_validate.py}"
CLIP_SUMMARIZE_SCRIPT="${CLIP_SUMMARIZE_SCRIPT:-$CONFIG_DIR/../scripts/clipping_summarize.py}"
MAG_MAP_BUILDER="${MAG_MAP_BUILDER:-$CONFIG_DIR/../scripts/prep_build_mag_contig_map.py}"
ASSEMBLER_MAPPING_MANIFEST="${ASSEMBLER_MAPPING_MANIFEST:-${MAPPING_DIR}/.cache/mapping_manifest.tsv}"

GTDB_RELEASE="${GTDB_RELEASE:-GTDB_R232}"
GTDB_GENE_TAX_DIR="${GTDB_GENE_TAX_DIR:-$PREP_DIR/.cache/gene_taxonomy}"
GTDB_GENE_TAX_TSV="${GTDB_GENE_TAX_TSV:-$PREP_DIR/taxonomy.gene.tsv}"
GTDB_CHIMERA_DIR="${GTDB_CHIMERA_DIR:-${CHIMERA_ROOT_DIR}/.cache}"
GTDB_CHIMERA_CONTIG="${GTDB_CHIMERA_CONTIG:-$GTDB_CHIMERA_DIR/contig_structural_chimera.${GTDB_RELEASE}.final.tsv}"
GTDB_CHIMERA_MAG_CONTIG="${GTDB_CHIMERA_MAG_CONTIG:-$GTDB_CHIMERA_DIR/MAG_contig_structural_chimera.${GTDB_RELEASE}.final.tsv}"
GTDB_CHIMERA_MAG="${GTDB_CHIMERA_MAG:-$GTDB_CHIMERA_DIR/MAG_structural_chimera.${GTDB_RELEASE}.final.tsv}"
MAG_GTDB_TAXONOMY_TSV="${MAG_GTDB_TAXONOMY_TSV:-$PREP_DIR/gtdbtk_taxonomy.MAG.tsv}"

CIRCULARITY_RESULT_DIR="${CIRCULARITY_RESULT_DIR:-$CIRCULARITY_DIR}"
FALSE_CLOSURE_RESULT_DIR="${FALSE_CLOSURE_RESULT_DIR:-$FALSE_CLOSURE_DIR}"

CIRCULARITY_ONLY_FINAL="${CIRCULARITY_ONLY_FINAL:-$CIRCULARITY_RESULT_DIR/.cache/circularity.MAG.tsv}"
CIRCULARITY_ONLY_AUDIT="${CIRCULARITY_ONLY_AUDIT:-$CIRCULARITY_RESULT_DIR/.cache/circularity.MAG.audit.tsv}"
CIRCULARITY_RUN_METADATA="${CIRCULARITY_RUN_METADATA:-$CIRCULARITY_RESULT_DIR/run_metadata.txt}"

FALSE_CLOSURE_FINAL="${FALSE_CLOSURE_FINAL:-$FALSE_CLOSURE_RESULT_DIR/.cache/false_closure.MAG.tsv}"
FALSE_CLOSURE_AUDIT="${FALSE_CLOSURE_AUDIT:-$FALSE_CLOSURE_RESULT_DIR/.cache/false_closure.MAG.audit.tsv}"
FALSE_CLOSURE_CONTIG_FINAL="${FALSE_CLOSURE_CONTIG_FINAL:-$FALSE_CLOSURE_RESULT_DIR/.cache/false_closure.contig.tsv}"
FALSE_CLOSURE_RUN_METADATA="${FALSE_CLOSURE_RUN_METADATA:-$FALSE_CLOSURE_RESULT_DIR/.cache/run_metadata.txt}"
FALSE_CLOSURE_CANDIDATE_READS="${FALSE_CLOSURE_CANDIDATE_READS:-$FALSE_CLOSURE_RESULT_DIR/.cache/false_closure.candidate_reads.tsv}"
FALSE_CLOSURE_EVIDENCE_TSV="${FALSE_CLOSURE_EVIDENCE_TSV:-$FALSE_CLOSURE_RESULT_DIR/.cache/false_closure.evidence.tsv}"
FALSE_CLOSURE_BREAKPOINTS_TSV="${FALSE_CLOSURE_BREAKPOINTS_TSV:-$FALSE_CLOSURE_RESULT_DIR/.cache/false_closure.breakpoint.tsv}"

ALLELE_ANVIO_DIR="${ALLELE_ANVIO_DIR:-${VARIANT_DIR}/.cache/hifi_all_mapped/anvio}"
ALLELE_BY_ASSEMBLER_DIR="${ALLELE_BY_ASSEMBLER_DIR:-$ALLELE_ANVIO_DIR/by_assembler}"
ALLELE_PROFILE_MANIFEST="${ALLELE_PROFILE_MANIFEST:-${VARIANT_DIR}/.cache/anvio_profiles.by_assembler.tsv}"
ALLELE_VARIANTS="${ALLELE_VARIANTS:-${VARIANT_DIR}/.cache/variants.site.tsv}"
ALLELE_CONTIG_ALL="${ALLELE_CONTIG_ALL:-${VARIANT_DIR}/.cache/contig_assembly_allele_support.all_contigs.tsv}"
ALLELE_SUMMARY_ALL="${ALLELE_SUMMARY_ALL:-${VARIANT_DIR}/.cache/assembly_allele_support.all_contigs.summary.tsv}"
ALLELE_ASSEMBLER_SUMMARY="${ALLELE_ASSEMBLER_SUMMARY:-${VARIANT_DIR}/.cache/assembly_allele_support.by_assembler.tsv}"
ALLELE_CONTIG="${ALLELE_CONTIG:-${VARIANT_DIR}/.cache/variants.contig.tsv}"
ALLELE_SUMMARY="${ALLELE_SUMMARY:-${VARIANT_DIR}/.cache/assembly_allele_support.summary.tsv}"
ALLELE_VARIANTS_MAG="${ALLELE_VARIANTS_MAG:-${VARIANT_DIR}/.cache/assembly_allele_disagreements.by_MAG.tsv}"
ALLELE_MAG="${ALLELE_MAG:-${VARIANT_DIR}/.cache/variants.MAG.tsv}"
ALLELE_CACHE_META="${ALLELE_CACHE_META:-${VARIANT_DIR}/.cache/variants.inputs.sha256}"

PHANTOM_CACHE_DIR="${PHANTOM_CACHE_DIR:-${PHANTOM_DIR}/.cache}"
PHANTOM_ZERO_BED="${PHANTOM_ZERO_BED:-${PHANTOM_CACHE_DIR}/zero_regions.bed}"
PHANTOM_ZERO_FA="${PHANTOM_ZERO_FA:-${PHANTOM_CACHE_DIR}/zero_regions.fa}"
PHANTOM_ZERO_TSV="${PHANTOM_ZERO_TSV:-${PHANTOM_CACHE_DIR}/zero_regions.tsv}"
PHANTOM_BLAST_DB="${PHANTOM_BLAST_DB:-${PHANTOM_CACHE_DIR}/raw_reads_blastdb/reads}"
PHANTOM_BLAST="${PHANTOM_BLAST:-${PHANTOM_CACHE_DIR}/zero_regions_vs_raw_reads.tsv}"
PHANTOM_FINAL="${PHANTOM_FINAL:-${PHANTOM_CACHE_DIR}/phantom.region.tsv}"
PHANTOM_CACHE_META="${PHANTOM_CACHE_META:-${PHANTOM_CACHE_DIR}/phantom.region.tsv.inputs.sha256}"
PHANTOM_MAG_SITES="${PHANTOM_MAG_SITES:-${PHANTOM_CACHE_DIR}/phantom.MAG_region.tsv}"
PHANTOM_MAG_SUMMARY="${PHANTOM_MAG_SUMMARY:-${PHANTOM_CACHE_DIR}/phantom.MAG.base.tsv}"
PHANTOM_MERYL_READ_DB="${PHANTOM_MERYL_READ_DB:-${PHANTOM_CACHE_DIR}/meryl/raw_reads.k21.meryl}"
PHANTOM_MERYL_READ_CACHE="${PHANTOM_MERYL_READ_CACHE:-${PHANTOM_MERYL_READ_DB}.inputs.sha256}"
PHANTOM_KMER_FINAL="${PHANTOM_KMER_FINAL:-${PHANTOM_CACHE_DIR}/phantom.kmer.tsv}"
PHANTOM_KMER_CACHE_META="${PHANTOM_KMER_CACHE_META:-${PHANTOM_KMER_FINAL}.inputs.sha256}"
PHANTOM_MAG_MERYL_FINAL="${PHANTOM_MAG_MERYL_FINAL:-${PHANTOM_CACHE_DIR}/phantom.MAG.tsv}"
PHANTOM_MAG_MERYL_CATALOGUE="${PHANTOM_MAG_MERYL_CATALOGUE:-${PHANTOM_CACHE_DIR}/phantom.MAG_catalogue.tsv}"
PHANTOM_MAG_MERYL_CACHE_META="${PHANTOM_MAG_MERYL_CACHE_META:-${PHANTOM_CACHE_DIR}/phantom.MAG.tsv.inputs.sha256}"
PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META="${PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META:-${PHANTOM_MAG_MERYL_CATALOGUE}.inputs.sha256}"
PHANTOM_DATASET="${PHANTOM_DATASET:-${PHANTOM_CACHE_DIR}/phantom.dataset.tsv}"
PHANTOM_CONTIG="${PHANTOM_CONTIG:-${PHANTOM_CACHE_DIR}/phantom.contig.tsv}"

REPEAT_SUM="${REPEAT_SUM:-${REPEAT_DIR}/.cache/repeats.contig.tsv}"
REPEAT_CACHE_META="${REPEAT_CACHE_META:-${REPEAT_DIR}/.cache/repeats.contig.tsv.inputs.sha256}"
REPEAT_MAG_SUM="${REPEAT_MAG_SUM:-${REPEAT_DIR}/.cache/repeats.MAG.tsv}"
REPEAT_DATASET_SUM="${REPEAT_DATASET_SUM:-${REPEAT_DIR}/.cache/repeats.dataset.tsv}"

FALSE_DUP_CACHE_DIR="${FALSE_DUP_CACHE_DIR:-${FALSE_DUP_DIR}/.cache}"
FALSE_DUP_CANDIDATE_PAF="${FALSE_DUP_CANDIDATE_PAF:-${FALSE_DUP_CACHE_DIR}/duplication_self_minimap2.paf.gz}"
FALSE_DUP_CANDIDATE_PROVENANCE="${FALSE_DUP_CANDIDATE_PROVENANCE:-${FALSE_DUP_CACHE_DIR}/duplication_candidate_search.provenance.tsv}"
FALSE_DUP_ALL_PAIR_SUM="${FALSE_DUP_ALL_PAIR_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication_pairs.all_contigs.final.tsv}"
FALSE_DUP_ALL_BOUNDARY_SUM="${FALSE_DUP_ALL_BOUNDARY_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication_boundaries.all_contigs.final.tsv}"
FALSE_DUP_ALL_CONTIG_SUM="${FALSE_DUP_ALL_CONTIG_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication.all_contigs.tsv}"
FALSE_DUP_ASSEMBLER_SUM="${FALSE_DUP_ASSEMBLER_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication.by_assembler.tsv}"
FALSE_DUP_PAIR_SUM="${FALSE_DUP_PAIR_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication.pair.tsv}"
FALSE_DUP_BOUNDARY_SUM="${FALSE_DUP_BOUNDARY_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication_boundaries.final.tsv}"
FALSE_DUP_CONTIG_SUM="${FALSE_DUP_CONTIG_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication.contig.tsv}"
FALSE_DUP_MAG_SUM="${FALSE_DUP_MAG_SUM:-${FALSE_DUP_CACHE_DIR}/false_duplication.MAG.tsv}"
FALSE_DUP_PROVENANCE="${FALSE_DUP_PROVENANCE:-${FALSE_DUP_CACHE_DIR}/false_duplication.provenance.tsv}"
FALSE_DUP_CLIPPING_BOUNDARY_SUM="${FALSE_DUP_CLIPPING_BOUNDARY_SUM:-${FALSE_DUP_CACHE_DIR}/clipping/false_duplication.boundaries.clipping.tsv}"
FALSE_DUP_CLIPPING_PAIR_SUM="${FALSE_DUP_CLIPPING_PAIR_SUM:-${FALSE_DUP_CACHE_DIR}/clipping/false_duplication.pairs.clipping.tsv}"
FALSE_DUP_CLIPPING_CONTIG_SUM="${FALSE_DUP_CLIPPING_CONTIG_SUM:-${FALSE_DUP_CACHE_DIR}/clipping/false_duplication.contig.clipping.tsv}"
FALSE_DUP_CLIPPING_MAG_SUM="${FALSE_DUP_CLIPPING_MAG_SUM:-${FALSE_DUP_CACHE_DIR}/clipping/false_duplication.MAG.clipping.tsv}"
FALSE_DUP_CLIPPING_PROVENANCE="${FALSE_DUP_CLIPPING_PROVENANCE:-${FALSE_DUP_CACHE_DIR}/clipping/false_duplication.clipping.provenance.tsv}"

SUMMARY_CONTIG="${SUMMARY_CONTIG:-${SUMMARY_DIR}/summary.contig.tsv}"
SUMMARY_MAG="${SUMMARY_MAG:-${SUMMARY_DIR}/summary.MAG.tsv}"

mkdir -p "$OUT/00_log" "$PREP_DIR"


# QC helper functions (previously lib_qc.sh).
qc_sq_count() {
    local bam="$1"
    samtools view -H "$bam" 2>/dev/null \
        | awk 'BEGIN{n=0} /^@SQ/{n++} END{print n+0}'
}

qc_bam_is_valid() {
    local bam="$1"
    local n_sq
    [[ -s "$bam" ]] || return 1
    samtools quickcheck "$bam" >/dev/null 2>&1 || return 1
    n_sq="$(qc_sq_count "$bam")"
    [[ "$n_sq" =~ ^[0-9]+$ ]] && (( n_sq > 0 ))
}

qc_bam_matches_fai() {
    local bam="$1"
    local fai="$2"

    qc_bam_is_valid "$bam" || return 1
    python3 - "$bam" "$fai" <<'PY_QC_REFS'
import subprocess
import sys

bam, fai = sys.argv[1:]
expected = {}
with open(fai) as handle:
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        if len(fields) >= 2:
            name, length = fields[0], int(fields[1])
            if name in expected:
                raise SystemExit(f'ERROR: duplicate FASTA sequence ID: {name}')
            expected[name] = length

observed = {}
header = subprocess.check_output(['samtools', 'view', '-H', bam], text=True)
for line in header.splitlines():
    if not line.startswith('@SQ\t'):
        continue
    tags = {}
    for field in line.split('\t')[1:]:
        if ':' in field:
            key, value = field.split(':', 1)
            tags[key] = value
    if 'SN' in tags and 'LN' in tags:
        name = tags['SN']
        if name in observed:
            raise SystemExit(f'ERROR: duplicate BAM @SQ SN: {name}')
        observed[name] = int(tags['LN'])

if observed != expected:
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    wrong = sorted(
        name for name in set(expected) & set(observed)
        if expected[name] != observed[name]
    )
    print(
        f'BAM/reference mismatch: missing={len(missing)}, '
        f'unexpected={len(unexpected)}, wrong_length={len(wrong)}',
        file=sys.stderr,
    )
    raise SystemExit(1)
PY_QC_REFS
}

qc_compute_mapping_fingerprint() {
    local assembly="$1"
    local reads="$2"
    local preset="$3"
    local index_chunk="$4"
    local max_secondary="${MINIMAP2_MAX_SECONDARY:-5}"

    python3 - "$assembly" "$reads" "$preset" "$index_chunk" "$max_secondary" <<'PY_QC_FP'
import hashlib
import json
import os
import subprocess
import sys

assembly, reads, preset, index_chunk, max_secondary = sys.argv[1:]

def stat_payload(path):
    st = os.stat(path)
    return {
        'realpath': os.path.realpath(path),
        'size': st.st_size,
        'mtime_ns': st.st_mtime_ns,
    }

fai = assembly + '.fai'
digest = hashlib.sha256()
with open(fai, 'rb') as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        digest.update(chunk)

def command_first_line(command):
    raw = subprocess.check_output(command, stderr=subprocess.STDOUT)
    lines = raw.splitlines()
    if not lines:
        raise SystemExit(f'ERROR: command returned no version output: {command[0]}')
    return lines[0].decode('utf-8', errors='replace').strip()

version = command_first_line(['minimap2', '--version'])
samtools_version = command_first_line(['samtools', '--version'])

payload = {
    'assembly': stat_payload(assembly),
    'assembly_fai_sha256': digest.hexdigest(),
    'reads': stat_payload(reads),
    'preset': preset,
    'index_chunk': index_chunk,
    'mapping_arguments': [
        '-ax', preset, '-p', '1', '--secondary-seq', '-I', index_chunk,
        'samtools view -F 4'
    ],
    'minimap2_defaults_assumed': {
        'secondary': 'yes',
        'max_secondary': int(max_secondary),
    },
    'minimap2_version': version,
    'samtools_version': samtools_version,
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_QC_FP
}

qc_hash_text() {
    python3 - "$@" <<'PY_QC_TEXT_HASH'
import hashlib
import sys
print(hashlib.sha256('\0'.join(sys.argv[1:]).encode()).hexdigest())
PY_QC_TEXT_HASH
}

qc_compute_mapping_stage_fingerprints() {
    local assembly="$1"
    local reads="$2"
    local preset="$3"
    local index_chunk="$4"
    local mapping_fp samtools_version sorted_fp primary_fp

    mapping_fp="$(qc_compute_mapping_fingerprint \
        "$assembly" "$reads" "$preset" "$index_chunk")" || return 1
    samtools_version="$(samtools --version | awk 'NR==1{print}')" || return 1
    sorted_fp="$(qc_hash_text \
        "$mapping_fp" "sort_mapped_only" "$samtools_version")" || return 1
    primary_fp="$(qc_hash_text \
        "$sorted_fp" "samtools view -F 2308")" || return 1
    printf '%s\t%s\t%s\n' "$mapping_fp" "$sorted_fp" "$primary_fp"
}

qc_cache_matches() {
    local path="$1"
    local expected="$2"
    [[ -s "$path" ]] && [[ "$(<"$path")" == "$expected" ]]
}

qc_hifi_bam_mapping_semantically_matches() {
    local bam="$1"
    local assembly="$2"
    local reads="$3"
    local preset="$4"
    local index_chunk="$5"
    local max_secondary="${6:-${MINIMAP2_MAX_SECONDARY:-5}}"

    [[ -s "$bam" && -s "$assembly" && -s "$reads" && -s "${assembly}.fai" ]] || return 1
    samtools quickcheck "$bam" >/dev/null 2>&1 || return 1
    qc_bam_matches_fai "$bam" "${assembly}.fai" || return 1

    python3 - "$bam" "$assembly" "$reads" "$preset" "$index_chunk" "$max_secondary" <<'PY_QC_HIFI_SEMANTIC'
import os
import re
import shlex
import subprocess
import sys

bam, assembly, reads, preset, expected_chunk, expected_n = sys.argv[1:]
expected_n = int(expected_n)

try:
    header = subprocess.check_output(
        ["samtools", "view", "-H", bam], stderr=subprocess.STDOUT, text=True
    )
except Exception as exc:
    print(f"HiFi BAM semantic audit failed: cannot read header: {exc}", file=sys.stderr)
    raise SystemExit(1)

hd = next((line for line in header.splitlines() if line.startswith("@HD\t")), "")
if "\tSO:coordinate" not in hd:
    print("HiFi BAM semantic audit failed: BAM is not coordinate-sorted", file=sys.stderr)
    raise SystemExit(1)

pg_lines = [line for line in header.splitlines() if line.startswith("@PG\t")]
commands = []
for line in pg_lines:
    fields = line.split("\t")[1:]
    rec = {}
    for field in fields:
        if ":" in field:
            key, value = field.split(":", 1)
            rec[key] = value
    if rec.get("PN") == "minimap2" or "minimap2" in rec.get("CL", ""):
        commands.append(rec.get("CL", ""))

if not commands:
    print("HiFi BAM semantic audit failed: minimap2 @PG command not found", file=sys.stderr)
    raise SystemExit(1)

def option_value(tokens, opt):
    for i, tok in enumerate(tokens):
        if tok == opt and i + 1 < len(tokens):
            return tokens[i + 1]
        if tok.startswith(opt + "="):
            return tok.split("=", 1)[1]
        if len(opt) == 2 and tok.startswith(opt) and tok != opt:
            return tok[len(opt):]
    return None

def parse_size(value):
    m = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([kKmMgGtT]?)", value or "")
    if not m:
        return None
    number = float(m.group(1))
    suffix = m.group(2).upper()
    scale = {"": 1, "K": 1_000, "M": 1_000_000, "G": 1_000_000_000, "T": 1_000_000_000_000}[suffix]
    return int(number * scale)

reference_bases = 0
with open(assembly + ".fai", "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            reference_bases += int(parts[1])

expected_chunk_size = parse_size(expected_chunk)
assembly_real = os.path.realpath(assembly)
reads_real = os.path.realpath(reads)

failures = []
for command in commands:
    try:
        tokens = shlex.split(command)
    except ValueError:
        continue

    try:
        mm_i = next(i for i, tok in enumerate(tokens) if os.path.basename(tok) == "minimap2")
        tokens = tokens[mm_i + 1:]
    except StopIteration:
        pass

    preset_value = option_value(tokens, "-ax")
    if preset_value is None:
        x_value = option_value(tokens, "-x")
        has_a = "-a" in tokens
        preset_value = x_value if has_a else None
    p_value = option_value(tokens, "-p")
    n_value = option_value(tokens, "-N")
    chunk_value = option_value(tokens, "-I")

    if preset_value != preset:
        failures.append(f"preset={preset_value!r}")
        continue
    try:
        if float(p_value) != 1.0:
            failures.append(f"-p={p_value!r}")
            continue
    except (TypeError, ValueError):
        failures.append(f"-p={p_value!r}")
        continue
    try:
        if int(n_value) != expected_n:
            failures.append(f"-N={n_value!r}")
            continue
    except (TypeError, ValueError):
        failures.append(f"-N={n_value!r}")
        continue
    if "--secondary-seq" not in tokens:
        failures.append("--secondary-seq absent")
        continue
    if "--secondary=no" in tokens or "--secondary=no" in command:
        failures.append("secondary alignments disabled")
        continue

    observed_chunk_size = parse_size(chunk_value)
    chunk_ok = chunk_value == expected_chunk
    if not chunk_ok and observed_chunk_size is not None and expected_chunk_size is not None:
        chunk_ok = observed_chunk_size >= reference_bases and expected_chunk_size >= reference_bases
    if not chunk_ok:
        failures.append(f"-I={chunk_value!r}")
        continue

    existing_paths = []
    for tok in tokens:
        if tok == "-":
            continue
        try:
            if os.path.exists(tok):
                existing_paths.append(os.path.realpath(tok))
        except OSError:
            pass
    if assembly_real not in existing_paths:
        failures.append("assembly path absent from minimap2 command")
        continue
    if reads_real not in existing_paths:
        failures.append("read path absent from minimap2 command")
        continue

    raise SystemExit(0)

reason = failures[-1] if failures else "no parseable compatible minimap2 command"
print(f"HiFi BAM semantic audit failed: {reason}", file=sys.stderr)
raise SystemExit(1)
PY_QC_HIFI_SEMANTIC
}

qc_write_cache_atomic() {
    local path="$1"
    local value="$2"
    local tmp="${path}.tmp.$$"
    mkdir -p "$(dirname "$path")"
    printf '%s\n' "$value" > "$tmp"
    mv -f "$tmp" "$path"
}

qc_file_sha256() {
    local path="$1"
    python3 - "$path" <<'PY_QC_FILE_HASH'
import hashlib
import sys

digest = hashlib.sha256()
with open(sys.argv[1], 'rb') as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        digest.update(chunk)
print(digest.hexdigest())
PY_QC_FILE_HASH
}

qc_result_cache_matches() {
    local meta="$1"
    local input_fingerprint="$2"
    local result="$3"
    local cached_fingerprint cached_hash extra actual_hash

    [[ -s "$meta" && -s "$result" ]] || return 1
    IFS=$'\t' read -r cached_fingerprint cached_hash extra < "$meta"
    [[ -z "${extra:-}" ]] || return 1
    [[ "$cached_fingerprint" == "$input_fingerprint" ]] || return 1
    [[ "$cached_hash" =~ ^[0-9a-fA-F]{64}$ ]] || return 1
    actual_hash="$(qc_file_sha256 "$result")" || return 1
    [[ "${actual_hash,,}" == "${cached_hash,,}" ]]
}

qc_write_result_cache_atomic() {
    local meta="$1"
    local input_fingerprint="$2"
    local result="$3"
    local result_hash tmp

    result_hash="$(qc_file_sha256 "$result")"
    tmp="${meta}.tmp.$$"
    mkdir -p "$(dirname "$meta")"
    printf '%s\t%s\n' "$input_fingerprint" "$result_hash" > "$tmp"
    mv -f "$tmp" "$meta"
}

qc_compute_anvio_clipping_raw_fingerprint() {
    local clip_bam_meta="$1"
    local default_ratio="$2"
    local min_dist_to_end="$3"
    local algorithm_id="$4"
    local executable="$5"
    local anvio_version="$6"

    [[ -s "$clip_bam_meta" ]] || return 1
    [[ -x "$executable" ]] || return 1
    qc_hash_text \
        "$(<"$clip_bam_meta")" \
        "$algorithm_id" \
        "$default_ratio" \
        "$min_dist_to_end" \
        "$(readlink -f "$executable")" \
        "$anvio_version" \
        "--just-do-it"
}

qc_compute_anvio_clipping_result_fingerprint() {
    local raw_cache_meta="$1"
    local raw_clipping="$2"
    local stringent_depth="$3"
    local stringent_ratio="$4"
    local normalization_algorithm="$5"
    local validation_algorithm="$6"
    local normalizer_script="$7"
    local validator_script="$8"

    [[ -s "$raw_cache_meta" ]] || return 1
    [[ -s "$raw_clipping" ]] || return 1
    [[ -s "$normalizer_script" ]] || return 1
    [[ -s "$validator_script" ]] || return 1
    qc_hash_text \
        "$(<"$raw_cache_meta")" \
        "$(qc_file_sha256 "$raw_clipping")" \
        "$normalization_algorithm" \
        "$stringent_depth" \
        "$stringent_ratio" \
        "$(qc_file_sha256 "$normalizer_script")" \
        "$validation_algorithm" \
        "$(qc_file_sha256 "$validator_script")"
}

qc_compute_source_specific_clipping_fingerprint() {
    local manifest="$1"
    local default_ratio="$2"
    local min_dist_to_end="$3"
    local terminal_window="$4"
    local stringent_depth="$5"
    local stringent_ratio="$6"
    local discovery_algorithm="$7"
    local normalization_algorithm="$8"
    local validation_algorithm="$9"
    shift 9
    local normalizer_script="$1"
    local validator_script="$2"
    local classifier_script="$3"
    local executable="$4"
    local anvio_version="$5"

    python3 - \
        "$manifest" "$default_ratio" "$min_dist_to_end" "$terminal_window" \
        "$stringent_depth" "$stringent_ratio" \
        "$discovery_algorithm" "$normalization_algorithm" "$validation_algorithm" \
        "$normalizer_script" "$validator_script" "$classifier_script" \
        "$executable" "$anvio_version" <<'PY_QC_SOURCE_CLIP_FP'
import csv
import hashlib
import json
import os
import sys

(
    manifest, default_ratio, min_dist, terminal_window, stringent_depth,
    stringent_ratio, discovery_algorithm, normalization_algorithm,
    validation_algorithm, normalizer, validator, classifier, executable,
    anvio_version,
) = sys.argv[1:]

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def manifest_sha256(path):
    # The first seven manifest columns only: the flagstat counts appended by map do not change results.
    with open(path, 'rb') as handle:
        lines = handle.read().split(b'\n')
    core = b'\n'.join(b'\t'.join(line.split(b'\t')[:7]) for line in lines)
    return hashlib.sha256(core).hexdigest()

def stat_payload(path):
    st = os.stat(path)
    return {
        'realpath': os.path.realpath(path),
        'size': st.st_size,
        'mtime_ns': st.st_mtime_ns,
    }

for path in (manifest, normalizer, validator, classifier, executable):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f'ERROR: clipping fingerprint input missing: {path}')

sources = []
with open(manifest, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'Assembler', 'Assembly_FASTA', 'HiFi_all_BAM'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'ERROR: mapping manifest lacks columns: {sorted(missing)}')
    for row in reader:
        assembler = (row.get('Assembler') or '').strip()
        assembly = (row.get('Assembly_FASTA') or '').strip()
        bam = (row.get('HiFi_all_BAM') or '').strip()
        if not assembler or not assembly or not bam:
            raise SystemExit('ERROR: blank mapping manifest field')
        fai = assembly + '.fai'
        bam_meta = os.path.join(os.path.dirname(bam), '.cache', os.path.basename(bam) + '.inputs.sha256')
        for path in (assembly, fai, bam, bam_meta):
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise SystemExit(f'ERROR: clipping input missing: {path}')
        sources.append({
            'assembler': assembler,
            'assembly': stat_payload(assembly),
            'assembly_fai_sha256': sha256(fai),
            'hifi_bam': stat_payload(bam),
            'hifi_bam_mapping_fingerprint': open(bam_meta, encoding='utf-8').read().strip(),
        })

sources.sort(key=lambda x: x['assembler'])
payload = {
    'manifest_sha256': manifest_sha256(manifest),
    'sources': sources,
    'candidate_discovery': {
        'clipping_ratio': float(default_ratio),
        'min_dist_to_end': int(min_dist),
        'terminal_window_bp': int(terminal_window),
        'algorithm': discovery_algorithm,
        'anvio_executable': os.path.realpath(executable),
        'anvio_version': anvio_version,
    },
    'high_confidence': {
        'min_depth': int(stringent_depth),
        'clipping_ratio': float(stringent_ratio),
    },
    'normalization_algorithm': normalization_algorithm,
    'validation_algorithm': validation_algorithm,
    'normalizer_sha256': sha256(normalizer),
    'validator_sha256': sha256(validator),
    'classifier_sha256': sha256(classifier),
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_QC_SOURCE_CLIP_FP
}

qc_validate_current_clipping_result() {
    local anvio_version expected_final_fp

    [[ -s "$ASSEMBLER_MAPPING_MANIFEST" ]] || {
        echo "ERROR: assembler mapping manifest missing; run 01_mapping.sh" >&2
        return 1
    }
    [[ -s "$CLIP_FINAL" && -s "$CLIP_FINAL_CACHE_META" ]] || {
        echo "ERROR: clipping result/cache missing; run 03_clipping.sh" >&2
        return 1
    }
    [[ -x "$ANVIO_FIND_MISASSEMBLIES" ]] || return 1
    for path in "$CLIP_VALIDATE_SCRIPT" "$CLIP_SUMMARIZE_SCRIPT"; do
        [[ -s "$path" ]] || return 1
    done

    while IFS=$'\t' read -r assembler assembly hifi_bam _rest; do
        [[ "$assembler" == "Assembler" ]] && continue
        [[ -n "$assembler" && -s "$assembly" && -s "${assembly}.fai" && -s "$hifi_bam" ]] || {
            echo "ERROR: incomplete clipping input for $assembler" >&2
            return 1
        }
        samtools quickcheck -v "$hifi_bam" >/dev/null 2>&1 || {
            echo "ERROR: HiFi BAM failed quickcheck: $hifi_bam" >&2
            return 1
        }
        qc_bam_matches_fai "$hifi_bam" "${assembly}.fai" || {
            echo "ERROR: HiFi BAM/reference mismatch: $assembler" >&2
            return 1
        }
    done < "$ASSEMBLER_MAPPING_MANIFEST"

    anvio_version="$($ANVIO_ENV_PREFIX/bin/python - <<'PY_QC_ANVIO_VERSION'
import anvio
print(anvio.__version__)
PY_QC_ANVIO_VERSION
)" || return 1

    expected_final_fp="$(qc_compute_source_specific_clipping_fingerprint \
        "$ASSEMBLER_MAPPING_MANIFEST" "$ANVIO_DEFAULT_CLIPPING_RATIO" \
        "$MIN_DIST_TO_END" "$TERMINAL_CLIPPING_WINDOW_BP" \
        "$STRINGENT_CLIP_MIN_DEPTH" "$STRINGENT_CLIP_RATIO" \
        "$ANVIO_CLIPPING_ALGORITHM_ID" "$CLIPPING_NORMALIZATION_ALGORITHM_ID" \
        "$CLIPPING_VALIDATION_ALGORITHM_ID" "$CLIP_VALIDATE_SCRIPT" \
        "$CLIP_VALIDATE_SCRIPT" "$CLIP_SUMMARIZE_SCRIPT" \
        "$ANVIO_FIND_MISASSEMBLIES" "$anvio_version")" || return 1

    qc_result_cache_matches \
        "$CLIP_FINAL_CACHE_META" "$expected_final_fp" "$CLIP_FINAL" || {
        echo "ERROR: final clipping cache is stale; run 03_clipping.sh" >&2
        return 1
    }

    python3 - "$CLIP_FINAL" <<'PY_QC_INTERNAL_CLIP'
import csv
import sys
path = sys.argv[1]
with open(path, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'Assembler', 'region_type', 'Clipping'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'ERROR: clipping table lacks new schema columns: {sorted(missing)}')
    bad = [row for row in reader if (row.get('region_type') or '').strip() != 'INTERNAL']
if bad:
    raise SystemExit(f'ERROR: internal clipping table must contain INTERNAL rows only: {path}')
PY_QC_INTERNAL_CLIP
}

qc_parallel_jobs() {
    local items="$1" jobs="${PARALLEL_JOBS:-$1}"
    (( jobs > items )) && jobs="$items"
    (( jobs < 1 )) && jobs=1
    echo "$jobs"
}

qc_duration() {
    local s="$1"
    if (( s >= 3600 )); then printf '%dh%02dm' $(( s / 3600 )) $(( s % 3600 / 60 ))
    elif (( s >= 60 )); then printf '%dm' $(( s / 60 ))
    else printf '%ds' "$s"; fi
}

qc_progress_watch() {
    local label="$1" total="$2" dir="$3" pattern="$4" parent="$$"
    local interval="${PROGRESS_INTERVAL:-300}" start="$SECONDS" finished elapsed remaining
    (( total > 0 )) || return 0
    while sleep "$interval"; do
        kill -0 "$parent" 2>/dev/null || break
        finished="$(find "$dir" -maxdepth 1 -name "$pattern" 2>/dev/null | wc -l)"
        elapsed=$(( SECONDS - start ))
        remaining="estimating"
        (( finished > 0 )) && remaining="$(qc_duration $(( elapsed * (total - finished) / finished )))"
        echo "[$label] $finished/$total ($(( 100 * finished / total ))%) | elapsed $(qc_duration "$elapsed") | remaining $remaining"
    done
}

qc_progress_stop() {
    kill "$1" 2>/dev/null || true
    wait "$1" 2>/dev/null || true
}

qc_split_threads() {
    local threads="$1" jobs="$2" each
    each=$(( threads / jobs ))
    (( each < 1 )) && each=1
    echo "$each"
}

qc_wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do
        sleep 1
    done
}

qc_wait_all() {
    local pid
    for pid in $(jobs -rp); do
        wait "$pid" 2>/dev/null || true
    done
}

qc_manifest_columns() {
    python3 - "$@" <<'PY_MANIFEST_COLUMNS'
import csv, sys
path, *columns = sys.argv[1:]
with open(path, encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    missing = [c for c in columns if c not in (reader.fieldnames or [])]
    if missing:
        sys.exit(f"ERROR: mapping manifest lacks columns {missing}: {path}")
    for row in reader:
        if (row.get("Assembler") or "").strip():
            print("\t".join((row.get(c) or "").strip() for c in columns))
PY_MANIFEST_COLUMNS
}

qc_bam_cache() {
    printf '%s/.cache/%s.inputs.sha256\n' "$(dirname "$1")" "$(basename "$1")"
}
