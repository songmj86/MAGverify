#!/usr/bin/env bash
set -euo pipefail
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QC_STAGE="${QC_STAGE:-$(basename "${BASH_SOURCE[1]:-lib_config.sh}" .sh)}"

# Each analysis stage accepts optional positional input and output paths.
# Additional inputs arrive through the CLI environment or optional shell variables.
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
                01_mapping) description='Map PacBio HiFi reads to each original source assembly'; input='HiFi FASTQ or FASTQ.GZ'; extras='ASSEMBLER_ASSEMBLY_DIR, ASSEMBLER_FASTAS'; output='OUTPUT_DIR/mapping_manifest.tsv' ;;
            03_clipping) description='Identify clipping and validate read support'; input='Combined assembly FASTA'; extras='MAG_DIR, ASSEMBLER_MAPPING_MANIFEST with source-specific HiFi BAMs from 01_mapping; MAG map generated automatically'; output='OUTPUT_DIR/clipping.all_regions.filtered.tsv' ;;
            07_chimera) description='Assess genome-wide GTDB taxonomic transitions'; input='GTDB R232 gene-taxonomy TSV'; extras='ASSEMBLY_FA, MAG_CONTIG_MAP, mapping manifest, 03_clipping result'; output='OUTPUT_DIR/GTDB_R232/clipping_defined/' ;;
            02_circularity) description='Validate end-to-start closure of single-contig MAGs'; input='Directory of MAG FASTA files'; extras='READS, ASSEMBLER_ASSEMBLY_DIR, ASSEMBLER_FASTAS'; output='OUTPUT_DIR/MAG_circularity.final.tsv' ;;
            04_false_closure) description='Evaluate false-closure candidates'; input='Combined assembly FASTA'; extras='READS, MAG_CONTIG_MAP, mapping, clipping, circularity results'; output='OUTPUT_DIR/MAG_false_closure.final.tsv' ;;
            09_assembly_variants) description='Quantify assembly allele support using anvi-o'; input='Combined assembly FASTA'; extras='READS, MAG_CONTIG_MAP, mapping manifest, active MAGverify Conda environment'; output='OUTPUT_DIR/contig_assembly_allele_support.final.tsv' ;;
            05_phantom) description='Assess zero-coverage sequences and read k-mer evidence'; input='Combined assembly FASTA'; extras='READS must be four-line FASTQ(.GZ), MAG_CONTIG_MAP, mapping manifest'; output='OUTPUT_DIR/phantom_sequences.final.tsv' ;;
            06_repeats) description='Quantify intra- and inter-contig repeats'; input='Combined assembly FASTA'; extras='MAG_CONTIG_MAP, source assembler FASTAs'; output='OUTPUT_DIR/contig_repeat_breadth.final.tsv' ;;
            08_false_duplication) description='Assess candidate false duplications'; input='Combined assembly FASTA'; extras='MAG_CONTIG_MAP, mapping manifest, clipping result'; output='OUTPUT_DIR/false_duplication.by_contig.tsv' ;;
            10_summary) description='Integrate and validate enabled QC results'; input='Combined assembly FASTA'; extras='MAG_CONTIG_MAP and outputs from all enabled analyses'; output='OUTPUT_DIR/integrated_QC_summary.by_MAG.tsv' ;;
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

# Project paths and optional software environments.
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

MAX_THREADS="$(detect_max_threads)"
THREADS="${THREADS:-$MAX_THREADS}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: THREADS must be a positive integer: $THREADS" >&2
    return 1 2>/dev/null || exit 1
}
if (( THREADS > MAX_THREADS )); then
    echo "ERROR: THREADS=$THREADS exceeds the detected allocation of $MAX_THREADS" >&2
    return 1 2>/dev/null || exit 1
fi
PRESET="map-hifi"

# Stage enablement for QC summaries and preflight checks.
RUN_CLIPPING="${RUN_CLIPPING:-Y}"
RUN_TAXONOMIC_CHIMERA="${RUN_TAXONOMIC_CHIMERA:-Y}"
RUN_CIRCULARITY="${RUN_CIRCULARITY:-Y}"
RUN_FALSE_CLOSURE="${RUN_FALSE_CLOSURE:-Y}"
RUN_ASSEMBLY_VARIANTS="${RUN_ASSEMBLY_VARIANTS:-Y}"
RUN_PHANTOM="${RUN_PHANTOM:-Y}"
RUN_REPEAT="${RUN_REPEAT:-Y}"
RUN_FALSE_DUPLICATION="${RUN_FALSE_DUPLICATION:-Y}"

# Full assembly FASTA generated from the same long reads.
ASSEMBLY_FA="${ASSEMBLY_FA:-}"

# Raw PacBio HiFi reads used for assembly.
READS="${READS:-}"


# One FASTA per MAG. Circularity validation is restricted to single-contig MAGs.
MAG_DIR="${MAG_DIR:-}"
MAG_EXTENSIONS="fa fasta fna"

# Original source-assembler co-assemblies. 01_mapping.sh and the canonical
# rotated circularity workflow use the same three references. In automatic
# circularity mode, each eligible single-contig MAG is assigned to exactly one
# source assembly by contig ID plus exact sequence identity.
ASSEMBLER_ASSEMBLY_DIR="${ASSEMBLER_ASSEMBLY_DIR:-}"
ASSEMBLER_FASTAS="${ASSEMBLER_FASTAS:-hifiasm_meta.fa metaflye.fa metaMDBG.fa}"

# Output directory.
OUT="${OUT:-$CONFIG_DIR/results}"

# Load all output directories, then apply only the active module's input overrides.
source "$CONFIG_DIR/stage_paths.sh"
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

# Positional arguments take precedence over per-module input overrides.
if [[ -n "$QC_CLI_INPUT" ]]; then
    case "$QC_STAGE" in
        01_mapping) READS="$QC_CLI_INPUT" ;;
        07_chimera) GTDB_GENE_TAX_TSV="$QC_CLI_INPUT" ;;
        02_circularity) MAG_DIR="$QC_CLI_INPUT" ;;
        *) ASSEMBLY_FA="$QC_CLI_INPUT" ;;
    esac
fi

# Optional but strongly recommended all-MAG contig ownership table.
# Format: MAG_ID<TAB>Contig. In addition to aggregation, this identifies which
# other MAG receives reads continuing from a premature-circularization site.
# When absent, MAGverify prepares a target-MAG map; that fallback cannot
# name non-target destination MAGs.
PREP_DIR="${PREP_DIR:-$OUT/00_prep}"
MAG_CONTIG_MAP="${MAG_CONTIG_MAP:-$PREP_DIR/MAG_contig_map.tsv}"

# Assembler-reported circularity table used by 04_false_closure.sh. The canonical
# rotated read-support module does not require or gate on assembler circularity.
# Recognized ID columns: ID, MAG_ID, genome, genome_id, bin, bin_id
# Recognized call columns: assembler_circularity, circularity, circular,
#                          is_circular.
CIRCULARITY_TSV="${CIRCULARITY_TSV:-}"


# Raw anvi'o HMM-hit table used by the Trigodet et al. premature-
# circularization proxy. lib_scg.sh writes the Bacteria_71 and Archaea_76
# ribosomal-protein hits with original contig names into this directory.
SCG_OUT_DIR="${SCG_OUT_DIR:-$PREP_DIR/SCG}"
SCG_TSV="${SCG_TSV:-$SCG_OUT_DIR/SCG_from_anvio.original_contig.tsv}"

# Shared mapping and clipping parameters

# Candidate discovery is performed by the actual anvi'o program. The documented
# default clipping-ratio threshold (0.80) is retained. The default 100-bp
# terminal exclusion is intentionally disabled so terminal clipping can be
# measured directly and compared with internal clipping per analyzed base.
ANVIO_DEFAULT_CLIPPING_RATIO=0.80
MIN_DIST_TO_END=0
TERMINAL_CLIPPING_WINDOW_BP=100

# A stringent candidate is a subset of the anvi'o default candidates. This is
# the exact-site depth/fraction tier used for the Trigodet et al. comparison.
STRINGENT_CLIP_MIN_DEPTH=10
STRINGENT_CLIP_RATIO=1.0

# Read-level validation is performed for every default candidate. Clipping
# candidate discovery and validation use ALL mapped alignment records: primary,
# secondary and supplementary. Validation is deduplicated by QNAME, so each
# HiFi read is counted once per candidate regardless of how many alignment
# records it has. No MAPQ, depth, QNAME-fraction, or fixed spanning-flank
# threshold is imposed on final read-level validation.

ANVIO_CLIPPING_ALGORITHM_ID="anvi_script_find_misassemblies_source_assembler_all_mapped_end0_v7"
CLIPPING_NORMALIZATION_ALGORITHM_ID="anvio_exact_default_and_stringent_source_specific_v2"
CLIPPING_VALIDATION_ALGORITHM_ID="qname_source_assembler_all_mapped_boundary_safe_exact_idn_length_audit_v9"

# Robust mapping / clipping runtime parameters

# samtools sort/view threads and memory per sort thread. Keep conservative when
# THREADS is high to avoid disk-I/O and temporary-file pressure.
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-$THREADS}"
# anvi'o discovery and read-level validation use the full CPU allocation visible
# to this job unless explicitly overridden.
CLIPPING_THREADS="${CLIPPING_THREADS:-$THREADS}"
SORT_MEM="2G"

# Robust minimap2 mode.
# Mapping uses a single minimap2 index part by default.
FORCE_SINGLE_PART_MINIMAP2="Y"
MINIMAP2_INDEX_CHUNK_SINGLE="900G"
MINIMAP2_INDEX_CHUNK="256G"
# Retain at most five secondary alignments per primary alignment. 01_mapping.sh
# and 02_circularity.sh both use -p 1 -N 5 --secondary-seq. Rotated circularity
# uses primary+supplementary as the main structural evidence; all-mapped records
# are retained only as an audit of secondary-alignment contributions.
MINIMAP2_MAX_SECONDARY=5


# Prefer CSI index for very large references/contigs.
INDEX_CSI="Y"

# Rebuild the full-assembly HiFi mapping even when the cache fingerprint
# matches. The environment-variable form supports a one-off forced rerun:
# FORCE_ASSEMBLY_REMAP=Y bash 01_mapping.sh
FORCE_ASSEMBLY_REMAP="${FORCE_ASSEMBLY_REMAP:-N}"


# anvi'o is used for clipping-candidate discovery and for modules that build
# contigs databases or profiles. Read-level clipping validation is performed
# directly from BAM alignments.
ANVIO_ENV_PREFIX="${ANVIO_ENV_PREFIX:-${CONDA_PREFIX:-}}"
if [[ -z "$ANVIO_ENV_PREFIX" ]] && command -v anvi-profile >/dev/null 2>&1; then
    ANVIO_ENV_PREFIX="$(cd "$(dirname "$(command -v anvi-profile)")/.." && pwd -P)"
fi
ANVIO_THREADS="${ANVIO_THREADS:-$THREADS}"

# Recompute clipping discovery and read-level validation even when
# the result fingerprint matches.
FORCE_RERUN_CLIPPING="${FORCE_RERUN_CLIPPING:-N}"
CLIPPING_POSTPROCESS_ONLY="${CLIPPING_POSTPROCESS_ONLY:-N}"
CLIPPING_MAG_SUMMARY_ONLY="${CLIPPING_MAG_SUMMARY_ONLY:-N}"


# Rotated END|START circularity validation

# Canonical circularity implementation. The reference is the corresponding original
# co-assembly after rotating each eligible single-contig MAG so
# that its original END|START adjacency lies at an internal coordinate. Rotation
# itself is not evidence; direct read spanning across that original adjacency is.
CIRCULARITY_SCRIPT="$CONFIG_DIR/02_circularity.sh"
CIRCULARITY_ROTATED_RESULT_DIR="${CIRCULARITY_ROTATED_RESULT_DIR:-$CIRCULARITY_DIR}"

# HiFi main evidence: primary + supplementary mapped alignments, no MAPQ cutoff,
# actual M/=X reference alignment >=2,000 bp on each side, and >=2 unique QNAMEs.
# Insertions at END|START are allowed. A CIGAR D or N crossing END|START invalidates
# direct topology support. All mapped records are retained only as a secondary-
# alignment audit and do not change the final circularity call.
CIRCULARITY_HIFI_ANCHOR_BP=2000
CIRCULARITY_MIN_HIFI_QNAMES=2

# Runtime values mirror 02_circularity.sh defaults.
CIRCULARITY_THREADS="${CIRCULARITY_THREADS:-$(( THREADS < 100 ? THREADS : 100 ))}"
CIRCULARITY_SORT_THREADS="${CIRCULARITY_SORT_THREADS:-$(( THREADS < 20 ? THREADS : 20 ))}"
CIRCULARITY_SORT_MEMORY="2G"
CIRCULARITY_MINIMAP2_INDEX_CHUNK="$MINIMAP2_INDEX_CHUNK_SINGLE"
CIRCULARITY_MINIMAP2_MAX_SECONDARY="$MINIMAP2_MAX_SECONDARY"

# Trigodet et al. (Nature Biotechnology, 2026) premature-circularization proxy:
# an assembler-reported circular contig strictly shorter than 500 kbp with at
# least three distinct ribosomal-protein HMM families (Bacteria_71/Archaea_76).
# This is a screening proxy only and never overrides the HiFi false-closure call.
PREMATURE_MAX_LENGTH_BP=500000
PREMATURE_MIN_RIBOSOMAL_FAMILIES=3
# False closure remains an independent analysis. At each candidate, at least two
# independent clipped QNAMEs on each side must show a primary/supplementary
# cross-contig connection; the two sides may connect to different contigs.
FALSE_CLOSURE_MIN_LINKED_READS_PER_SIDE=2

# Assembly-allele support analysis.

# anvi-profile sample label. Keep to letters, numbers, dot, underscore or dash.
ASSEMBLY_VARIANT_SAMPLE_ID="${ASSEMBLY_VARIANT_SAMPLE_ID:-HiFi_reads}"
# Rebuild the anvi'o contigs/profile databases and all exported summaries.
FORCE_RERUN_ASSEMBLY_VARIANTS="${FORCE_RERUN_ASSEMBLY_VARIANTS:-N}"

# HiFi anvi'o profiling is performed once per complete assembler-specific
# co-assembly. Assemblers are processed serially with no contig batching. The
# profile thread count is derived from a memory budget in 09_assembly_variants.sh.
# Profile threads are bounded by the available memory budget.
ASSEMBLY_VARIANT_MEMORY_BUDGET_GB="${ASSEMBLY_VARIANT_MEMORY_BUDGET_GB:-900}"
ASSEMBLY_VARIANT_MEMORY_RESERVE_GB="${ASSEMBLY_VARIANT_MEMORY_RESERVE_GB:-64}"
ALLELE_ANVIO_MEMORY_PER_THREAD_GB="${ALLELE_ANVIO_MEMORY_PER_THREAD_GB:-16}"
ALLELE_ANVIO_MEMORY_POLL_SECONDS="${ALLELE_ANVIO_MEMORY_POLL_SECONDS:-5}"
ALLELE_ANVIO_FULL_PROFILE_THREADS="${ALLELE_ANVIO_FULL_PROFILE_THREADS:-0}"  # 0 = memory-derived auto
ALLELE_ANVIO_PREPROCESS_THREADS="${ALLELE_ANVIO_PREPROCESS_THREADS:-70}"
# The exported SNV/INDEL tables and split map are the persistent analysis inputs.
# PROFILE.db is removed after successful export unless explicitly retained.
ASSEMBLY_VARIANT_KEEP_FULL_PROFILE="${ASSEMBLY_VARIANT_KEEP_FULL_PROFILE:-N}"
ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES="${ASSEMBLY_VARIANT_REMOVE_LEGACY_BATCHES:-Y}"
# Best-effort POSIX_FADV_DONTNEED is applied to large per-assembler files after
# use so Linux can release their clean page-cache pages without root privileges.
ASSEMBLY_VARIANT_DROP_FILE_CACHE="${ASSEMBLY_VARIANT_DROP_FILE_CACHE:-Y}"
# Keep gene calling by default so gene-level columns/counts remain identical to
# the previous workflow. Set Y only if gene-level variant summaries are dispensable.
ASSEMBLY_VARIANT_SKIP_GENE_CALLING="${ASSEMBLY_VARIANT_SKIP_GENE_CALLING:-N}"

# HiFi classification follows the Trigodet et al. reproducible workflow.
# SNV: unsupported when assembly support is zero; otherwise minority when the
# dominant non-assembly nucleotide exceeds 50% of nucleotide coverage.
# INDEL: event count == anvi'o event coverage -> unsupported; otherwise
# event count / coverage >0.5 -> minority. Exact ties are not called.

# Phantom sequence parameters. Independent module.

# Paper Methods searched exact zero-coverage regions longer than 500 bp against
# the original long reads with BLASTN -dust no.
ZERO_MIN_BP=501
PHANTOM_BLAST_THREADS="$THREADS"
PHANTOM_KMER_K=21
FORCE_RERUN_PHANTOM="N"

# Excessive repeat parameters. Independent module.

# Repeat alignments: >200 bp and >=80% nucleotide identity. Intra-contig
# repeats are assessed for all contigs; inter-contig repeats are additionally
# assessed between different contigs within each multi-contig MAG.
REPEAT_MIN_ALIGN_PAPER=200
REPEAT_PIDENT_PAPER=80
# Both thresholds are 70%, but are kept separate because they define different
# calls: contig-level excessive repeat and genome-wide excessive repeat.
REPEAT_CONTIG_HIGH=0.70
REPEAT_GENOME_HIGH=0.70
# Legacy alias retained for older scripts until 10_summary.sh is revised.
REPEAT_HIGH="$REPEAT_CONTIG_HIGH"
REPEAT_CHUNK_BP=10000000
REPEAT_BLAST_THREADS_PER_JOB=2
REPEAT_JOBS="${REPEAT_JOBS:-$(( THREADS / REPEAT_BLAST_THREADS_PER_JOB ))}"
(( REPEAT_JOBS >= 1 )) || REPEAT_JOBS=1
REPEAT_SORT_THREADS="$THREADS"
REPEAT_SORT_MEMORY="50%"
REPEAT_MIN_CONTIG_LEN=0
FORCE_RERUN_REPEAT="N"


# False-duplication screening

# 08_false_duplication.sh is independent of 06_repeats.sh. Candidate regions
# are obtained by Minimap2 assembly self-alignment and evaluated with
# source-assembler HiFi boundaries. Intra-contig tandem/interspersed and
# direct/inverted arrangements are recorded; inter-contig structural analysis
# is restricted to contig pairs assigned to the same MAG.
FORCE_RERUN_FALSE_DUPLICATION="N"
FALSE_DUP_MAPPING_MANIFEST="${FALSE_DUP_MAPPING_MANIFEST:-${ASSEMBLER_MAPPING_MANIFEST:-$MAPPING_DIR/mapping_manifest.tsv}}"
FALSE_DUP_MIN_COPY_BP="${FALSE_DUP_MIN_COPY_BP:-1000}"
FALSE_DUP_MIN_PIDENT="${FALSE_DUP_MIN_PIDENT:-95}"
# Final false-duplication boundary anchor is fixed at 500 bp.
FALSE_DUP_HIFI_ANCHOR_BP=500
FALSE_DUP_HIFI_MIN_SPANNING_READS="${FALSE_DUP_HIFI_MIN_SPANNING_READS:-2}"
FALSE_DUP_HIFI_MIN_REACHING_READS="${FALSE_DUP_HIFI_MIN_REACHING_READS:-2}"


# Derived paths

# Read evidence is routed through mapping_manifest.tsv to the complete
# co-assembly produced by the assembler that generated each evaluated contig.
# The concatenated ASSEMBLY_FA remains available for catalogue/sequence analyses
# that need the union of contigs, but it is not a read-mapping reference.

# Clipping candidate discovery and read-level validation use source-assembler
# HiFi_all_BAM files from mapping_manifest.tsv, retaining all mapped
# primary, secondary and supplementary alignment records.

if [[ -z "${ANVIO_FIND_MISASSEMBLIES:-}" ]]; then
    if [[ -n "$ANVIO_ENV_PREFIX" && -x "$ANVIO_ENV_PREFIX/bin/anvi-script-find-misassemblies" ]]; then
        ANVIO_FIND_MISASSEMBLIES="$ANVIO_ENV_PREFIX/bin/anvi-script-find-misassemblies"
    else
        ANVIO_FIND_MISASSEMBLIES="$(command -v anvi-script-find-misassemblies || true)"
    fi
fi
ANVIO_CLIP_PREFIX="${ANVIO_CLIP_PREFIX:-${CLIPPING_DIR}/anvio_find_misassemblies}"
ANVIO_CLIP_RAW="${ANVIO_CLIP_RAW:-${ANVIO_CLIP_PREFIX}-clipping.txt}"
ANVIO_ZERO_RAW="${ANVIO_ZERO_RAW:-${ANVIO_CLIP_PREFIX}-zero_cov.txt}"
CLIP_DEFAULT="${CLIP_DEFAULT:-${CLIPPING_DIR}/clipping.default.tsv}"
CLIP_STRINGENT="${CLIP_STRINGENT:-${CLIPPING_DIR}/clipping.stringent.tsv}"
CLIP_FINAL="${CLIP_FINAL:-${CLIPPING_DIR}/clipping.filtered.tsv}"
CLIP_POSITIVE="${CLIP_POSITIVE:-${CLIPPING_DIR}/clipping.Y.tsv}"
CLIP_MAG_SITE="${CLIP_MAG_SITE:-${CLIPPING_DIR}/MAG_clipping_sites.final.tsv}"
CLIP_MAG_POSITIVE_SITE="${CLIP_MAG_POSITIVE_SITE:-${CLIPPING_DIR}/MAG_clipping_sites.Y.tsv}"
CLIP_MAG_SUMMARY="${CLIP_MAG_SUMMARY:-${CLIPPING_DIR}/MAG_clipping.final.tsv}"
CLIP_MAG_POSITIVE="${CLIP_MAG_POSITIVE:-${CLIPPING_DIR}/MAG_clipping.Y.tsv}"
CLIP_CONTIG_SUMMARY="${CLIP_CONTIG_SUMMARY:-${CLIPPING_DIR}/contig_clipping.final.tsv}"
CLIP_FINAL_CACHE_META="${CLIP_FINAL_CACHE_META:-${CLIP_FINAL}.inputs.sha256}"
CLIP_CRITERIA="${CLIP_CRITERIA:-${CLIPPING_DIR}/clipping.criteria.tsv}"

CLIP_NORMALIZER_SCRIPT="${CLIP_NORMALIZER_SCRIPT:-$CONFIG_DIR/scripts/clipping_normalize_anvio_clipping.py}"
CLIP_SPANNING_VALIDATOR="${CLIP_SPANNING_VALIDATOR:-$CONFIG_DIR/scripts/clipping_validate_clipping_spanning.py}"
CLIP_MAG_SUMMARIZER_SCRIPT="${CLIP_MAG_SUMMARIZER_SCRIPT:-$CONFIG_DIR/scripts/clipping_summarize_clipping_by_mag.py}"
CLIP_REGION_CLASSIFIER_SCRIPT="${CLIP_REGION_CLASSIFIER_SCRIPT:-$CONFIG_DIR/scripts/clipping_classify_clipping_regions.py}"
CLIP_REGION_SUMMARIZER_SCRIPT="${CLIP_REGION_SUMMARIZER_SCRIPT:-$CONFIG_DIR/scripts/clipping_summarize_clipping_regions.py}"
MAG_MAP_BUILDER="${MAG_MAP_BUILDER:-$CONFIG_DIR/scripts/prep_build_mag_contig_map.py}"
ASSEMBLER_MAPPING_MANIFEST="${ASSEMBLER_MAPPING_MANIFEST:-${MAPPING_DIR}/mapping_manifest.tsv}"
VALIDATION_VIEW_EXPORTER="${VALIDATION_VIEW_EXPORTER:-$CONFIG_DIR/scripts/common_export_validation_view.py}"

# Gene-taxonomy generation remains a separate utility. The chimera module itself
# uses GTDB R232 only.

# GTDB R232 taxonomic-chimera inputs and outputs.
GTDB_RELEASE="${GTDB_RELEASE:-GTDB_R232}"
GTDB_GENE_TAX_DIR="${GTDB_GENE_TAX_DIR:-$PREP_DIR/gene_taxonomy/GTDB_R232}"
GTDB_GENE_TAX_TSV="${GTDB_GENE_TAX_TSV:-$GTDB_GENE_TAX_DIR/gene_taxonomy.mmseqs2.GTDB_R232.all_contigs.tsv}"
GTDB_CHIMERA_DIR="${GTDB_CHIMERA_DIR:-${CHIMERA_ROOT_DIR}/$GTDB_RELEASE}"
GTDB_CHIMERA_CONTIG="${GTDB_CHIMERA_CONTIG:-$GTDB_CHIMERA_DIR/contig_structural_chimera.${GTDB_RELEASE}.final.tsv}"
GTDB_CHIMERA_MAG_CONTIG="${GTDB_CHIMERA_MAG_CONTIG:-$GTDB_CHIMERA_DIR/MAG_contig_structural_chimera.${GTDB_RELEASE}.final.tsv}"
GTDB_CHIMERA_MAG="${GTDB_CHIMERA_MAG:-$GTDB_CHIMERA_DIR/MAG_structural_chimera.${GTDB_RELEASE}.final.tsv}"
MAG_GTDB_TAXONOMY_TSV="${MAG_GTDB_TAXONOMY_TSV:-}"

# Circularity and false closure are separate pipeline modules.
CIRCULARITY_RESULT_DIR="${CIRCULARITY_RESULT_DIR:-$CIRCULARITY_DIR}"
FALSE_CLOSURE_RESULT_DIR="${FALSE_CLOSURE_RESULT_DIR:-$FALSE_CLOSURE_DIR}"

CIRCULARITY_ONLY_FINAL="${CIRCULARITY_ONLY_FINAL:-$CIRCULARITY_RESULT_DIR/MAG_circularity.final.tsv}"
CIRCULARITY_ONLY_AUDIT="${CIRCULARITY_ONLY_AUDIT:-$CIRCULARITY_RESULT_DIR/MAG_circularity.audit.tsv}"
CIRCULARITY_RUN_METADATA="${CIRCULARITY_RUN_METADATA:-$CIRCULARITY_RESULT_DIR/run_metadata.txt}"

FALSE_CLOSURE_FINAL="${FALSE_CLOSURE_FINAL:-$FALSE_CLOSURE_RESULT_DIR/MAG_false_closure.final.tsv}"
FALSE_CLOSURE_AUDIT="${FALSE_CLOSURE_AUDIT:-$FALSE_CLOSURE_RESULT_DIR/MAG_false_closure.audit.tsv}"
FALSE_CLOSURE_CONTIG_FINAL="${FALSE_CLOSURE_CONTIG_FINAL:-$FALSE_CLOSURE_RESULT_DIR/contig_false_closure.final.tsv}"
FALSE_CLOSURE_RUN_METADATA="${FALSE_CLOSURE_RUN_METADATA:-$FALSE_CLOSURE_RESULT_DIR/run_metadata.txt}"
FALSE_CLOSURE_CANDIDATE_READS="${FALSE_CLOSURE_CANDIDATE_READS:-$FALSE_CLOSURE_RESULT_DIR/false_closure_candidate_reads.tsv}"
FALSE_CLOSURE_EVIDENCE_TSV="${FALSE_CLOSURE_EVIDENCE_TSV:-$FALSE_CLOSURE_RESULT_DIR/false_closure.tsv}"
FALSE_CLOSURE_BREAKPOINTS_TSV="${FALSE_CLOSURE_BREAKPOINTS_TSV:-$FALSE_CLOSURE_RESULT_DIR/false_closure.breakpoints.tsv}"

# Assembly-variant profiling uses the all-mapped HiFi population. Keep this
# cache separate from earlier primary+supplementary or secondary-excluded runs.
ALLELE_ANVIO_DIR="${ALLELE_ANVIO_DIR:-${VARIANT_DIR}/hifi_all_mapped/anvio}"
# Assembly variants are profiled separately across every contig of each
# assembler-specific complete co-assembly. Final-MAG outputs are generated by
# subsetting the resulting assembly-wide classifications, not by re-profiling.
ALLELE_BY_ASSEMBLER_DIR="${ALLELE_BY_ASSEMBLER_DIR:-$ALLELE_ANVIO_DIR/by_assembler}"
ALLELE_PROFILE_MANIFEST="${ALLELE_PROFILE_MANIFEST:-${VARIANT_DIR}/anvio_profiles.by_assembler.tsv}"
ALLELE_VARIANTS="${ALLELE_VARIANTS:-${VARIANT_DIR}/assembly_allele_disagreements.tsv}"
ALLELE_CONTIG_ALL="${ALLELE_CONTIG_ALL:-${VARIANT_DIR}/contig_assembly_allele_support.all_contigs.tsv}"
ALLELE_SUMMARY_ALL="${ALLELE_SUMMARY_ALL:-${VARIANT_DIR}/assembly_allele_support.all_contigs.summary.tsv}"
ALLELE_ASSEMBLER_SUMMARY="${ALLELE_ASSEMBLER_SUMMARY:-${VARIANT_DIR}/assembly_allele_support.by_assembler.tsv}"
# Existing MAG-level names are retained for downstream summaries.
ALLELE_CONTIG="${ALLELE_CONTIG:-${VARIANT_DIR}/contig_assembly_allele_support.final.tsv}"
ALLELE_SUMMARY="${ALLELE_SUMMARY:-${VARIANT_DIR}/assembly_allele_support.summary.tsv}"
ALLELE_VARIANTS_MAG="${ALLELE_VARIANTS_MAG:-${VARIANT_DIR}/assembly_allele_disagreements.by_MAG.tsv}"
ALLELE_MAG="${ALLELE_MAG:-${VARIANT_DIR}/assembly_allele_support.by_MAG.tsv}"
ALLELE_CACHE_META="${ALLELE_CACHE_META:-${VARIANT_DIR}/assembly_allele_support.inputs.sha256}"

PHANTOM_ZERO_BED="${PHANTOM_ZERO_BED:-${PHANTOM_DIR}/zero_regions.bed}"
PHANTOM_ZERO_FA="${PHANTOM_ZERO_FA:-${PHANTOM_DIR}/zero_regions.fa}"
PHANTOM_ZERO_TSV="${PHANTOM_ZERO_TSV:-${PHANTOM_DIR}/zero_regions.tsv}"
PHANTOM_BLAST_DB="${PHANTOM_BLAST_DB:-${PHANTOM_DIR}/raw_reads_blastdb/reads}"
PHANTOM_BLAST="${PHANTOM_BLAST:-${PHANTOM_DIR}/zero_regions_vs_raw_reads.tsv}"
PHANTOM_FINAL="${PHANTOM_FINAL:-${PHANTOM_DIR}/phantom_sequences.final.tsv}"
PHANTOM_CACHE_META="${PHANTOM_CACHE_META:-${PHANTOM_FINAL}.inputs.sha256}"
PHANTOM_MAG_SITES="${PHANTOM_MAG_SITES:-${PHANTOM_DIR}/MAG_phantom_sequences.final.tsv}"
PHANTOM_MAG_SUMMARY="${PHANTOM_MAG_SUMMARY:-${PHANTOM_DIR}/MAG_phantom_summary.final.tsv}"
PHANTOM_MERYL_READ_DB="${PHANTOM_MERYL_READ_DB:-${PHANTOM_DIR}/meryl/raw_reads.k21.meryl}"
PHANTOM_MERYL_READ_CACHE="${PHANTOM_MERYL_READ_CACHE:-${PHANTOM_MERYL_READ_DB}.inputs.sha256}"
PHANTOM_KMER_FINAL="${PHANTOM_KMER_FINAL:-${PHANTOM_DIR}/phantom_kmer_confirmation.final.tsv}"
PHANTOM_KMER_CACHE_META="${PHANTOM_KMER_CACHE_META:-${PHANTOM_KMER_FINAL}.inputs.sha256}"
PHANTOM_MAG_MERYL_FINAL="${PHANTOM_MAG_MERYL_FINAL:-${PHANTOM_DIR}/MAG_phantom_meryl.final.tsv}"
PHANTOM_MAG_MERYL_CATALOGUE="${PHANTOM_MAG_MERYL_CATALOGUE:-${PHANTOM_DIR}/MAG_phantom_meryl.catalogue.final.tsv}"
PHANTOM_MAG_MERYL_CACHE_META="${PHANTOM_MAG_MERYL_CACHE_META:-${PHANTOM_MAG_MERYL_FINAL}.inputs.sha256}"
PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META="${PHANTOM_MAG_MERYL_CATALOGUE_CACHE_META:-${PHANTOM_MAG_MERYL_CATALOGUE}.inputs.sha256}"

REPEAT_SUM="${REPEAT_SUM:-${REPEAT_DIR}/contig_repeat_breadth.final.tsv}"
REPEAT_CACHE_META="${REPEAT_CACHE_META:-${REPEAT_SUM}.inputs.sha256}"
REPEAT_SELF_BLAST="${REPEAT_SELF_BLAST:-${REPEAT_DIR}/self_blast.paper_same_contig.filtered.tsv}"
REPEAT_MAG_SUM="${REPEAT_MAG_SUM:-${REPEAT_DIR}/MAG_repeat_breadth.final.tsv}"
REPEAT_DATASET_SUM="${REPEAT_DATASET_SUM:-${REPEAT_DIR}/dataset_repeat_breadth.final.tsv}"
REPEAT_REGION_BED="${REPEAT_REGION_BED:-${REPEAT_DIR}/repeat_regions.final.bed}"
REPEAT_REGION_TSV="${REPEAT_REGION_TSV:-${REPEAT_DIR}/repeat_regions.final.tsv}"

FALSE_DUP_CANDIDATE_PAF="${FALSE_DUP_CANDIDATE_PAF:-${FALSE_DUP_DIR}/duplication_self_minimap2.paf.gz}"
FALSE_DUP_CANDIDATE_PROVENANCE="${FALSE_DUP_CANDIDATE_PROVENANCE:-${FALSE_DUP_DIR}/duplication_candidate_search.provenance.tsv}"
# Assembly-wide outputs across every contig in the three source co-assemblies.
FALSE_DUP_ALL_PAIR_SUM="${FALSE_DUP_ALL_PAIR_SUM:-${FALSE_DUP_DIR}/false_duplication_pairs.all_contigs.final.tsv}"
FALSE_DUP_ALL_BOUNDARY_SUM="${FALSE_DUP_ALL_BOUNDARY_SUM:-${FALSE_DUP_DIR}/false_duplication_boundaries.all_contigs.final.tsv}"
FALSE_DUP_ALL_CONTIG_SUM="${FALSE_DUP_ALL_CONTIG_SUM:-${FALSE_DUP_DIR}/false_duplication.all_contigs.tsv}"
FALSE_DUP_ASSEMBLER_SUM="${FALSE_DUP_ASSEMBLER_SUM:-${FALSE_DUP_DIR}/false_duplication.by_assembler.tsv}"
# MAG-only outputs; the contig table retains compatibility fields used by 10_summary.sh.
FALSE_DUP_PAIR_SUM="${FALSE_DUP_PAIR_SUM:-${FALSE_DUP_DIR}/false_duplication_pairs.final.tsv}"
FALSE_DUP_BOUNDARY_SUM="${FALSE_DUP_BOUNDARY_SUM:-${FALSE_DUP_DIR}/false_duplication_boundaries.final.tsv}"
FALSE_DUP_CONTIG_SUM="${FALSE_DUP_CONTIG_SUM:-${FALSE_DUP_DIR}/false_duplication.by_contig.tsv}"
FALSE_DUP_MAG_SUM="${FALSE_DUP_MAG_SUM:-${FALSE_DUP_DIR}/false_duplication.by_MAG.tsv}"
FALSE_DUP_PROVENANCE="${FALSE_DUP_PROVENANCE:-${FALSE_DUP_DIR}/false_duplication.provenance.tsv}"
FALSE_DUP_CLIPPING_BOUNDARY_SUM="${FALSE_DUP_CLIPPING_BOUNDARY_SUM:-${FALSE_DUP_DIR}/validation/clipping/false_duplication.boundaries.clipping.tsv}"
FALSE_DUP_CLIPPING_PAIR_SUM="${FALSE_DUP_CLIPPING_PAIR_SUM:-${FALSE_DUP_DIR}/validation/clipping/false_duplication.pairs.clipping.tsv}"
FALSE_DUP_CLIPPING_CONTIG_SUM="${FALSE_DUP_CLIPPING_CONTIG_SUM:-${FALSE_DUP_DIR}/validation/clipping/false_duplication.by_contig.clipping.tsv}"
FALSE_DUP_CLIPPING_MAG_SUM="${FALSE_DUP_CLIPPING_MAG_SUM:-${FALSE_DUP_DIR}/validation/clipping/false_duplication.by_MAG.clipping.tsv}"
FALSE_DUP_CLIPPING_PROVENANCE="${FALSE_DUP_CLIPPING_PROVENANCE:-${FALSE_DUP_DIR}/validation/clipping/false_duplication.clipping.provenance.tsv}"

SUMMARY_CONTIG="${SUMMARY_CONTIG:-${SUMMARY_DIR}/integrated_QC_summary.total_contigs.tsv}"
SUMMARY_MAG_CONTIG="${SUMMARY_MAG_CONTIG:-${SUMMARY_DIR}/integrated_QC_summary.MAG_contigs.tsv}"
SUMMARY_MAG="${SUMMARY_MAG:-${SUMMARY_DIR}/integrated_QC_summary.by_MAG.tsv}"

mkdir -p "$OUT/00_log" "$PREP_DIR"
