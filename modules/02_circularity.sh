#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_NAME="02_circularity.sh"
SCRIPT_ID="02_circularity"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"
SCRIPT_PATH="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"

usage() {
    cat <<'USAGE'
Validate END-to-START closure of single-contig long-read MAGs by mapping
HiFi reads to all contigs of the assembler that produced each MAG, with the
single-contig MAG contigs rotated.

Input and output options
  bash 02_circularity.sh MAG_DIR OUT_DIR   Direct positional input/output
  --assembly FILE              Optional single-assembler override
  --assembler-label TEXT       Label for --assembly single-assembler mode
  --assembler-assembly-dir DIR Directory containing the assembler FASTAs [configured]
  --assembler-fastas TEXT      Space-separated assembler FASTA basenames [configured]
  --mag-dir DIR                Directory containing final MAG FASTA files [configured]
  --hifi FILE                  HiFi FASTQ/FASTA; repeat for multiple files [configured]
  --outdir DIR                 New or empty output directory [configured]

Optional arguments
  --threads INT                Mapping/index-building threads [100]
  --sort-threads INT           samtools view/sort/index threads [20]
  --sort-memory STR            samtools sort memory per thread [2G]
  --hifi-anchor INT            Actual M/=X-aligned bp per original END|START side [2000]
  --min-hifi-reads INT         Minimum unique junction-spanning HiFi QNAMEs [2]
  --rerun-existing                Recompute stages even when their final outputs exist
  --help                       Show this message

Important interpretation
  Rotation is not circularity evidence. It relocates the original END|START
  adjacency to an internal reference coordinate without adding, deleting, or
  duplicating sequence. Only raw-read alignments across that coordinate are
  counted as closure evidence. The newly introduced linear reference break is
  never evaluated as an assembly error in this workflow.

  With no arguments, the script finds each single-contig MAG contig by exact sequence
  matching against hifiasm_meta.fa, metaflye.fa, and metaMDBG.fa. It then runs
  each assembler-specific group sequentially and merges the final results.

  --assembly is retained only for an optional single-assembler run. It must be
  one assembler's own assembly, never a concatenation of multiple assemblers.
USAGE
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    [[ -f "$1" ]] || die "File not found: $1"
    [[ -s "$1" ]] || die "File is empty: $1"
}

require_dir() {
    [[ -d "$1" ]] || die "Directory not found: $1"
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

all_nonempty() {
    local path
    for path in "$@"; do
        [[ -s "$path" ]] || return 1
    done
    return 0
}

tsv_has_column() {
    local path=$1
    local column=$2
    [[ -s "$path" ]] || return 1
    awk -F '\t' -v c="$column" 'NR==1 {for (i=1; i<=NF; i++) if ($i==c) found=1; exit} END {exit(found?0:1)}' "$path"
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

ASSEMBLY=""
ASSEMBLER_LABEL=""
WORKER_MODE="N"
DEFAULT_HIFI="$READS"
OUTDIR="$CIRCULARITY_ROTATED_RESULT_DIR"
THREADS="$CIRCULARITY_THREADS"
SORT_THREADS="$CIRCULARITY_SORT_THREADS"
SORT_MEMORY="$CIRCULARITY_SORT_MEMORY"
HIFI_ANCHOR="$CIRCULARITY_HIFI_ANCHOR_BP"
MIN_HIFI_READS="$CIRCULARITY_MIN_HIFI_QNAMES"
CIRCULARITY_FINALIZE="$SCRIPT_DIR/../scripts/circularity_finalize.py"
REUSE_EXISTING="Y"
MINIMAP2_INDEX_CHUNK="$CIRCULARITY_MINIMAP2_INDEX_CHUNK"
MINIMAP2_MAX_SECONDARY="$CIRCULARITY_MINIMAP2_MAX_SECONDARY"
HIFI_FILES=()

if (( $# >= 2 )) && [[ "$1" != -* ]]; then
    shift 2
fi

while (($#)); do
    case "$1" in
        --assembly) ASSEMBLY=${2:-}; shift 2 ;;
        --assembler-label) ASSEMBLER_LABEL=${2:-}; shift 2 ;;
        --assembler-assembly-dir) ASSEMBLER_ASSEMBLY_DIR=${2:-}; shift 2 ;;
        --assembler-fastas) ASSEMBLER_FASTAS=${2:-}; shift 2 ;;
        --mag-dir) MAG_DIR=${2:-}; shift 2 ;;
        --hifi) HIFI_FILES+=("${2:-}"); shift 2 ;;
        --outdir) OUTDIR=${2:-}; shift 2 ;;
        --threads) THREADS=${2:-}; shift 2 ;;
        --sort-threads) SORT_THREADS=${2:-}; shift 2 ;;
        --sort-memory) SORT_MEMORY=${2:-}; shift 2 ;;
        --hifi-anchor) HIFI_ANCHOR=${2:-}; shift 2 ;;
        --min-hifi-reads) MIN_HIFI_READS=${2:-}; shift 2 ;;
        --rerun-existing) REUSE_EXISTING="N"; shift ;;
        --group-worker) WORKER_MODE="Y"; shift ;;
        --help|-h) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

if ((${#HIFI_FILES[@]} == 0)); then
    HIFI_FILES=("$DEFAULT_HIFI")
fi

# Keep one module-level log for all assembler groups. Workers inherit the
# parent's stdout/stderr so they must not reopen or overwrite this file.
if [[ "$WORKER_MODE" == "N" ]]; then
    mkdir -p "$OUT/00_log"
    LOG="$OUT/00_log/02_circularity.log"
    exec > >(tee "$LOG") 2>&1
fi

if [[ "$WORKER_MODE" == "N" ]]; then
    CIRCULARITY_TSV="${CIRCULARITY_TSV:-$OUT/00_prep/assembler_circularity.MAG.tsv}"
    if [[ ! -s "$CIRCULARITY_TSV" ]]; then
        die "Missing assembler circularity table: $CIRCULARITY_TSV. Run MAGverify prep-circularity first."
    fi
fi
require_dir "$MAG_DIR"
require_file "$CIRCULARITY_FINALIZE"
for file in "${HIFI_FILES[@]}"; do require_file "$file"; done

for value in "$THREADS" "$SORT_THREADS" "$HIFI_ANCHOR" "$MIN_HIFI_READS"; do
    is_positive_integer "$value" || die "Expected a positive integer, observed: $value"
done

require_command python3
require_command samtools

PREFERRED_MINIMAP2_BIN="${MINIMAP2_BIN:-}"
if [[ -x "$PREFERRED_MINIMAP2_BIN" ]]; then
    MINIMAP2_CMD="$PREFERRED_MINIMAP2_BIN"
else
    MINIMAP2_CMD="$(command -v minimap2 2>/dev/null || true)"
fi
[[ -n "$MINIMAP2_CMD" && -x "$MINIMAP2_CMD" ]] \
    || die "Required command not found: minimap2"
MINIMAP2_VERSION="$($MINIMAP2_CMD --version 2>/dev/null | head -n 1 || true)"
[[ -n "$MINIMAP2_VERSION" ]] || MINIMAP2_VERSION="unknown"

echo "[$SCRIPT_NAME] minimap2 executable: $MINIMAP2_CMD"
echo "[$SCRIPT_NAME] minimap2 version: $MINIMAP2_VERSION"

migrate_circularity_layout() {
    local cache="$OUTDIR/.cache" pair old new
    mkdir -p "$cache"
    for pair in \
        "00_source_plan|.cache/assembler_plan" "00_assembler_plan|.cache/assembler_plan" \
        "01_by_source|.cache/by_assembler" "01_by_assembler|.cache/by_assembler" \
        "MAG_END_START_closure.final.tsv|.cache/circularity.MAG.junction.tsv" \
        "circularity.MAG.junction.tsv|.cache/circularity.MAG.junction.tsv" \
        "MAG_circularity.audit.tsv|.cache/circularity.MAG.audit.tsv" \
        "circularity.MAG.audit.tsv|.cache/circularity.MAG.audit.tsv" \
        "run_parameters.txt|.cache/run_parameters.txt" \
        "source_group_outputs.tsv|.cache/assembler_group_outputs.tsv" \
        "assembler_group_outputs.tsv|.cache/assembler_group_outputs.tsv"; do
        old="$OUTDIR/${pair%%|*}"
        new="$OUTDIR/${pair#*|}"
        if [[ -e "$old" && ! -e "$new" ]]; then
            mv "$old" "$new"
        elif [[ -e "$old" && -f "$old" ]]; then
            rm -f "$old"
        fi
    done
    rm -f "$OUTDIR/MAG_circularity.final.tsv"
}

run_assembler_groups() {
    migrate_circularity_layout
    local plan_dir="$OUTDIR/.cache/assembler_plan"
    local grouped_mag_root="$plan_dir/group_MAGs"
    local plan_table="$plan_dir/assembler_execution_plan.tsv"
    local group_table="$plan_dir/assembler_groups.tsv"
    local validation_table="$plan_dir/assembler_matching_validation.tsv"
    local grouped_output_root="$OUTDIR/.cache/by_assembler"

    require_dir "$ASSEMBLER_ASSEMBLY_DIR"
    local -a assembler_fasta_names=()
    read -r -a assembler_fasta_names <<< "$ASSEMBLER_FASTAS"
    ((${#assembler_fasta_names[@]} > 0)) || die "Configured assembler FASTA list is empty"
    local fasta_name
    for fasta_name in "${assembler_fasta_names[@]}"; do
        [[ "$fasta_name" != */* ]] || die "Assembler FASTA must be a basename: $fasta_name"
        require_file "$ASSEMBLER_ASSEMBLY_DIR/$fasta_name"
    done

    local reuse_plan="N"
    if [[ "$REUSE_EXISTING" == "Y" ]] && all_nonempty "$plan_table" "$group_table" "$validation_table"; then
        reuse_plan="Y"
        local check_group_id check_assembler check_assembly check_count check_mag_dir
        while IFS=$'\t' read -r check_group_id check_assembler check_assembly check_count check_mag_dir; do
            [[ "$check_group_id" == "Group_ID" ]] && continue
            if [[ ! -d "$check_mag_dir" ]]; then
                reuse_plan="N"
                break
            fi
        done < "$group_table"
    fi

    if [[ "$reuse_plan" == "Y" ]]; then
        echo "[$SCRIPT_NAME] SKIP: existing MAG-to-assembler matching plan is complete."
    else
        rm -rf "$plan_dir"
        mkdir -p "$plan_dir" "$grouped_mag_root"
        echo "[$SCRIPT_NAME] Matching each single-contig MAG to its contig in the assembler FASTAs."

        python3 "$SCRIPT_DIR/../scripts/circularity_reference.py" match "$ASSEMBLER_ASSEMBLY_DIR" "$ASSEMBLER_FASTAS" "$MAG_DIR" \
            "$grouped_mag_root" "$plan_table" "$group_table" \
            "$validation_table" "$HIFI_ANCHOR"
    fi

    mkdir -p "$grouped_output_root"
    local group_id assembler assembly mag_count group_mag_dir
    local n_groups worker_jobs worker_threads worker_sort_threads
    n_groups="$(awk -F '\t' 'NR > 1 && $1 != "" {n++} END {print n+0}' "$group_table")"
    worker_jobs=1
    worker_threads="$(qc_split_threads "$THREADS" "$worker_jobs")"
    worker_sort_threads="$(qc_split_threads "$SORT_THREADS" "$worker_jobs")"
    local worker_sort_memory="$SORT_MEMORY" sort_mb
    if [[ -n "${MEMORY_GB:-}" ]]; then
        sort_mb=$(( MEMORY_GB * 1024 / 2 / (worker_jobs * worker_sort_threads) ))
        (( sort_mb < 256 )) && sort_mb=256
        worker_sort_memory="${sort_mb}M"
    fi
    echo "[$SCRIPT_NAME] Assembler groups run one at a time (threads $worker_threads, sort threads $worker_sort_threads x $worker_sort_memory)"
    while IFS=$'\t' read -r group_id assembler assembly mag_count group_mag_dir; do
        [[ "$group_id" == "Group_ID" ]] && continue
        [[ -n "$group_id" && -n "$assembly" && -n "$group_mag_dir" ]] \
            || die "Malformed row in $group_table"
        echo "[$SCRIPT_NAME] Running assembler group: $assembler ($mag_count MAGs)"
        local -a worker_command=(
            bash "$SCRIPT_PATH"
            --group-worker
            --assembly "$assembly"
            --assembler-label "$assembler"
            --mag-dir "$group_mag_dir"
            --hifi-anchor "$HIFI_ANCHOR"
            --min-hifi-reads "$MIN_HIFI_READS"
            --outdir "$grouped_output_root/$group_id"
            --threads "$worker_threads"
            --sort-threads "$worker_sort_threads"
            --sort-memory "$worker_sort_memory"
        )
        if [[ "$REUSE_EXISTING" == "N" ]]; then
            worker_command+=(--rerun-existing)
        fi
        local hifi_file
        for hifi_file in "${HIFI_FILES[@]}"; do
            worker_command+=(--hifi "$hifi_file")
        done
        mkdir -p "$grouped_output_root/$group_id"
        rm -f "$grouped_output_root/$group_id/.worker_done"
        ( "${worker_command[@]}" && : > "$grouped_output_root/$group_id/.worker_done" ) &
        qc_wait_for_slot "$worker_jobs"
    done < "$group_table"
    qc_wait_all
    while IFS=$'\t' read -r group_id assembler assembly mag_count group_mag_dir; do
        [[ "$group_id" == "Group_ID" || -z "$group_id" ]] && continue
        [[ -e "$grouped_output_root/$group_id/.worker_done" ]] \
            || die "circularity failed for assembler group $group_id ($assembler)"
    done < "$group_table"

    local merged_final="$OUTDIR/.cache/circularity.MAG.junction.tsv"
    local merged_parameters="$OUTDIR/.cache/run_parameters.txt"
    local group_outputs="$OUTDIR/.cache/assembler_group_outputs.tsv"

    python3 "$CIRCULARITY_FINALIZE" merge "$plan_table" "$group_table" "$grouped_output_root" \
        "$validation_table" "$merged_final" "$group_outputs"

    local circularity_compat_final="$OUTDIR/.cache/circularity.MAG.tsv"
    local circularity_compat_audit="$OUTDIR/.cache/circularity.MAG.audit.tsv"
    python3 "$CIRCULARITY_FINALIZE" refresh \
        --canonical "$merged_final" \
        --assembler-circularity "$CIRCULARITY_TSV" \
        --final "$circularity_compat_final" \
        --audit "$circularity_compat_audit"
    python3 "$CIRCULARITY_FINALIZE" add-junction \
        "$circularity_compat_final" "$merged_final"
    local assembler_contig_calls=""
    if [[ "$CIRCULARITY_TSV" == *.MAG.tsv && -s "${CIRCULARITY_TSV%.MAG.tsv}.contig.tsv" ]]; then
        assembler_contig_calls="${CIRCULARITY_TSV%.MAG.tsv}.contig.tsv"
    fi
    python3 "$CIRCULARITY_FINALIZE" contig \
        --mag-table "$circularity_compat_final" \
        --mag-dir "$MAG_DIR" \
        --assembly-dir "$ASSEMBLER_ASSEMBLY_DIR" \
        --assembly-fastas "$ASSEMBLER_FASTAS" \
        --assembler-circularity "$assembler_contig_calls" \
        --output "$OUTDIR/.cache/circularity.contig.tsv"

    {
        echo "Script: $SCRIPT_NAME"
        echo "Module: $SCRIPT_ID"
        echo "Execution_mode: automatic_assembler_matched_by_exact_contig_sequence"
        echo "Assembler_assembly_directory: $ASSEMBLER_ASSEMBLY_DIR"
        echo "Assembler_FASTAs: $ASSEMBLER_FASTAS"
        echo "Assembler_execution_plan: $plan_table"
        echo "Assembler_matching_validation: $validation_table"
        echo "MAG_directory: $MAG_DIR"
        echo "Integrated_final_output: $merged_final"
                echo "Assembler_group_outputs: $group_outputs"
        echo "Output_layout: final_files_at_module_root; intermediate_analyses_in_numbered_directories"
        printf 'HiFi_files:'; printf ' %s' "${HIFI_FILES[@]}"; echo
        echo "Reuse_existing_completed_stages: $REUSE_EXISTING"
        echo "Junction_support_policy: INS_allowed_DN_crossing_excluded_PS_main_v1"
        echo "Reference_policy: exact_assembler_total_assembly_only"
        echo "Multi_assembler_concatenated_reference_used: N"
        echo "Analyzed_MAG_scope: single_contig_only"
        echo "Circularity_analyzable_rule: Y only for single-contig MAGs with >= minimum HiFi QNAMEs reaching both END and START anchors; otherwise N"
        echo "Final_output_scope: all_MAGs; multi_contig_MAGs_reported_as_non_analyzable_NA_rows"
        echo "Mapping_threads: $THREADS"
        echo "Samtools_sort_threads: $SORT_THREADS"
        echo "Samtools_sort_memory_per_thread: $SORT_MEMORY"
        echo "HiFi_anchor_bp_per_side: $HIFI_ANCHOR"
        echo "Minimum_junction_spanning_HiFi_reads: $MIN_HIFI_READS"
        echo "HiFi_main_alignment_classes: primary_plus_supplementary"
        echo "HiFi_sensitivity_alignment_classes: all_mapped_primary_secondary_supplementary"
        echo "HiFi_minimap2_secondary_policy: --secondary=yes_-p_1_-N_${MINIMAP2_MAX_SECONDARY}_--secondary-seq"
        echo "HiFi_minimap2_executable: $MINIMAP2_CMD"
        echo "HiFi_minimap2_version: $MINIMAP2_VERSION"
        echo "HiFi_MAPQ_cutoff_used: N"
                        echo "Circularity_supported_rule: CIRCULARITY_SUPPORTED equals HIFI_PRIMARY_SUPPLEMENTARY_END_START_SUPPORTED_Y_N; ALL-mapped and secondary-only contributions are sensitivity/audit outputs"
        echo "MAG_circularity_evaluable_rule: Circularity_evaluable mirrors canonical CIRCULARITY_ANALYZABLE"
        echo "MAG_read_supported_topology_rule: Read_supported_circular_topology preserves canonical Y/N/NA"
        echo "High_confidence_circular_rule: assembler-reported circularity=Y and read-supported circular topology=Y"
        echo "Additional_circular_candidate_rule: assembler-reported circularity=N and read-supported circular topology=Y"
        echo "MAG_Circular_compatibility_alias: Circular equals High_confidence_circular"
        } > "$merged_parameters"

    echo "[$SCRIPT_NAME] PASS: all assembler groups completed."
    echo "Final merged summary: $merged_final"
    echo "Assembler matching validation: $validation_table"
}

if [[ -e "$OUTDIR" ]] && find "$OUTDIR" -mindepth 1 -print -quit 2>/dev/null | grep -q .; then
    if [[ "$REUSE_EXISTING" == "Y" ]]; then
        echo "[$SCRIPT_NAME] Existing output directory detected; completed stages will be skipped: $OUTDIR"
    else
        echo "[$SCRIPT_NAME] --rerun-existing requested; existing stage outputs may be replaced: $OUTDIR"
    fi
else
    mkdir -p "$OUTDIR"
fi

if [[ "$WORKER_MODE" == "N" && -z "$ASSEMBLY" ]]; then
    run_assembler_groups
    exit 0
fi

[[ -n "$ASSEMBLY" ]] || die "Internal group-worker error: --assembly is empty"
require_file "$ASSEMBLY"
if [[ -z "$ASSEMBLER_LABEL" ]]; then
    ASSEMBLER_LABEL="$(basename "$ASSEMBLY")"
    ASSEMBLER_LABEL="${ASSEMBLER_LABEL%.gz}"
    ASSEMBLER_LABEL="${ASSEMBLER_LABEL%.*}"
fi

REFERENCE_DIR="$OUTDIR/00_rotated_reference"
VALIDATION_ROOT="$OUTDIR/validation"
HIFI_DIR="$VALIDATION_ROOT/hifi"
mkdir -p "$REFERENCE_DIR" "$HIFI_DIR"

JUNCTION_POLICY_VERSION="INS_allowed_DN_crossing_excluded_PS_main_indel_metrics"
JUNCTION_POLICY_MARKER="$OUTDIR/.junction_support_policy"
JUNCTION_POLICY_CURRENT="N"
if [[ "$REUSE_EXISTING" == "Y" && -s "$JUNCTION_POLICY_MARKER" ]] && \
   [[ "$(cat "$JUNCTION_POLICY_MARKER")" == "$JUNCTION_POLICY_VERSION" ]]; then
    JUNCTION_POLICY_CURRENT="Y"
fi
if [[ "$REUSE_EXISTING" == "Y" && "$JUNCTION_POLICY_CURRENT" == "N" ]]; then
    echo "[$SCRIPT_NAME] Junction-support policy changed: existing mappings will be reused, but END|START evidence will be recalculated."
fi

ROTATED_REFERENCE="$REFERENCE_DIR/assembly.rotated_MAG_contigs.fa"
ROTATION_MANIFEST="$REFERENCE_DIR/rotation_manifest.tsv"
for legacy_pair in \
    "$REFERENCE_DIR/full_assembly_with_rotated_single_contig_MAGs.fa|$ROTATED_REFERENCE" \
    "$REFERENCE_DIR/full_assembly_with_rotated_single_contig_MAGs.fa.fai|$ROTATED_REFERENCE.fai" \
    "$REFERENCE_DIR/full_assembly_with_rotated_single_contig_MAGs.mmi|$REFERENCE_DIR/assembly.rotated_MAG_contigs.mmi" \
    "$HIFI_DIR/hifi_to_rotated_full_assembly.sorted.bam|$HIFI_DIR/hifi_to_assembly.rotated_MAG_contigs.sorted.bam" \
    "$HIFI_DIR/hifi_to_rotated_full_assembly.sorted.bam.bai|$HIFI_DIR/hifi_to_assembly.rotated_MAG_contigs.sorted.bam.bai" \
    "$HIFI_DIR/hifi_to_rotated_full_assembly.sorted.bam.csi|$HIFI_DIR/hifi_to_assembly.rotated_MAG_contigs.sorted.bam.csi"; do
    legacy_file="${legacy_pair%%|*}"
    current_file="${legacy_pair#*|}"
    if [[ -e "$legacy_file" && ! -e "$current_file" ]]; then
        mv "$legacy_file" "$current_file"
    fi
done
CANDIDATE_SELECTION="$REFERENCE_DIR/candidate_selection.tsv"
CANDIDATE_BED="$REFERENCE_DIR/single_contig_MAG_contigs.bed"
JUNCTION_BED="$REFERENCE_DIR/original_END_START_junction_regions.bed"

ROTATED_REFERENCE_REUSED="N"
if [[ "$REUSE_EXISTING" == "Y" ]] && all_nonempty \
    "$ROTATED_REFERENCE" "$ROTATION_MANIFEST" "$CANDIDATE_SELECTION" \
    "$CANDIDATE_BED" "$JUNCTION_BED"; then
    ROTATED_REFERENCE_REUSED="Y"
    echo "[$SCRIPT_NAME] SKIP: rotated reference and manifests already exist."
else
    echo "[$SCRIPT_NAME] Building the reference: all $ASSEMBLER_LABEL contigs, single-contig MAG contigs rotated."

    python3 "$SCRIPT_DIR/../scripts/circularity_reference.py" rotate "$ASSEMBLY" "$MAG_DIR" "$ROTATED_REFERENCE" \
    "$ROTATION_MANIFEST" "$CANDIDATE_SELECTION" "$CANDIDATE_BED" \
    "$JUNCTION_BED" "$HIFI_ANCHOR"
fi

if [[ "$ROTATED_REFERENCE_REUSED" == "Y" ]]; then
    if [[ -s "${ROTATED_REFERENCE}.fai" ]]; then
        echo "[$SCRIPT_NAME] SKIP: rotated-reference FASTA index already exists."
    else
        echo "[$SCRIPT_NAME] SKIP: rotated reference is reused; the cached circularity outputs do not require rebuilding the .fai."
    fi
else
    samtools faidx "$ROTATED_REFERENCE"
fi

HIFI_INDEX="$REFERENCE_DIR/assembly.rotated_MAG_contigs.mmi"
HIFI_BAM="$HIFI_DIR/hifi_to_assembly.rotated_MAG_contigs.sorted.bam"
HIFI_LOG="$HIFI_DIR/minimap2.log"
HIFI_SPAN_DETAILS="$HIFI_DIR/hifi_spanning_read_details.tsv"
HIFI_SUMMARY="$HIFI_DIR/hifi_END_START_summary.tsv"
HIFI_EVIDENCE_COMPLETE="N"
if [[ "$REUSE_EXISTING" == "Y" && "$JUNCTION_POLICY_CURRENT" == "Y" ]] && all_nonempty "$HIFI_SPAN_DETAILS" "$HIFI_SUMMARY"; then
    HIFI_EVIDENCE_COMPLETE="Y"
fi

if [[ "$HIFI_EVIDENCE_COMPLETE" == "Y" ]]; then
    echo "[$SCRIPT_NAME] SKIP: existing HiFi evidence is complete; minimap2 HiFi index is not needed."
elif [[ "$REUSE_EXISTING" == "Y" && -s "$HIFI_INDEX" ]]; then
    echo "[$SCRIPT_NAME] SKIP: minimap2 rotated-reference index already exists."
else
    echo "[$SCRIPT_NAME] Building one-part minimap2 index for global competitive mapping."
    if ! "$MINIMAP2_CMD" -t "$THREADS" -I "$MINIMAP2_INDEX_CHUNK" -d "$HIFI_INDEX" "$ROTATED_REFERENCE" \
        2> "$REFERENCE_DIR/minimap2_index.log"; then
        echo "ERROR: minimap2 index construction failed. Last log lines:" >&2
        tail -n 50 "$REFERENCE_DIR/minimap2_index.log" >&2 || true
        exit 1
    fi
fi

if [[ "$HIFI_EVIDENCE_COMPLETE" == "Y" ]]; then
    echo "[$SCRIPT_NAME] SKIP: existing HiFi evidence is complete; HiFi mapping is not repeated."
elif [[ "$REUSE_EXISTING" == "Y" && -s "$HIFI_BAM" ]]; then
    echo "[$SCRIPT_NAME] SKIP: existing HiFi rotated-reference BAM will be reused."
else
    echo "[$SCRIPT_NAME] Mapping HiFi reads to all $ASSEMBLER_LABEL contigs (MAG contigs rotated)."
    if ! python3 "$SCRIPT_DIR/../scripts/core_progress_reader.py" "${HIFI_FILES[@]}" \
        --label "$SCRIPT_NAME][$ASSEMBLER_LABEL mapping" --interval "${PROGRESS_INTERVAL:-300}" |
        "$MINIMAP2_CMD" -t "$THREADS" -ax map-hifi \
        --secondary=yes \
        -p 1 \
        -N "$MINIMAP2_MAX_SECONDARY" \
        --secondary-seq \
        "$HIFI_INDEX" - 2> "$HIFI_LOG" |
        samtools view -@ "$SORT_THREADS" -b -F 4 - |
        samtools sort -@ "$SORT_THREADS" -m "$SORT_MEMORY" -o "$HIFI_BAM" -; then
        echo "ERROR: HiFi minimap2/samtools mapping pipeline failed. Last minimap2 log lines:" >&2
        tail -n 50 "$HIFI_LOG" >&2 || true
        exit 1
    fi
fi
if [[ "$HIFI_EVIDENCE_COMPLETE" == "Y" ]]; then
    echo "[$SCRIPT_NAME] SKIP: existing HiFi evidence is complete; HiFi BAM index is not needed."
elif [[ "$REUSE_EXISTING" == "Y" && ( -s "${HIFI_BAM}.bai" || -s "${HIFI_BAM}.csi" ) ]]; then
    echo "[$SCRIPT_NAME] SKIP: HiFi BAM index already exists."
else
    samtools index -@ "$SORT_THREADS" "$HIFI_BAM"
fi

if [[ "$REUSE_EXISTING" == "Y" && "$JUNCTION_POLICY_CURRENT" == "Y" ]] && all_nonempty "$HIFI_SPAN_DETAILS" "$HIFI_SUMMARY"; then
    echo "[$SCRIPT_NAME] SKIP: existing HiFi END|START evidence tables will be reused."
else
    echo "[$SCRIPT_NAME] Detecting junction-spanning HiFi alignments across the internalized original END|START."

    python3 "$SCRIPT_DIR/../scripts/circularity_junction.py" summarize "$HIFI_BAM" "$JUNCTION_BED" "$ROTATION_MANIFEST" \
    "$HIFI_ANCHOR" "$MIN_HIFI_READS" \
    "$HIFI_SPAN_DETAILS" "$HIFI_SUMMARY"
fi

FINAL_SUMMARY="$OUTDIR/circularity.MAG.junction.tsv"
METHODS_SUMMARY="$OUTDIR/run_parameters.txt"

if [[ "$REUSE_EXISTING" == "Y" && "$JUNCTION_POLICY_CURRENT" == "Y" ]] && \
   tsv_has_column "$FINAL_SUMMARY" "HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N" && \
   tsv_has_column "$FINAL_SUMMARY" "HIFI_END_START_SUPPORTED" && \
   tsv_has_column "$FINAL_SUMMARY" "CIRCULARITY_SUPPORTED"; then
    echo "[$SCRIPT_NAME] SKIP: HiFi-only final circularity table is current."
else
    echo "[$SCRIPT_NAME] Building HiFi-only final circularity table."
    python3 "$SCRIPT_DIR/../scripts/circularity_junction.py" integrate "$CANDIDATE_SELECTION" "$ROTATION_MANIFEST" "$HIFI_SUMMARY" "$FINAL_SUMMARY"
fi

{
    echo "Script: $SCRIPT_NAME"
    echo "Module: $SCRIPT_ID"
    echo "Execution_mode: assembler_group_worker"
    echo "Assembler_label: $ASSEMBLER_LABEL"
    echo "Assembly: $ASSEMBLY"
    echo "MAG_directory: $MAG_DIR"
    printf 'HiFi_files:'; printf ' %s' "${HIFI_FILES[@]}"; echo
    echo "Reuse_existing_completed_stages: $REUSE_EXISTING"
    echo "Junction_support_policy: $JUNCTION_POLICY_VERSION"
    echo "Integrated_final_output: $FINAL_SUMMARY"
    echo "Output_layout: final_files_at_assembler_group_root; intermediate_analyses_in_numbered_directories"
    echo "Mapping_threads: $THREADS"
    echo "Samtools_sort_threads: $SORT_THREADS"
    echo "Samtools_sort_memory_per_thread: $SORT_MEMORY"
    echo "HiFi_anchor_bp_per_side: $HIFI_ANCHOR"
    echo "Minimum_junction_spanning_HiFi_reads: $MIN_HIFI_READS"
    echo "Reference_definition: all contigs of the assembler that produced the tested MAG contigs"
    echo "Analyzed_MAG_scope: single_contig_only"
    echo "Circularity_analyzable_rule: Y only for single-contig MAGs with >= minimum HiFi QNAMEs reaching both END and START anchors; otherwise N"
    echo "Final_output_scope: all_MAGs_present_in_worker_MAG_directory; multi_contig_MAGs_reported_as_non_analyzable_NA_rows"
    echo "HiFi_main_alignment_classes: primary_plus_supplementary"
    echo "HiFi_sensitivity_alignment_classes: all_mapped_primary_secondary_supplementary"
    echo "HiFi_minimap2_secondary_policy: --secondary=yes_-p_1_-N_${MINIMAP2_MAX_SECONDARY}_--secondary-seq"
    echo "HiFi_minimap2_executable: $MINIMAP2_CMD"
    echo "HiFi_minimap2_version: $MINIMAP2_VERSION"
    echo "HiFi_MAPQ_cutoff_used: N"
    echo "HiFi_support_definition: main call uses unique QNAMEs from primary+supplementary alignments with >= required actual M/=X-aligned reference bases on both sides of internalized END|START; no HiFi MAPQ cutoff; junction insertion is allowed; D/N crossing END|START invalidates that alignment as direct support"
    echo "HiFi_all_mapped_sensitivity_definition: same topology rule with primary+secondary+supplementary; reported separately and does not replace the P+S main call"
    echo "HiFi_secondary_contribution_rule: HIFI_SECONDARY_ADDED_END_START_SUPPORT_Y_N=Y when ALL-mapped support is Y but P+S main support is not Y"
    echo "Clean_junction_definition: HIFI_PRIMARY_SUPPLEMENTARY_CLEAN_JUNCTION_Y_N=Y when >= $MIN_HIFI_READS distinct P+S topology-supporting HiFi reads have no I/D/N operation crossing END|START; annotation only"
    echo "Junction_event_audit: I/D/N crossing counts and maximum event lengths are reported separately for topology-supporting HiFi reads"
    echo "Proper_pair_only_support_used: N"
    echo "Circularity_supported_rule: CIRCULARITY_SUPPORTED equals HIFI_PRIMARY_SUPPLEMENTARY_END_START_SUPPORTED_Y_N; ALL-mapped and secondary-only contributions are sensitivity/audit outputs"
    echo "Rotation_definition: coordinate change only; no base added, removed, or duplicated"
    echo "New_reference_break_policy: excluded from closure and assembly-error assessment"
} > "$METHODS_SUMMARY"
printf '%s\n' "$JUNCTION_POLICY_VERSION" > "$JUNCTION_POLICY_MARKER"

echo "[$SCRIPT_NAME] PASS"
echo "Final summary: $FINAL_SUMMARY"
echo "HiFi analysis summary: $HIFI_SUMMARY"
echo "Rotation manifest: $ROTATION_MANIFEST"
