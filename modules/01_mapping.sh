#!/usr/bin/env bash
set -Eeuo pipefail

# Map HiFi reads to each source-assembler co-assembly.
# Output BAMs preserve primary, secondary and supplementary alignments.

SCRIPT_ID="01_mapping"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

: "${OUT:?}" "${READS:?}" "${THREADS:?}" "${PRESET:?}"


MAP_DIR="$MAPPING_DIR"
mkdir -p "$MAP_DIR" "$OUT/00_log"
LOG="$OUT/00_log/01_mapping.log"
exec > >(tee -a "$LOG") 2>&1

MM2_THREADS="${MM2_THREADS:-$THREADS}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-50}"
(( SAMTOOLS_THREADS > THREADS )) && SAMTOOLS_THREADS="$THREADS"
SORT_MEM="${SORT_MEM:-2G}"
INDEX_CSI="${INDEX_CSI:-Y}"
# Independent assembler-specific pooled co-assemblies. BAMs, indexes and
# flagstat files are written directly to the mapping directory.
RUN_PER_ASSEMBLER_MAPPING="${RUN_PER_ASSEMBLER_MAPPING:-Y}"
ASSEMBLER_ASSEMBLY_DIR="${ASSEMBLER_ASSEMBLY_DIR:-}"
ASSEMBLER_FASTAS="${ASSEMBLER_FASTAS:-hifiasm_meta.fa metaflye.fa metaMDBG.fa}"
PER_ASSEMBLER_MAP_ROOT="${PER_ASSEMBLER_MAP_ROOT:-$MAP_DIR}"
HIFI_MAPPING_DIR="${HIFI_MAPPING_DIR:-$MAP_DIR}"
PER_ASSEMBLER_MANIFEST="${PER_ASSEMBLER_MANIFEST:-$ASSEMBLER_MAPPING_MANIFEST}"
FORCE_PER_ASSEMBLER_HIFI_REMAP="${FORCE_PER_ASSEMBLER_HIFI_REMAP:-N}"
FORCE_PER_ASSEMBLER_HIFI_PROFILES="${FORCE_PER_ASSEMBLER_HIFI_PROFILES:-N}"


MINIMAP2_INDEX_CHUNK="${MINIMAP2_INDEX_CHUNK:-256G}"
if [[ "${FORCE_SINGLE_PART_MINIMAP2:-Y}" == "Y" ]]; then
    MINIMAP2_INDEX_CHUNK="${MINIMAP2_INDEX_CHUNK_SINGLE:-900G}"
fi

TMP_PER_ASSEMBLER_FILES=()
TMP_PER_ASSEMBLER_PREFIXES=()
cleanup() {
    local path
    for path in "${TMP_PER_ASSEMBLER_FILES[@]}"; do
        [[ -n "$path" && -e "$path" ]] && rm -f -- "$path" || true
    done
    for path in "${TMP_PER_ASSEMBLER_PREFIXES[@]}"; do
        [[ -n "$path" ]] && rm -f -- "${path}".* 2>/dev/null || true
    done
}
trap cleanup EXIT


for tool in minimap2 samtools awk python3; do
    command -v "$tool" >/dev/null 2>&1 \
        || { echo "ERROR: missing tool: $tool"; exit 1; }
done

MINIMAP2_SECONDARY_SEQ_ARGS=()
_MINIMAP2_HELP_TEXT="$(minimap2 2>&1 || true)"
if grep -q -- '--secondary-seq' <<< "$_MINIMAP2_HELP_TEXT"; then
    MINIMAP2_SECONDARY_SEQ_ARGS=(--secondary-seq)
    MINIMAP2_SECONDARY_SEQ_STATUS="Y"
else
    MINIMAP2_SECONDARY_SEQ_STATUS="N"
fi
unset _MINIMAP2_HELP_TEXT
if [[ "$RUN_PER_ASSEMBLER_MAPPING" == "Y" ]]; then
    [[ -s "$READS" ]] || { echo "ERROR: missing READS: $READS"; exit 1; }
fi
for variable in \
    RUN_PER_ASSEMBLER_MAPPING \
    FORCE_PER_ASSEMBLER_HIFI_REMAP \
    FORCE_PER_ASSEMBLER_HIFI_PROFILES; do
    [[ "${!variable}" == "Y" || "${!variable}" == "N" ]] || {
        echo "ERROR: $variable must be Y or N" >&2
        exit 1
    }
done
if [[ "$RUN_PER_ASSEMBLER_MAPPING" == "Y" ]]; then
    [[ -d "$ASSEMBLER_ASSEMBLY_DIR" ]] || {
        echo "ERROR: assembler assembly directory missing: $ASSEMBLER_ASSEMBLY_DIR" >&2
        exit 1
    }
    read -r -a ASSEMBLY_FILES <<< "$ASSEMBLER_FASTAS"
    ((${#ASSEMBLY_FILES[@]} > 0)) || {
        echo "ERROR: ASSEMBLER_FASTAS is empty" >&2
        exit 1
    }
    declare -A SEEN_ASSEMBLERS=()
    for fasta_name in "${ASSEMBLY_FILES[@]}"; do
        [[ "$fasta_name" != */* ]] || {
            echo "ERROR: ASSEMBLER_FASTAS must contain basenames only: $fasta_name" >&2
            exit 1
        }
        [[ "$fasta_name" =~ ^[A-Za-z0-9_.-]+\.(fa|fasta|fna)$ ]] || {
            echo "ERROR: unsupported assembler FASTA filename: $fasta_name" >&2
            exit 1
        }
        assembler="${fasta_name%.*}"
        [[ -z "${SEEN_ASSEMBLERS[$assembler]+x}" ]] || {
            echo "ERROR: duplicate assembler label: $assembler" >&2
            exit 1
        }
        SEEN_ASSEMBLERS[$assembler]=1
        [[ -s "$ASSEMBLER_ASSEMBLY_DIR/$fasta_name" ]] || {
            echo "ERROR: assembler FASTA missing or empty: $ASSEMBLER_ASSEMBLY_DIR/$fasta_name" >&2
            exit 1
        }
    done
fi
[[ "${MINIMAP2_MAX_SECONDARY:-5}" == "5" ]] || {
    echo "ERROR: downstream ambiguity checks assume minimap2 default maximum secondary count of 5"
    exit 1
}

SAMTOOLS_VERSION="$(samtools --version | awk 'NR==1{print}')"

ensure_index() {
    local bam="$1"
    if [[ -s "$bam.bai" || -s "$bam.csi" ]]; then
        if samtools idxstats "$bam" >/dev/null 2>&1; then
            return 0
        fi
        echo "Removing invalid/stale BAM index: $bam"
        rm -f "$bam.bai" "$bam.csi"
    fi
    if [[ "$INDEX_CSI" == "Y" ]]; then
        samtools index -@ "$SAMTOOLS_THREADS" -c "$bam"
    else
        samtools index -@ "$SAMTOOLS_THREADS" "$bam" \
            || samtools index -@ "$SAMTOOLS_THREADS" -c "$bam"
    fi
}


per_assembler_bam_cache_current() {
    local bam="$1" meta="$2" fingerprint="$3" fai="$4"
    [[ -s "$bam" && -s "$meta" ]] \
        && [[ "$(tr -d '\r\n' < "$meta")" == "$fingerprint" ]] \
        && samtools quickcheck "$bam" >/dev/null 2>&1 \
        && qc_bam_matches_fai "$bam" "$fai"
}

remove_per_assembler_bam() {
    local path="$1"
    case "$path" in
        "$PER_ASSEMBLER_MAP_ROOT"/*)
            rm -f -- "$path" "$path.bai" "$path.csi"
            ;;
        *)
            echo "ERROR: refusing to remove output outside PER_ASSEMBLER_MAP_ROOT: $path" >&2
            exit 1
            ;;
    esac
}

finalize_per_assembler_bam() {
    local temporary="$1" final="$2" meta="$3" fingerprint="$4" fai="$5"
    samtools quickcheck -v "$temporary"
    qc_bam_matches_fai "$temporary" "$fai" || {
        echo "ERROR: assembler-specific BAM/reference mismatch: $temporary" >&2
        exit 1
    }
    mv -f -- "$temporary" "$final"
    rm -f -- "$final.bai" "$final.csi"
    ensure_index "$final"
    qc_write_cache_atomic "$meta" "$fingerprint"
}

build_per_assembler_hifi_profile() {
    local assembler="$1" source_bam="$2" fai="$3"
    local label="$4" final="$5" required_flags="$6" excluded_flags="$7"
    local source_fp="$8" meta="${final}.inputs.sha256"
    local fingerprint temporary

    fingerprint="$(qc_hash_text \
        "$source_fp" \
        "assembler=$assembler" \
        "samtools_view_required_flags=$required_flags" \
        "samtools_view_excluded_flags=$excluded_flags")"

    if [[ "$FORCE_PER_ASSEMBLER_HIFI_PROFILES" != "Y" ]] \
       && per_assembler_bam_cache_current \
            "$final" "$meta" "$fingerprint" "$fai"; then
        ensure_index "$final"
        echo "    Reusing HiFi $label BAM: $final"
        return 0
    fi

    echo "    Creating HiFi $label BAM"
    remove_per_assembler_bam "$final"
    rm -f -- "$meta"
    temporary="${final}.tmp.$$"
    TMP_PER_ASSEMBLER_FILES+=("$temporary")
    if (( required_flags > 0 )); then
        samtools view -@ "$SAMTOOLS_THREADS" -b \
            -f "$required_flags" -F "$excluded_flags" \
            -o "$temporary" "$source_bam"
    else
        samtools view -@ "$SAMTOOLS_THREADS" -b \
            -F "$excluded_flags" -o "$temporary" "$source_bam"
    fi
    finalize_per_assembler_bam \
        "$temporary" "$final" "$meta" "$fingerprint" "$fai"
}


# Independent mapping to each assembler-specific pooled co-assembly

if [[ "$RUN_PER_ASSEMBLER_MAPPING" == "Y" ]]; then
    mkdir -p "$PER_ASSEMBLER_MAP_ROOT" "$HIFI_MAPPING_DIR"
    manifest_tmp="${PER_ASSEMBLER_MANIFEST}.tmp.$$"
    TMP_PER_ASSEMBLER_FILES+=("$manifest_tmp")
    mkdir -p "$(dirname "$PER_ASSEMBLER_MANIFEST")"
    # Downstream modules take the HiFi reads from this column when --reads is omitted.
    READS_REALPATH="$(cd "$(dirname "$READS")" && pwd)/$(basename "$READS")"
    printf '%s\n' $'Assembler\tAssembly_FASTA\tHiFi_all_BAM\tHiFi_primary_BAM\tHiFi_secondary_BAM\tHiFi_supplementary_BAM\tHiFi_primary_supplementary_BAM\tOutput_directory\tReads' \
        > "$manifest_tmp"

    echo
    echo "[Per-assembler mapping] Started: $(date)"
    echo "Assembler assembly directory : $ASSEMBLER_ASSEMBLY_DIR"
    echo "Assembler FASTAs            : $ASSEMBLER_FASTAS"
    echo "Pooled HiFi reads           : $READS"
    echo "Output root                 : $PER_ASSEMBLER_MAP_ROOT"

    for fasta_name in "${ASSEMBLY_FILES[@]}"; do
        assembly="$ASSEMBLER_ASSEMBLY_DIR/$fasta_name"
        assembler="${fasta_name%.*}"
        assembler_dir="$PER_ASSEMBLER_MAP_ROOT"
        hifi_dir="$HIFI_MAPPING_DIR"
        mkdir -p "$hifi_dir"

        echo
        echo "[$assembler] Reference: $assembly"
        samtools faidx "$assembly"
        fai="${assembly}.fai"

        hifi_all="$hifi_dir/LR_to_${assembler}.mapped.sorted.bam"
        hifi_all_meta="${hifi_all}.inputs.sha256"
        hifi_primary="$hifi_dir/LR_to_${assembler}.primary.sorted.bam"
        hifi_secondary="$hifi_dir/LR_to_${assembler}.secondary.sorted.bam"
        hifi_supplementary="$hifi_dir/LR_to_${assembler}.supplementary.sorted.bam"
        hifi_primary_supplementary="$hifi_dir/LR_to_${assembler}.primary_supplementary.sorted.bam"
        hifi_flagstat="$hifi_dir/LR_to_${assembler}.mapped.flagstat.txt"

        hifi_reference_fp="$(qc_compute_mapping_fingerprint \
            "$assembly" "$READS" "$PRESET" "$MINIMAP2_INDEX_CHUNK")"
        hifi_mapping_fp="$(qc_hash_text \
            "$hifi_reference_fp" \
            "assembler=$assembler" \
            "minimap2=-ax_${PRESET}_-p_1_-N_${MINIMAP2_MAX_SECONDARY}_--secondary-seq" \
            "retained=mapped_records_only" \
            "sort=coordinate")"

        hifi_remapped="N"
        if [[ "$FORCE_PER_ASSEMBLER_HIFI_REMAP" != "Y" ]] \
           && per_assembler_bam_cache_current \
                "$hifi_all" "$hifi_all_meta" "$hifi_mapping_fp" "$fai"; then
            echo "    Reusing pooled HiFi mapped BAM"
            ensure_index "$hifi_all"
        elif [[ "$FORCE_PER_ASSEMBLER_HIFI_REMAP" != "Y" ]] \
             && qc_hifi_bam_mapping_semantically_matches \
                    "$hifi_all" "$assembly" "$READS" "$PRESET" \
                    "$MINIMAP2_INDEX_CHUNK" "$MINIMAP2_MAX_SECONDARY"; then
            echo "    Reusing semantically compatible pooled HiFi BAM; refreshing cache metadata only"
            ensure_index "$hifi_all"
            qc_write_cache_atomic "$hifi_all_meta" "$hifi_mapping_fp"
        else
            echo "    Mapping pooled HiFi reads"
            remove_per_assembler_bam "$hifi_all"
            rm -f -- "$hifi_all_meta" "$hifi_flagstat"
            hifi_tmp="${hifi_all}.tmp.$$"
            hifi_sort_prefix="$hifi_dir/.samtools_sort.LR_to_${assembler}.$$"
            TMP_PER_ASSEMBLER_FILES+=("$hifi_tmp")
            TMP_PER_ASSEMBLER_PREFIXES+=("$hifi_sort_prefix")
            minimap2 \
                -t "$MM2_THREADS" \
                -ax "$PRESET" \
                -p 1 \
                -N "$MINIMAP2_MAX_SECONDARY" \
                --secondary-seq \
                -I "$MINIMAP2_INDEX_CHUNK" \
                "$assembly" "$READS" \
              | samtools view -@ "$SAMTOOLS_THREADS" -b -F 4 - \
              | samtools sort -@ "$SAMTOOLS_THREADS" -m "$SORT_MEM" -O BAM \
                    -T "$hifi_sort_prefix" -o "$hifi_tmp" -
            finalize_per_assembler_bam \
                "$hifi_tmp" "$hifi_all" "$hifi_all_meta" \
                "$hifi_mapping_fp" "$fai"
            hifi_remapped="Y"
        fi

        build_per_assembler_hifi_profile \
            "$assembler" "$hifi_all" "$fai" \
            "primary-only" "$hifi_primary" 0 2308 "$hifi_mapping_fp"
        build_per_assembler_hifi_profile \
            "$assembler" "$hifi_all" "$fai" \
            "secondary-flagged" "$hifi_secondary" 256 4 "$hifi_mapping_fp"
        build_per_assembler_hifi_profile \
            "$assembler" "$hifi_all" "$fai" \
            "supplementary-flagged" "$hifi_supplementary" 2048 4 \
            "$hifi_mapping_fp"
        build_per_assembler_hifi_profile \
            "$assembler" "$hifi_all" "$fai" \
            "primary+supplementary" "$hifi_primary_supplementary" 0 260 \
            "$hifi_mapping_fp"

        if [[ "$hifi_remapped" == "Y" || ! -s "$hifi_flagstat" ]]; then
            hifi_flagstat_tmp="${hifi_flagstat}.tmp.$$"
            TMP_PER_ASSEMBLER_FILES+=("$hifi_flagstat_tmp")
            samtools flagstat -@ "$SAMTOOLS_THREADS" "$hifi_all" \
                > "$hifi_flagstat_tmp"
            mv -f -- "$hifi_flagstat_tmp" "$hifi_flagstat"
        fi


        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$assembler" "$assembly" "$hifi_all" "$hifi_primary" \
            "$hifi_secondary" "$hifi_supplementary" \
            "$hifi_primary_supplementary" "$assembler_dir" "$READS_REALPATH" \
            >> "$manifest_tmp"

        echo "    HiFi mapped BAM: $hifi_all"
    done

    mv -f -- "$manifest_tmp" "$PER_ASSEMBLER_MANIFEST"
    echo
    echo "[Per-assembler mapping] Finished: $(date)"
    echo "Mapping manifest: $PER_ASSEMBLER_MANIFEST"
else
    echo "Per-assembler mapping: disabled"
fi


echo "[01_mapping] Finished: $(date)"
