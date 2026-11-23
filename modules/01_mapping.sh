#!/usr/bin/env bash
set -Eeuo pipefail

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
TMP_PER_ASSEMBLER_DIRS=()
cleanup() {
    local path
    for path in "${TMP_PER_ASSEMBLER_FILES[@]}"; do
        [[ -n "$path" && -e "$path" ]] && rm -f -- "$path" || true
    done
    for path in "${TMP_PER_ASSEMBLER_PREFIXES[@]}"; do
        [[ -n "$path" ]] && rm -f -- "${path}".* 2>/dev/null || true
    done
    for path in "${TMP_PER_ASSEMBLER_DIRS[@]}"; do
        [[ -n "$path" && -d "$path" ]] && rm -rf -- "$path" || true
    done
}
trap cleanup EXIT

for tool in minimap2 samtools awk python3; do
    command -v "$tool" >/dev/null 2>&1 \
        || { echo "ERROR: missing tool: $tool"; exit 1; }
done

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
    mkdir -p "$(dirname "$meta")"
    qc_write_cache_atomic "$meta" "$fingerprint"
}

build_per_assembler_hifi_profile() {
    local assembler="$1" source_bam="$2" fai="$3"
    local label="$4" final="$5" required_flags="$6" excluded_flags="$7"
    local source_fp="$8" meta
    meta="$(qc_bam_cache "$final")"
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
        echo "[$assembler] Reusing HiFi $label BAM: $final"
        return 0
    fi

    echo "[$assembler] Creating HiFi $label BAM"
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

rename_legacy_mapping_outputs() {
    local dir="$1" assembler="$2" kind old new suffix
    mkdir -p "$dir/.cache"
    for kind in mapped primary primary_supplementary; do
        old="$dir/LR_to_${assembler}.${kind}.sorted.bam"
        new="$dir/${assembler}.${kind}.bam"
        if [[ -e "$old" && ! -e "$new" ]]; then
            for suffix in "" .csi .bai .inputs.sha256; do
                [[ -e "$old$suffix" ]] && mv -f -- "$old$suffix" "$new$suffix"
            done
            echo "[$assembler] Renamed $(basename "$old") -> $(basename "$new")"
        fi
        if [[ -e "$new.inputs.sha256" ]]; then
            mv -f -- "$new.inputs.sha256" "$(qc_bam_cache "$new")"
        fi
    done
    new="$dir/.cache/${assembler}.flagstat.txt"
    for old in "$dir/LR_to_${assembler}.mapped.flagstat.txt" "$dir/${assembler}.flagstat.txt"; do
        if [[ -e "$old" && ! -e "$new" ]]; then
            mkdir -p "$dir/.cache"
            mv -f -- "$old" "$new"
        fi
        rm -f -- "$old"
    done
    for kind in secondary supplementary; do
        for old in "$dir/LR_to_${assembler}.${kind}.sorted.bam" "$dir/${assembler}.${kind}.bam"; do
            rm -f -- "$old" "$old.csi" "$old.bai" "$old.inputs.sha256"
        done
    done
    return 0
}

map_one_assembler() {
    local fasta_name="$1" row_file="$2"
    local assembly assembler assembler_dir hifi_dir fai hifi_all hifi_all_meta
    local hifi_primary hifi_primary_supplementary
    local hifi_flagstat hifi_reference_fp hifi_mapping_fp hifi_remapped
    local hifi_tmp hifi_sort_prefix hifi_flagstat_tmp
    assembly="$ASSEMBLER_ASSEMBLY_DIR/$fasta_name"
    assembler="${fasta_name%.*}"
    assembler_dir="$PER_ASSEMBLER_MAP_ROOT"
    hifi_dir="$HIFI_MAPPING_DIR"
    mkdir -p "$hifi_dir"

    echo
    echo "[$assembler] Reference: $assembly"
    samtools faidx "$assembly"
    fai="${assembly}.fai"

    hifi_all="$hifi_dir/${assembler}.mapped.bam"
    hifi_all_meta="$(qc_bam_cache "$hifi_all")"
    hifi_primary="$hifi_dir/${assembler}.primary.bam"
    hifi_primary_supplementary="$hifi_dir/${assembler}.primary_supplementary.bam"
    hifi_flagstat="$hifi_dir/.cache/${assembler}.flagstat.txt"
    mkdir -p "$hifi_dir/.cache"
    rename_legacy_mapping_outputs "$hifi_dir" "$assembler"

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
        echo "[$assembler] Reusing pooled HiFi mapped BAM"
        ensure_index "$hifi_all"
    elif [[ "$FORCE_PER_ASSEMBLER_HIFI_REMAP" != "Y" ]] \
         && qc_hifi_bam_mapping_semantically_matches \
                "$hifi_all" "$assembly" "$READS" "$PRESET" \
                "$MINIMAP2_INDEX_CHUNK" "$MINIMAP2_MAX_SECONDARY"; then
        echo "[$assembler] Reusing semantically compatible pooled HiFi BAM; refreshing cache metadata only"
        ensure_index "$hifi_all"
        mkdir -p "$(dirname "$hifi_all_meta")"
        qc_write_cache_atomic "$hifi_all_meta" "$hifi_mapping_fp"
    else
        echo "[$assembler] Mapping pooled HiFi reads"
        remove_per_assembler_bam "$hifi_all"
        rm -f -- "$hifi_all_meta" "$hifi_flagstat"
        hifi_tmp="${hifi_all}.tmp.$$"
        hifi_sort_prefix="$hifi_dir/.samtools_sort.${assembler}.$$"
        TMP_PER_ASSEMBLER_FILES+=("$hifi_tmp")
        TMP_PER_ASSEMBLER_PREFIXES+=("$hifi_sort_prefix")
        python3 "$SCRIPT_DIR/../scripts/core_progress_reader.py" "$READS" \
            --label "01_mapping][$assembler" --interval "${PROGRESS_INTERVAL:-300}" \
          | minimap2 \
            -t "$MM2_THREADS" \
            -ax "$PRESET" \
            -p 1 \
            -N "$MINIMAP2_MAX_SECONDARY" \
            --secondary-seq \
            -I "$MINIMAP2_INDEX_CHUNK" \
            "$assembly" - \
          | samtools view -@ "$SAMTOOLS_THREADS" -b -F 4 - \
          | samtools sort -@ "$SAMTOOLS_THREADS" -m "$SORT_MEM" -O BAM \
                -T "$hifi_sort_prefix" -o "$hifi_tmp" -
        finalize_per_assembler_bam \
            "$hifi_tmp" "$hifi_all" "$hifi_all_meta" \
            "$hifi_mapping_fp" "$fai"
        hifi_remapped="Y"
    fi

    echo "[$assembler] Mapping done; building alignment-class BAMs while the next assembler maps"
    (
        TMP_PER_ASSEMBLER_FILES=()
        TMP_PER_ASSEMBLER_PREFIXES=()
        TMP_PER_ASSEMBLER_DIRS=()
        trap cleanup EXIT
        SAMTOOLS_THREADS=$(( SAMTOOLS_THREADS / 2 )); (( SAMTOOLS_THREADS < 1 )) && SAMTOOLS_THREADS=1
        local profile_pids=() pid
        build_per_assembler_hifi_profile \
            "$assembler" "$hifi_all" "$fai" \
            "primary-only" "$hifi_primary" 0 2308 "$hifi_mapping_fp" & profile_pids+=($!)
        build_per_assembler_hifi_profile \
            "$assembler" "$hifi_all" "$fai" \
            "primary+supplementary" "$hifi_primary_supplementary" 0 260 \
            "$hifi_mapping_fp" & profile_pids+=($!)
        for pid in "${profile_pids[@]}"; do
            wait "$pid" || { echo "ERROR: [$assembler] alignment-class BAM failed" >&2; exit 1; }
        done

        if [[ "$hifi_remapped" == "Y" || ! -s "$hifi_flagstat" ]]; then
            hifi_flagstat_tmp="${hifi_flagstat}.tmp.$$"
            TMP_PER_ASSEMBLER_FILES+=("$hifi_flagstat_tmp")
            samtools flagstat -@ "$SAMTOOLS_THREADS" "$hifi_all" \
                > "$hifi_flagstat_tmp"
            mv -f -- "$hifi_flagstat_tmp" "$hifi_flagstat"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$assembler" "$assembly" "$hifi_all" "$hifi_primary" \
            "$hifi_primary_supplementary" "$assembler_dir" "$READS_REALPATH" \
            > "$row_file"
    ) &
    POST_PIDS+=($!)

    echo "[$assembler] Done: $hifi_all"
}

if [[ "$RUN_PER_ASSEMBLER_MAPPING" == "Y" ]]; then
    mkdir -p "$PER_ASSEMBLER_MAP_ROOT" "$HIFI_MAPPING_DIR"
    for stale in "${PER_ASSEMBLER_MANIFEST}".tmp.*; do
        [[ -e "$stale" ]] || continue
        kill -0 "${stale##*.}" 2>/dev/null || rm -f -- "$stale"
    done
    manifest_tmp="${PER_ASSEMBLER_MANIFEST}.tmp.$$"
    TMP_PER_ASSEMBLER_FILES+=("$manifest_tmp")
    mkdir -p "$(dirname "$PER_ASSEMBLER_MANIFEST")"
    READS_REALPATH="$(cd "$(dirname "$READS")" && pwd)/$(basename "$READS")"
    printf '%s\n' $'Assembler\tAssembly_FASTA\tHiFi_all_BAM\tHiFi_primary_BAM\tHiFi_primary_supplementary_BAM\tOutput_directory\tReads' \
        > "$manifest_tmp"

    echo
    echo "[Per-assembler mapping] Started: $(date)"
    echo "Assembler assembly directory : $ASSEMBLER_ASSEMBLY_DIR"
    echo "Assembler FASTAs            : $ASSEMBLER_FASTAS"
    echo "Pooled HiFi reads           : $READS"
    echo "Output root                 : $PER_ASSEMBLER_MAP_ROOT"

    if [[ -n "${MEMORY_GB:-}" ]]; then
        sort_mb=$(( MEMORY_GB * 1024 / 2 / SAMTOOLS_THREADS ))
        (( sort_mb < 256 )) && sort_mb=256
        SORT_MEM="${sort_mb}M"
    fi
    echo "Assemblers mapped one at a time (minimap2 -t $MM2_THREADS, samtools -@ $SAMTOOLS_THREADS -m $SORT_MEM); alignment-class BAMs built in parallel"

    ROW_DIR="$(mktemp -d "$PER_ASSEMBLER_MAP_ROOT/.manifest_rows.XXXXXX")"
    TMP_PER_ASSEMBLER_DIRS+=("$ROW_DIR")
    POST_PIDS=()
    job_index=0
    for fasta_name in "${ASSEMBLY_FILES[@]}"; do
        job_index=$(( job_index + 1 ))
        map_one_assembler "$fasta_name" "$ROW_DIR/$(printf '%03d' "$job_index").tsv"
    done
    failed=0
    for pid in "${POST_PIDS[@]}"; do
        wait "$pid" || failed=1
    done
    job_index=0
    for fasta_name in "${ASSEMBLY_FILES[@]}"; do
        job_index=$(( job_index + 1 ))
        row_file="$ROW_DIR/$(printf '%03d' "$job_index").tsv"
        [[ -s "$row_file" ]] || { echo "ERROR: mapping failed for $fasta_name" >&2; failed=1; continue; }
        cat "$row_file" >> "$manifest_tmp"
    done
    (( failed == 0 )) || exit 1

    python3 - "$manifest_tmp" "$HIFI_MAPPING_DIR/.cache" <<'PY_FLAGSTAT'
import os
import re
import sys

manifest, flagstat_dir = sys.argv[1:]
labels = [
    ('in total (QC-passed reads + QC-failed reads)', 'Alignments_total'),
    ('primary', 'Primary'), ('secondary', 'Secondary'), ('supplementary', 'Supplementary'),
    ('duplicates', 'Duplicates'), ('mapped', 'Mapped'), ('primary mapped', 'Primary_mapped'),
]
columns = [name for _, name in labels] + ['Mapped_percent', 'Primary_mapped_percent']
with open(manifest, encoding='utf-8') as handle:
    lines = handle.read().splitlines()
out = [lines[0] + '\t' + '\t'.join(columns)]
for line in lines[1:]:
    if not line:
        continue
    assembler = line.split('\t')[0]
    values = {name: 'NA' for name in columns}
    path = os.path.join(flagstat_dir, f'{assembler}.flagstat.txt')
    if os.path.isfile(path):
        with open(path, encoding='utf-8') as stats:
            for raw in stats:
                match = re.match(r'(\d+) \+ \d+ (.+?)(?: \((.*?)%? :.*\))?$', raw.strip())
                if not match:
                    continue
                count, label, percent = match.groups()
                for text, name in labels:
                    if label == text:
                        values[name] = count
                        if name == 'Mapped' and percent:
                            values['Mapped_percent'] = percent.rstrip('%')
                        if name == 'Primary_mapped' and percent:
                            values['Primary_mapped_percent'] = percent.rstrip('%')
    out.append(line + '\t' + '\t'.join(values[name] for name in columns))
with open(manifest, 'w', encoding='utf-8') as handle:
    handle.write('\n'.join(out) + '\n')
PY_FLAGSTAT
    mv -f -- "$manifest_tmp" "$PER_ASSEMBLER_MANIFEST"
    echo
    echo "[Per-assembler mapping] Finished: $(date)"
    echo "Mapping manifest: $PER_ASSEMBLER_MANIFEST"
else
    echo "Per-assembler mapping: disabled"
fi

echo "[01_mapping] Finished: $(date)"
