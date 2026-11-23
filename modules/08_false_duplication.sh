#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_ID="08_false_duplication"
CANDIDATE_ALGORITHM_VERSION="false_duplication_candidates_nested_representative"
CLIPPING_ANNOTATION_VERSION="false_duplication_clipping_annotation_v1"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

FALSE_DUP_MINIMAP_PRESET="${FALSE_DUP_MINIMAP_PRESET:-asm5}"
FALSE_DUP_MINIMAP_INDEX_CHUNK="${FALSE_DUP_MINIMAP_INDEX_CHUNK:-${MINIMAP2_INDEX_CHUNK_SINGLE:-900G}}"
FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP="${FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP:-0.90}"

: "${OUT:?ERROR: OUT is not set}"
: "${MAG_CONTIG_MAP:?ERROR: MAG_CONTIG_MAP is not set}"
: "${FALSE_DUP_MAPPING_MANIFEST:?ERROR: FALSE_DUP_MAPPING_MANIFEST is not set}"
: "${FALSE_DUP_CANDIDATE_PAF:?ERROR: FALSE_DUP_CANDIDATE_PAF is not set}"
: "${FALSE_DUP_CANDIDATE_PROVENANCE:?ERROR: FALSE_DUP_CANDIDATE_PROVENANCE is not set}"
: "${FALSE_DUP_ALL_PAIR_SUM:?ERROR: FALSE_DUP_ALL_PAIR_SUM is not set}"
: "${FALSE_DUP_ALL_BOUNDARY_SUM:?ERROR: FALSE_DUP_ALL_BOUNDARY_SUM is not set}"
: "${FALSE_DUP_ALL_CONTIG_SUM:?ERROR: FALSE_DUP_ALL_CONTIG_SUM is not set}"
: "${FALSE_DUP_ASSEMBLER_SUM:?ERROR: FALSE_DUP_ASSEMBLER_SUM is not set}"
: "${FALSE_DUP_PAIR_SUM:?ERROR: FALSE_DUP_PAIR_SUM is not set}"
: "${FALSE_DUP_BOUNDARY_SUM:?ERROR: FALSE_DUP_BOUNDARY_SUM is not set}"
: "${FALSE_DUP_CONTIG_SUM:?ERROR: FALSE_DUP_CONTIG_SUM is not set}"
: "${FALSE_DUP_MAG_SUM:?ERROR: FALSE_DUP_MAG_SUM is not set}"
: "${FALSE_DUP_PROVENANCE:?ERROR: FALSE_DUP_PROVENANCE is not set}"
: "${FALSE_DUP_CLIPPING_BOUNDARY_SUM:?ERROR: FALSE_DUP_CLIPPING_BOUNDARY_SUM is not set}"
: "${FALSE_DUP_CLIPPING_PAIR_SUM:?ERROR: FALSE_DUP_CLIPPING_PAIR_SUM is not set}"
: "${FALSE_DUP_CLIPPING_CONTIG_SUM:?ERROR: FALSE_DUP_CLIPPING_CONTIG_SUM is not set}"
: "${FALSE_DUP_CLIPPING_MAG_SUM:?ERROR: FALSE_DUP_CLIPPING_MAG_SUM is not set}"
: "${FALSE_DUP_CLIPPING_PROVENANCE:?ERROR: FALSE_DUP_CLIPPING_PROVENANCE is not set}"

LOG="$OUT/00_log/08_false_duplication.log"
FALSE_DUP_IDN_EVENTS="${FALSE_DUP_CACHE_DIR}/false_duplication_IDN.events.tsv"
FALSE_DUP_IDN_SUMMARY="${FALSE_DUP_CACHE_DIR}/false_duplication_IDN.summary.tsv"
FALSE_DUP_CANDIDATE_PREPROCESSING="${FALSE_DUP_CACHE_DIR}/false_duplication_candidates.preprocessing.tsv"

mkdir -p "$OUT/00_log" "${FALSE_DUP_DIR}" "$(dirname "$FALSE_DUP_CLIPPING_MAG_SUM")"

for pair in \
    "duplication_self_minimap2.paf.gz|$FALSE_DUP_CANDIDATE_PAF" \
    "duplication_candidate_search.provenance.tsv|$FALSE_DUP_CANDIDATE_PROVENANCE"; do
    old="${FALSE_DUP_DIR}/${pair%%|*}"; new="${pair#*|}"
    if [[ -e "$old" && ! -e "$new" ]]; then mkdir -p "$(dirname "$new")"; mv "$old" "$new"; fi
done
for legacy in false_duplication_HiFi_sensitivity.tsv false_duplication_HiFi_anchor_sensitivity.tsv \
    false_duplication_HiFi_anchor_sensitivity.pairs.tsv duplication_self_minimap2.paf.gz \
    duplication_candidate_search.provenance.tsv false_duplication_pairs.all_contigs.final.tsv \
    false_duplication_boundaries.all_contigs.final.tsv false_duplication.all_contigs.tsv \
    false_duplication.by_assembler.tsv false_duplication_pairs.final.tsv \
    false_duplication_boundaries.final.tsv false_duplication.provenance.tsv \
    false_duplication_IDN.events.tsv false_duplication_IDN.summary.tsv \
    false_duplication_candidates.preprocessing.tsv; do
    rm -f "${FALSE_DUP_DIR}/$legacy"
done
rm -rf "${FALSE_DUP_DIR}/validation"

exec > >(tee -a "$LOG") 2>&1

TMP_ROOT=""
cleanup() { [[ -z "$TMP_ROOT" || ! -d "$TMP_ROOT" ]] || rm -rf "$TMP_ROOT"; }
on_error() {
    local rc=$?
    echo "ERROR: 08_false_duplication.sh failed at line ${BASH_LINENO[0]}: ${BASH_COMMAND} (exit ${rc})"
    exit "$rc"
}
trap cleanup EXIT
trap on_error ERR

for tool in python3 samtools minimap2 gzip sha256sum stat awk cat; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: required tool not found: $tool"; exit 1; }
done

DETECTOR="$SCRIPT_DIR/../scripts/false_duplication.py"
[[ -s "$DETECTOR" ]] || { echo "ERROR: detector missing: $DETECTOR"; exit 1; }
CLIPPING_ANNOTATOR="$SCRIPT_DIR/../scripts/false_duplication.py"
[[ -s "$CLIPPING_ANNOTATOR" ]] || { echo "ERROR: clipping annotator missing: $CLIPPING_ANNOTATOR"; exit 1; }
[[ -s "$MAG_CONTIG_MAP" ]] || { echo "ERROR: MAG-contig map missing: $MAG_CONTIG_MAP"; exit 1; }
[[ -s "$FALSE_DUP_MAPPING_MANIFEST" ]] || {
    echo "ERROR: assembler mapping manifest missing: $FALSE_DUP_MAPPING_MANIFEST"
    echo "Run the current 01_mapping.sh first."
    exit 1
}

FALSE_DUP_CLIPPING_INPUT="NA"
if [[ "${RUN_CLIPPING:-N}" == "Y" ]]; then
    [[ -s "$CLIP_FINAL" ]] || {
        echo "ERROR: validated clipping result missing: $CLIP_FINAL"
        echo "Run the current 03_clipping.sh before 08_false_duplication.sh."
        exit 1
    }
    qc_validate_current_clipping_result || exit 1
    FALSE_DUP_CLIPPING_INPUT="$CLIP_FINAL"
fi

for name in \
    FALSE_DUP_MIN_COPY_BP FALSE_DUP_HIFI_ANCHOR_BP \
    FALSE_DUP_HIFI_MIN_SPANNING_READS FALSE_DUP_HIFI_MIN_REACHING_READS; do
    value="${!name}"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: $name must be a positive integer: $value"; exit 1; }
done
python3 - "$FALSE_DUP_MIN_PIDENT" "$FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP" <<'PY'
import sys
minimum = float(sys.argv[1])
redundancy = float(sys.argv[2])
if not 0 < minimum <= 100:
    raise SystemExit('ERROR: FALSE_DUP_MIN_PIDENT must be >0 and <=100')
if not 0 < redundancy <= 1:
    raise SystemExit('ERROR: FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP must be >0 and <=1')
PY

manifest_records() {
    python3 - "$FALSE_DUP_MAPPING_MANIFEST" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as h:
    r = csv.DictReader(h, delimiter='\t')
    fields = set(r.fieldnames or [])
    required = {'Assembler', 'Assembly_FASTA', 'HiFi_all_BAM'}
    missing = required - fields
    if missing:
        raise SystemExit('ERROR: mapping manifest lacks columns: ' + ','.join(sorted(missing)))
    rows = list(r)
if not rows:
    raise SystemExit('ERROR: mapping manifest contains no assembler rows')
for row in rows:
    print('\t'.join([
        row['Assembler'].strip(),
        row['Assembly_FASTA'].strip(),
        row['HiFi_all_BAM'].strip(),
    ]))
PY
}

ensure_bam_index() {
    local bam="$1"
    if [[ ! -s "${bam}.bai" && ! -s "${bam}.csi" ]]; then
        echo "Creating BAM index: $bam"
        samtools index -@ "$THREADS" "$bam"
    fi
}

while IFS=$'\t' read -r assembler assembly hifi_bam; do
    [[ -n "$assembler" && -s "$assembly" && -s "${assembly}.fai" ]] || {
        echo "ERROR: assembler FASTA/FAI missing for $assembler: $assembly"; exit 1;
    }
    [[ -s "$hifi_bam" ]] || { echo "ERROR: HiFi BAM missing for $assembler: $hifi_bam"; exit 1; }
    samtools quickcheck -v "$hifi_bam" || { echo "ERROR: HiFi BAM failed quickcheck: $hifi_bam"; exit 1; }
    qc_bam_matches_fai "$hifi_bam" "${assembly}.fai" || { echo "ERROR: HiFi BAM/reference mismatch for $assembler"; exit 1; }
    ensure_bam_index "$hifi_bam"
done < <(manifest_records)

echo "[08_false_duplication] Started: $(date)"
echo "SCRIPT_ID: $SCRIPT_ID"
echo "Candidate self-alignment: minimap2 -x ${FALSE_DUP_MINIMAP_PRESET} -DP"
echo "Candidate definition: alignment >=${FALSE_DUP_MIN_COPY_BP} bp; identity >=${FALSE_DUP_MIN_PIDENT}%"
echo "Candidate preprocessing: intrinsic same-contig copy overlap excluded"
echo "Redundancy collapse: nested at both loci -> outer representative; partial overlap requires reciprocal overlap >=${FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP} at both loci"
echo "Intra-contig classification: tandem/interspersed; direct/inverted orientation recorded"
echo "Tandem definition: exact copy adjacency"
echo "Direct tandem validation: two outer boundaries supported + shared internal junction unsupported => false duplication"
echo "Inter-contig candidates: evaluated only when both contigs belong to the same MAG"
echo "HiFi structural main: primary + supplementary; anchor=${FALSE_DUP_HIFI_ANCHOR_BP} bp; min spanning=${FALSE_DUP_HIFI_MIN_SPANNING_READS}; min reaching/side=${FALSE_DUP_HIFI_MIN_REACHING_READS}"
echo "Clipping annotation: validated HiFi clipping is reported as co-occurrence only"

CANDIDATE_FP="$({
    while IFS=$'\t' read -r assembler assembly hifi_bam; do
        printf '%s\n' "$assembler"
        stat -c '%n\t%s\t%Y' "$assembly"
        sha256sum "${assembly}.fai"
    done < <(manifest_records)
    minimap2 --version
    printf '%s\n' \
        "preset=${FALSE_DUP_MINIMAP_PRESET}" \
        "index_chunk=${FALSE_DUP_MINIMAP_INDEX_CHUNK}" \
        "$CANDIDATE_ALGORITHM_VERSION"
} | sha256sum | awk '{print $1}')"

CANDIDATE_REUSE="N"
if [[ "${FORCE_RERUN_FALSE_DUPLICATION:-N}" == "N" && -s "$FALSE_DUP_CANDIDATE_PAF" && -s "$FALSE_DUP_CANDIDATE_PROVENANCE" ]]; then
    OLD="$(awk -F '\t' '$1=="fingerprint"{print $2; exit}' "$FALSE_DUP_CANDIDATE_PROVENANCE")"
    [[ "$OLD" == "$CANDIDATE_FP" ]] && CANDIDATE_REUSE="Y"
fi

if [[ "$CANDIDATE_REUSE" == "Y" ]]; then
    echo "Reusing Minimap2 assembly self-alignment."
else
    TMP_ROOT="$(mktemp -d "${FALSE_DUP_DIR}/.candidate.tmp.XXXXXX")"
    paf_parts=()
    self_rows=()
    while IFS=$'\t' read -r assembler assembly hifi_bam; do
        self_rows+=("$assembler"$'\t'"$assembly")
    done < <(manifest_records)
    self_jobs="$(qc_parallel_jobs "${#self_rows[@]}")"
    self_threads="$(qc_split_threads "$THREADS" "$self_jobs")"
    for row in "${self_rows[@]}"; do
        IFS=$'\t' read -r assembler assembly <<< "$row"
        safe="$(printf '%s' "$assembler" | tr -cs 'A-Za-z0-9_.-' '_')"
        part="$TMP_ROOT/${safe}.self.paf.gz"
        echo "Self-aligning $assembler with Minimap2 (${self_threads} threads; ${self_jobs} assemblers in parallel)"
        (
            minimap2 -x "$FALSE_DUP_MINIMAP_PRESET" -DP \
                -t "$self_threads" -I "$FALSE_DUP_MINIMAP_INDEX_CHUNK" \
                "$assembly" "$assembly" 2> "$TMP_ROOT/${safe}.minimap2.log" \
                | gzip -c > "$part.tmp"
            mv -f "$part.tmp" "$part"
        ) &
        qc_wait_for_slot "$self_jobs"
        paf_parts+=("$part")
    done
    qc_wait_all
    for part in "${paf_parts[@]}"; do
        [[ -s "$part" ]] || { echo "ERROR: Minimap2 self-alignment output missing: $part"; exit 1; }
    done
    TMP_PAF="$TMP_ROOT/duplication_self_minimap2.paf.gz"
    cat "${paf_parts[@]}" > "$TMP_PAF"
    mv -f "$TMP_PAF" "$FALSE_DUP_CANDIDATE_PAF"

    TMP_PROV="${FALSE_DUP_CANDIDATE_PROVENANCE}.tmp.$$"
    {
        printf 'key\tvalue\n'
        printf 'fingerprint\t%s\n' "$CANDIDATE_FP"
        printf 'candidate_method\t%s\n' 'Minimap2 assembly self-alignment'
        printf 'minimap2_preset\t%s\n' "$FALSE_DUP_MINIMAP_PRESET"
        printf 'minimap2_options\t%s\n' '-DP'
        printf 'index_chunk\t%s\n' "$FALSE_DUP_MINIMAP_INDEX_CHUNK"
        printf 'scope\t%s\n' 'all intra-contig self-alignments plus same-MAG inter-contig alignments; direct and inverted orientations parsed downstream'
    } > "$TMP_PROV"
    mv -f "$TMP_PROV" "$FALSE_DUP_CANDIDATE_PROVENANCE"
    rm -rf "$TMP_ROOT"
    TMP_ROOT=""
fi

FP_INPUT="$({
    sha256sum "$FALSE_DUP_CANDIDATE_PAF" "$MAG_CONTIG_MAP" "$DETECTOR"
    printf '%s  %s\n' "$(cut -f1-7 "$FALSE_DUP_MAPPING_MANIFEST" | sha256sum | cut -d' ' -f1)" "$FALSE_DUP_MAPPING_MANIFEST"
    while IFS=$'\t' read -r assembler assembly hifi_bam; do
        sha256sum "${assembly}.fai"
        stat -c '%n\t%s\t%Y' "$hifi_bam"
    done < <(manifest_records)
    printf '%s\n' \
        "min_copy_bp=${FALSE_DUP_MIN_COPY_BP}" \
        "min_pident=${FALSE_DUP_MIN_PIDENT}" \
        "redundant_reciprocal_overlap=${FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP}" \
        "hifi_anchor=${FALSE_DUP_HIFI_ANCHOR_BP}" \
        "hifi_min_spanning=${FALSE_DUP_HIFI_MIN_SPANNING_READS}" \
        "hifi_min_reaching=${FALSE_DUP_HIFI_MIN_REACHING_READS}" \
        "$SCRIPT_ID"
} | sha256sum | awk '{print $1}')"

REUSE="N"
if [[ "${FORCE_RERUN_FALSE_DUPLICATION:-N}" == "N" \
      && -s "$FALSE_DUP_ALL_PAIR_SUM" && -s "$FALSE_DUP_ALL_BOUNDARY_SUM" \
      && -s "$FALSE_DUP_ALL_CONTIG_SUM" && -s "$FALSE_DUP_ASSEMBLER_SUM" \
      && -s "$FALSE_DUP_PAIR_SUM" && -s "$FALSE_DUP_BOUNDARY_SUM" \
      && -s "$FALSE_DUP_CONTIG_SUM" && -s "$FALSE_DUP_MAG_SUM" \
      && -s "$FALSE_DUP_IDN_EVENTS" && -s "$FALSE_DUP_IDN_SUMMARY" \
      && -s "$FALSE_DUP_CANDIDATE_PREPROCESSING" \
      && -s "$FALSE_DUP_PROVENANCE" ]]; then
    OLD="$(awk -F '\t' '$1=="fingerprint"{print $2; exit}' "$FALSE_DUP_PROVENANCE")"
    [[ "$OLD" == "$FP_INPUT" ]] && REUSE="Y"
fi

if [[ "$REUSE" == "Y" ]]; then
    echo "Reusing current false-duplication results."
else
    python3 "$DETECTOR" detect \
        --candidate-paf "$FALSE_DUP_CANDIDATE_PAF" \
        --mapping-manifest "$FALSE_DUP_MAPPING_MANIFEST" \
        --mag-map "$MAG_CONTIG_MAP" \
        --min-copy-bp "$FALSE_DUP_MIN_COPY_BP" \
        --min-pident "$FALSE_DUP_MIN_PIDENT" \
        --redundant-reciprocal-overlap "$FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP" \
        --hifi-anchor-bp "$FALSE_DUP_HIFI_ANCHOR_BP" \
        --hifi-min-spanning-reads "$FALSE_DUP_HIFI_MIN_SPANNING_READS" \
        --hifi-min-reaching-reads "$FALSE_DUP_HIFI_MIN_REACHING_READS" \
        --out-all-pairs "$FALSE_DUP_ALL_PAIR_SUM" \
        --out-all-boundaries "$FALSE_DUP_ALL_BOUNDARY_SUM" \
        --out-all-contigs "$FALSE_DUP_ALL_CONTIG_SUM" \
        --out-assemblers "$FALSE_DUP_ASSEMBLER_SUM" \
        --out-pairs "$FALSE_DUP_PAIR_SUM" \
        --out-boundaries "$FALSE_DUP_BOUNDARY_SUM" \
        --out-contigs "$FALSE_DUP_CONTIG_SUM" \
        --out-mags "$FALSE_DUP_MAG_SUM" \
        --out-idn-events "$FALSE_DUP_IDN_EVENTS" \
        --out-idn-summary "$FALSE_DUP_IDN_SUMMARY" \
        --out-candidate-preprocessing "$FALSE_DUP_CANDIDATE_PREPROCESSING"

    TMP_PROV="${FALSE_DUP_PROVENANCE}.tmp.$$"
    {
        printf 'key\tvalue\n'
        printf 'fingerprint\t%s\n' "$FP_INPUT"
        printf 'script_id\t%s\n' "$SCRIPT_ID"
        printf 'candidate_policy\t%s\n' "Minimap2 asm5 -DP; alignment>=${FALSE_DUP_MIN_COPY_BP}bp; identity>=${FALSE_DUP_MIN_PIDENT}%"
        printf 'candidate_arrangements\t%s\n' 'assembly-level analysis uses intra-contig tandem/interspersed candidates; same-MAG inter-contig candidates are additionally evaluated at MAG level; direct/inverted orientations retained'
        printf 'overlap_preprocessing_policy\t%s\n' 'same-contig self-alignments whose two copy intervals intrinsically overlap are excluded from the false-duplication candidate denominator'
        printf 'redundancy_policy\t%s\n' "exact reciprocal duplicates reduced; nested alignments at both duplicated loci collapse to the wider outer representative; partial overlaps collapse only within the same assembler/contig pair/orientation/scope when both loci have reciprocal overlap >=${FALSE_DUP_REDUNDANT_RECIPROCAL_OVERLAP}"
        printf 'inter_contig_policy\t%s\n' 'only contig pairs assigned to the same MAG are structurally evaluated'
        printf 'hifi_main_policy\t%s\n' "primary+supplementary alignments; unique QNAME per boundary; anchor=${FALSE_DUP_HIFI_ANCHOR_BP}; min_spanning=${FALSE_DUP_HIFI_MIN_SPANNING_READS}; min_reaching_per_side=${FALSE_DUP_HIFI_MIN_REACHING_READS}"
        printf 'tandem_policy\t%s\n' 'exactly adjacent intra-contig copies are tandem; direct and inverted orientations are reported separately; direct tandem false duplication requires supported outer boundaries and an unsupported shared internal junction; inverted tandem remains unresolved unless its internal junction is directly supported'
        printf 'final_call_policy\t%s\n' 'HiFi Y/N/NA; positive if a false-duplication pattern is established; negative only when the relevant candidate boundaries are resolved without a false-duplication pattern; otherwise unresolved'
    } > "$TMP_PROV"
    mv -f "$TMP_PROV" "$FALSE_DUP_PROVENANCE"
fi

python3 "$CLIPPING_ANNOTATOR" annotate \
    --pairs "$FALSE_DUP_ALL_PAIR_SUM" \
    --boundaries "$FALSE_DUP_ALL_BOUNDARY_SUM" \
    --contigs "$FALSE_DUP_ALL_CONTIG_SUM" \
    --clipping "$FALSE_DUP_CLIPPING_INPUT" \
    --out-boundaries "$FALSE_DUP_CLIPPING_BOUNDARY_SUM" \
    --out-pairs "$FALSE_DUP_CLIPPING_PAIR_SUM" \
    --out-contigs "$FALSE_DUP_CLIPPING_CONTIG_SUM" \
    --out-mags "$FALSE_DUP_CLIPPING_MAG_SUM"

CLIPPING_ANNOTATION_FP="$({
    sha256sum \
        "$FALSE_DUP_ALL_PAIR_SUM" "$FALSE_DUP_ALL_BOUNDARY_SUM" \
        "$FALSE_DUP_ALL_CONTIG_SUM" "$CLIPPING_ANNOTATOR"
    if [[ "$FALSE_DUP_CLIPPING_INPUT" != "NA" ]]; then
        sha256sum "$FALSE_DUP_CLIPPING_INPUT"
    else
        printf '%s\n' 'clipping_input=NA'
    fi
    printf '%s\n' \
        "annotation_version=${CLIPPING_ANNOTATION_VERSION}" \
        'match_rule=1based_boundary_minus_1_or_boundary' \
        'interpretation=co_occurrence_only_no_change_to_false_duplication_call'
} | sha256sum | awk '{print $1}')"

TMP_CLIPPING_PROV="${FALSE_DUP_CLIPPING_PROVENANCE}.tmp.$$"
{
    printf 'key\tvalue\n'
    printf 'fingerprint\t%s\n' "$CLIPPING_ANNOTATION_FP"
    printf 'annotation_version\t%s\n' "$CLIPPING_ANNOTATION_VERSION"
    printf 'clipping_data_available\t%s\n' "$([[ "$FALSE_DUP_CLIPPING_INPUT" == "NA" ]] && printf N || printf Y)"
    printf 'clipping_input\t%s\n' "$FALSE_DUP_CLIPPING_INPUT"
    printf 'coordinate_match_rule\t%s\n' '1-based clipping position equals boundary-1 or boundary'
    printf 'biological_interpretation\t%s\n' 'validated clipping is a co-occurring signal from the same HiFi mapping and does not alter the false-duplication call'
} > "$TMP_CLIPPING_PROV"
mv -f "$TMP_CLIPPING_PROV" "$FALSE_DUP_CLIPPING_PROVENANCE"

for output in \
    "$FALSE_DUP_CANDIDATE_PAF" "$FALSE_DUP_CANDIDATE_PROVENANCE" \
    "$FALSE_DUP_ALL_PAIR_SUM" "$FALSE_DUP_ALL_BOUNDARY_SUM" "$FALSE_DUP_ALL_CONTIG_SUM" "$FALSE_DUP_ASSEMBLER_SUM" \
    "$FALSE_DUP_PAIR_SUM" "$FALSE_DUP_BOUNDARY_SUM" "$FALSE_DUP_CONTIG_SUM" "$FALSE_DUP_MAG_SUM" \
    "$FALSE_DUP_IDN_EVENTS" "$FALSE_DUP_IDN_SUMMARY" \
    "$FALSE_DUP_CANDIDATE_PREPROCESSING" "$FALSE_DUP_PROVENANCE" \
    "$FALSE_DUP_CLIPPING_BOUNDARY_SUM" "$FALSE_DUP_CLIPPING_PAIR_SUM" \
    "$FALSE_DUP_CLIPPING_CONTIG_SUM" "$FALSE_DUP_CLIPPING_MAG_SUM" \
    "$FALSE_DUP_CLIPPING_PROVENANCE"; do
    [[ -s "$output" ]] || { echo "ERROR: expected output missing: $output"; exit 1; }
done

echo "Candidate self-alignment    : $FALSE_DUP_CANDIDATE_PAF"
echo "All-contig pair output      : $FALSE_DUP_ALL_PAIR_SUM"
echo "All-contig boundary output  : $FALSE_DUP_ALL_BOUNDARY_SUM"
echo "All-contig summary          : $FALSE_DUP_ALL_CONTIG_SUM"
echo "Assembler summary           : $FALSE_DUP_ASSEMBLER_SUM"
echo "MAG pair output             : $FALSE_DUP_PAIR_SUM"
echo "MAG boundary output         : $FALSE_DUP_BOUNDARY_SUM"
echo "MAG-contig output           : $FALSE_DUP_CONTIG_SUM"
echo "MAG output                  : $FALSE_DUP_MAG_SUM"
echo "I/D/N event QC              : $FALSE_DUP_IDN_EVENTS"
echo "I/D/N summary QC            : $FALSE_DUP_IDN_SUMMARY"
echo "Candidate preprocessing     : $FALSE_DUP_CANDIDATE_PREPROCESSING"
echo "Clipping boundary overlap   : $FALSE_DUP_CLIPPING_BOUNDARY_SUM"
echo "Clipping pair summary       : $FALSE_DUP_CLIPPING_PAIR_SUM"
echo "Clipping contig summary     : $FALSE_DUP_CLIPPING_CONTIG_SUM"
echo "Clipping MAG summary        : $FALSE_DUP_CLIPPING_MAG_SUM"
echo "[08_false_duplication] Done: $(date)"
