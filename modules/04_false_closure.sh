#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_ID="04_false_closure"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

FALSE_CLOSURE_SCRIPT="$SCRIPT_DIR/../scripts/false_closure.py"
CLIP_ALL_FINAL="${CLIP_ALL_FINAL:-${CLIPPING_DIR}/.cache/clipping.site.tsv}"
FALSE_CLOSURE_LARGE_INSERTION_MIN_BP="${FALSE_CLOSURE_LARGE_INSERTION_MIN_BP:-100}"
FALSE_CLOSURE_LARGE_INSERTION_MIN_QNAMES="${FALSE_CLOSURE_LARGE_INSERTION_MIN_QNAMES:-2}"
FALSE_CLOSURE_LARGE_INSERTION_TSV="${FALSE_CLOSURE_LARGE_INSERTION_TSV:-$FALSE_CLOSURE_RESULT_DIR/.cache/false_closure.large_insertion_sensitivity.tsv}"
FALSE_CLOSURE_CONTIG_CIRCULARITY_TSV="${FALSE_CLOSURE_CONTIG_CIRCULARITY_TSV:-$FALSE_CLOSURE_RESULT_DIR/.cache/contig_circularity.resolved.tsv}"
FALSE_CLOSURE_FLYE_ASSEMBLY_INFO="${FALSE_CLOSURE_FLYE_ASSEMBLY_INFO:-}"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_files() {
    local path
    for path in "$@"; do
        [[ -s "$path" ]] || fail "required false-closure input missing: $path"
    done
}

require_tools() {
    local tool
    for tool in "$@"; do
        command -v "$tool" >/dev/null 2>&1 || fail "missing tool: $tool"
    done
}

validate_source_mapping_manifest() {
    python3 - "$ASSEMBLER_MAPPING_MANIFEST" <<'PY_VALIDATE_MANIFEST'
import csv
import os
import sys
path = sys.argv[1]
with open(path, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'Assembler', 'Assembly_FASTA', 'HiFi_all_BAM', 'HiFi_primary_supplementary_BAM'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'ERROR: assembler mapping manifest lacks columns: {sorted(missing)}')
    rows = list(reader)
if not rows:
    raise SystemExit('ERROR: assembler mapping manifest contains no data rows')
for row in rows:
    assembler = (row.get('Assembler') or '').strip()
    assembly = (row.get('Assembly_FASTA') or '').strip()
    bam = (row.get('HiFi_all_BAM') or '').strip()
    ps_bam = (row.get('HiFi_primary_supplementary_BAM') or '').strip()
    for p, label in ((assembly, 'assembly'), (assembly + '.fai', 'FAI'), (bam, 'HiFi_all_BAM'), (ps_bam, 'HiFi_primary_supplementary_BAM')):
        if not p or not os.path.isfile(p) or os.path.getsize(p) == 0:
            raise SystemExit(f'ERROR: {assembler} {label} missing: {p}')
PY_VALIDATE_MANIFEST

    while IFS=$'\t' read -r assembler assembly hifi_all hifi_primary_supplementary; do
        samtools quickcheck -v "$hifi_all" >/dev/null 2>&1 \
            || fail "HiFi all-mapped BAM failed quickcheck: $hifi_all"
        samtools quickcheck -v "$hifi_primary_supplementary" >/dev/null 2>&1 \
            || fail "HiFi primary+supplementary BAM failed quickcheck: $hifi_primary_supplementary"
        qc_bam_matches_fai "$hifi_all" "${assembly}.fai" \
            || fail "HiFi all-mapped BAM/reference mismatch: $assembler"
        qc_bam_matches_fai "$hifi_primary_supplementary" "${assembly}.fai" \
            || fail "HiFi primary+supplementary BAM/reference mismatch: $assembler"
    done < <(qc_manifest_columns "$ASSEMBLER_MAPPING_MANIFEST" Assembler Assembly_FASTA HiFi_all_BAM HiFi_primary_supplementary_BAM)
}
[[ "$FALSE_CLOSURE_LARGE_INSERTION_MIN_BP" =~ ^[1-9][0-9]*$ ]] \
    || fail "FALSE_CLOSURE_LARGE_INSERTION_MIN_BP must be a positive integer"
[[ "$FALSE_CLOSURE_LARGE_INSERTION_MIN_QNAMES" =~ ^[1-9][0-9]*$ ]] \
    || fail "FALSE_CLOSURE_LARGE_INSERTION_MIN_QNAMES must be a positive integer"
[[ "$PREMATURE_MAX_LENGTH_BP" =~ ^[1-9][0-9]*$ ]] \
    || fail "PREMATURE_MAX_LENGTH_BP must be a positive integer"
[[ "$PREMATURE_MIN_RIBOSOMAL_FAMILIES" =~ ^[1-9][0-9]*$ ]] \
    || fail "PREMATURE_MIN_RIBOSOMAL_FAMILIES must be a positive integer"

[[ -s "$SCG_TSV" ]] \
    || fail "SCG table missing: $SCG_TSV. Run MAGverify prep-scg -a DIR -o DIR (same -o) first"

require_files \
    "$FALSE_CLOSURE_SCRIPT" \
    "$ASSEMBLY_FA" "$MAG_CONTIG_MAP" "$CIRCULARITY_TSV" "$SCG_TSV" \
    "$ASSEMBLER_MAPPING_MANIFEST" "$CLIP_ALL_FINAL" "$CLIP_FINAL"
FLYE_INFO_ARGS=()
if [[ -n "$FALSE_CLOSURE_FLYE_ASSEMBLY_INFO" ]]; then
    require_files "$FALSE_CLOSURE_FLYE_ASSEMBLY_INFO"
    FLYE_INFO_ARGS=(--flye-assembly-info "$FALSE_CLOSURE_FLYE_ASSEMBLY_INFO")
fi
mkdir -p "$OUT/00_log" "$FALSE_CLOSURE_RESULT_DIR/.cache"
for legacy in false_closure.MAG_contig.tsv contig_circularity.resolved.tsv \
    assembly_contig_circularity.resolved.tsv false_closure.cross_assembler_sensitivity.tsv \
    false_closure.MAG.audit.tsv run_metadata.txt false_closure_candidate_reads.tsv \
    false_closure.tsv false_closure.breakpoints.tsv false_closure.large_insertion_sensitivity.tsv \
    false_closure.breakpoint.tsv; do
    rm -f "$FALSE_CLOSURE_RESULT_DIR/$legacy"
done
rm -rf "$FALSE_CLOSURE_RESULT_DIR/validation"
require_tools python3 samtools

LOG="$OUT/00_log/04_false_closure.log"
exec > >(tee -a "$LOG") 2>&1

echo "[$SCRIPT_ID] Started: $(date)"

samtools faidx "$ASSEMBLY_FA"
qc_validate_current_clipping_result \
    || fail "clipping result/cache is stale; rerun 03_clipping.sh"
validate_source_mapping_manifest

python3 - "$CLIP_ALL_FINAL" "$FALSE_CLOSURE_MIN_LINKED_READS_PER_SIDE" <<'PY_VALIDATE'
import csv
import sys
path, linked_text = sys.argv[1:]
linked = int(linked_text)
if linked < 1:
    raise SystemExit('ERROR: FALSE_CLOSURE_MIN_LINKED_READS_PER_SIDE must be at least 1')
with open(path, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    fields = set(reader.fieldnames or [])
    required = {
        'Assembler', 'region_type', 'reference_boundary',
        'Clipping_analyzable', 'Clipping', 'validated_clipping',
        'n_unique_all_clipped_qnames',
        'n_unique_spanning_with_I_qnames', 'spanning_I_length_bp_min',
        'spanning_I_length_bp_median', 'spanning_I_length_bp_max',
    }
    missing = required - fields
    if missing:
        raise SystemExit(f'ERROR: all-region clipping result lacks columns: {sorted(missing)}')
    spanning_column = (
        'n_unique_exact_spanning_qnames'
        if 'n_unique_exact_spanning_qnames' in fields
        else ('n_unique_clean_spanning_qnames' if 'n_unique_clean_spanning_qnames' in fields else None)
    )
    if spanning_column is None:
        raise SystemExit(
            'ERROR: all-region clipping result lacks a spanning-QNAME column; '
            'expected n_unique_exact_spanning_qnames'
        )
    bad_region = []
    for row in reader:
        region = (row.get('region_type') or '').strip()
        if region not in {'INTERNAL', 'TERMINAL'}:
            bad_region.append(region)
            break
        analyzable = (row.get('Clipping_analyzable') or '').strip().upper()
        clipping = (row.get('Clipping') or '').strip().upper()
        validated = (row.get('validated_clipping') or '').strip().upper()
        if analyzable not in {'Y', 'N'}:
            raise SystemExit(f'ERROR: invalid Clipping_analyzable value: {analyzable!r}')
        if analyzable == 'N':
            if clipping != 'NA':
                raise SystemExit('ERROR: Clipping_analyzable=N requires Clipping=NA')
            continue
        if clipping not in {'Y', 'N'}:
            raise SystemExit(f'ERROR: Clipping_analyzable=Y requires Clipping=Y/N, observed {clipping!r}')
        if clipping == 'Y' and validated != 'Y':
            raise SystemExit('ERROR: Clipping=Y row is not validated_clipping=Y')
        if clipping == 'Y' and int(row.get(spanning_column) or 0) != 0:
            raise SystemExit('ERROR: Clipping=Y row has non-zero spanning QNAMEs')
if bad_region:
    raise SystemExit(f'ERROR: invalid region_type in all-region clipping result: {bad_region[0]!r}')
PY_VALIDATE

echo "Candidate clipping table: $CLIP_ALL_FINAL"
echo "Candidate admission: every final HiFi Clipping=Y site across all assembly contigs"
echo "Clipping regions: TERMINAL + INTERNAL; region_type is annotation only"
echo "Exact physical reference boundary: excluded as BOUNDARY_UNRESOLVED"
echo "HiFi mapping context: per-assembler original assembly"
echo "HiFi structural alignment input: HiFi_primary_supplementary_BAM (primary + supplementary only)"
echo "Structural connection: split alignment involving >=1 supplementary segment"
echo "Secondary alignments in structural evidence: excluded"
echo "Read counting: unique QNAME"
echo "Minimum cross-contig connected QNAMEs per side: $FALSE_CLOSURE_MIN_LINKED_READS_PER_SIDE"
echo "Same target contig on both sides required: N"
echo "False-closure circularity gate: assembler contig-level circularity; MAG-level circularity is single-contig fallback/audit only"
echo "Assembler circularity metadata: metaMDBG FASTA circular=yes/no; hifiasm-meta contig suffix c/l; metaFlye assembly_info.txt when available"
echo "Contig-level scope: all assembler-reported circular contigs, independent of MAG membership"
echo "MAG-level scope: final MAG subset summarized after contig-level classification"
echo "Premature-circularization candidate MAG: single-contig MAG only"
echo "Large-insertion sensitivity: >=${FALSE_CLOSURE_LARGE_INSERTION_MIN_QNAMES} exact-spanning HiFi QNAMEs with per-QNAME insertion length >=${FALSE_CLOSURE_LARGE_INSERTION_MIN_BP} bp"
echo "Large-insertion sensitivity changes final HiFi FALSE_CLOSURE call: N"
echo "Additional HiFi MAPQ/depth/identity/anchor filter: none"
echo "Paper premature-circularization proxy: assembler circular AND length < ${PREMATURE_MAX_LENGTH_BP} bp AND >=${PREMATURE_MIN_RIBOSOMAL_FAMILIES} distinct ribosomal-protein families"
echo "Paper proxy changes final HiFi FALSE_CLOSURE call: N"

python3 "$FALSE_CLOSURE_SCRIPT" resolve \
    --mag-contig-map "$MAG_CONTIG_MAP" \
    --mapping-manifest "$ASSEMBLER_MAPPING_MANIFEST" \
    --mag-circularity "$CIRCULARITY_TSV" \
    "${FLYE_INFO_ARGS[@]}" \
    --output "$FALSE_CLOSURE_CONTIG_CIRCULARITY_TSV"

python3 "$FALSE_CLOSURE_SCRIPT" detect \
    --assemblies-fai "${ASSEMBLY_FA}.fai" \
    --mapping-manifest "$ASSEMBLER_MAPPING_MANIFEST" \
    --clipping-all-final-tsv "$CLIP_ALL_FINAL" \
    --contig-circularity "$FALSE_CLOSURE_CONTIG_CIRCULARITY_TSV" \
    --mag-contig-map "$MAG_CONTIG_MAP" \
    --scg-tsv "$SCG_TSV" \
    --premature-max-length-bp "$PREMATURE_MAX_LENGTH_BP" \
    --premature-min-ribosomal-families "$PREMATURE_MIN_RIBOSOMAL_FAMILIES" \
    --min-linked-reads-per-side "$FALSE_CLOSURE_MIN_LINKED_READS_PER_SIDE" \
    --large-insertion-min-bp "$FALSE_CLOSURE_LARGE_INSERTION_MIN_BP" \
    --large-insertion-min-qnames "$FALSE_CLOSURE_LARGE_INSERTION_MIN_QNAMES" \
    --large-insertion-details "$FALSE_CLOSURE_LARGE_INSERTION_TSV" \
    --candidate-read-output "$FALSE_CLOSURE_CANDIDATE_READS" \
    --output "$FALSE_CLOSURE_EVIDENCE_TSV" \
    --details "$FALSE_CLOSURE_BREAKPOINTS_TSV"

[[ -s "$FALSE_CLOSURE_CONTIG_CIRCULARITY_TSV" && -s "$FALSE_CLOSURE_EVIDENCE_TSV" && -s "$FALSE_CLOSURE_BREAKPOINTS_TSV" && -s "$FALSE_CLOSURE_LARGE_INSERTION_TSV" ]] \
    || fail "false-closure detection output was not created"

python3 "$FALSE_CLOSURE_SCRIPT" mag \
    --false-closure-tsv "$FALSE_CLOSURE_EVIDENCE_TSV" \
    --false-final "$FALSE_CLOSURE_FINAL" \
    --false-audit "$FALSE_CLOSURE_AUDIT" \
    --metadata "$FALSE_CLOSURE_RUN_METADATA"

python3 "$FALSE_CLOSURE_SCRIPT" contig \
    --mag-contig-map "$MAG_CONTIG_MAP" \
    --assemblies-fai "${ASSEMBLY_FA}.fai" \
    --contig-circularity "$FALSE_CLOSURE_CONTIG_CIRCULARITY_TSV" \
    --clipping-all-final-tsv "$CLIP_ALL_FINAL" \
    --details-tsv "$FALSE_CLOSURE_BREAKPOINTS_TSV" \
    --scg-tsv "$SCG_TSV" \
    --premature-max-length-bp "$PREMATURE_MAX_LENGTH_BP" \
    --premature-min-ribosomal-families "$PREMATURE_MIN_RIBOSOMAL_FAMILIES" \
    --output "$FALSE_CLOSURE_CONTIG_FINAL"

cat >> "$FALSE_CLOSURE_RUN_METADATA" <<EOF_METADATA
script_id\t$SCRIPT_ID
candidate_table\t$CLIP_ALL_FINAL
clipping_regions\tTERMINAL_AND_INTERNAL_ANALYZED_IDENTICALLY
candidate_rule\tall_final_HiFi_Clipping_Y_sites_across_all_assembler_contigs
exact_reference_boundary\tBOUNDARY_UNRESOLVED_EXCLUDED
reference\tall_contigs_of_the_producing_assembler
alignment_population_main\tHiFi_primary_supplementary_BAM
structural_connection_rule\tcross_contig_split_alignment_involving_supplementary
secondary_alignments_structural_evidence\texcluded
read_counting\tunique_QNAME
minimum_connected_QNAMEs_per_side\t$FALSE_CLOSURE_MIN_LINKED_READS_PER_SIDE
false_closure_analyzability\tpositive_structural_evidence_is_definitive; otherwise contig_circularity_Y_AND_no_unresolved_clipping_candidate_on_that_contig
target_contig_matching_between_sides\tnot_required
contig_level_scope\tall_assembler_reported_circular_source_contigs
MAG_level_scope\tfinal_MAG_subset_after_contig_level_call
premature_circularization_MAG_rule\tsingle_contig_MAG_only
large_insertion_sensitivity_min_bp\t$FALSE_CLOSURE_LARGE_INSERTION_MIN_BP
large_insertion_sensitivity_min_QNAMEs\t$FALSE_CLOSURE_LARGE_INSERTION_MIN_QNAMES
large_insertion_sensitivity_rule\texact_spanning_HiFi_QNAMEs_with_candidate_local_read_side_insertion_MIN_length_threshold
large_insertion_sensitivity_changes_false_closure_call\tN
contig_circularity_table\t$FALSE_CLOSURE_CONTIG_CIRCULARITY_TSV
additional_depth_filter\tnone
additional_clipping_fraction_filter\tnone
identity_filter\tnone
MAPQ_filter\tnone
anchor_filter\tnone
EOF_METADATA

echo "Outputs:"
echo "  $FALSE_CLOSURE_FINAL"
echo "  $FALSE_CLOSURE_CONTIG_FINAL"
echo "[$SCRIPT_ID] Finished: $(date)"
