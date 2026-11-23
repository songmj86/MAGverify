#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_ID="10_summary"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

GTDB_CHIMERA_REPORT_CONTIG="${GTDB_CHIMERA_REPORT_CONTIG:-$CHIMERA_ROOT_DIR/.cache/chimera.contig.tsv}"

mkdir -p "$OUT/00_log" "$(dirname "$SUMMARY_CONTIG")" "$(dirname "$SUMMARY_MAG")"
LOG="$OUT/00_log/10_summary.log"
exec > >(tee -a "$LOG") 2>&1

for tool in samtools python3; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing tool: $tool"; exit 1; }
done

[[ "$RUN_CLIPPING" != "Y" || -s "$CLIP_FINAL" ]] \
    || { echo "ERROR: enabled clipping result is missing: $CLIP_FINAL"; exit 1; }
[[ "$RUN_TAXONOMIC_CHIMERA" != "Y" || -s "$GTDB_CHIMERA_REPORT_CONTIG" ]] \
    || { echo "ERROR: enabled clipping-defined GTDB taxonomic-chimera result is missing: $GTDB_CHIMERA_REPORT_CONTIG"; exit 1; }
[[ "$RUN_ASSEMBLY_VARIANTS" != "Y" || -s "$ALLELE_CONTIG" ]] \
    || { echo "ERROR: enabled variant result is missing: $ALLELE_CONTIG"; exit 1; }
[[ "$RUN_PHANTOM" != "Y" || -s "$PHANTOM_FINAL" ]] \
    || { echo "ERROR: enabled phantom result is missing: $PHANTOM_FINAL"; exit 1; }
[[ "$RUN_REPEAT" != "Y" || -s "$REPEAT_SUM" ]] \
    || { echo "ERROR: enabled repeat result is missing: $REPEAT_SUM"; exit 1; }
[[ "$RUN_REPEAT" != "Y" || -s "$REPEAT_MAG_SUM" ]] \
    || { echo "ERROR: enabled MAG repeat result is missing: $REPEAT_MAG_SUM"; exit 1; }
[[ "$RUN_FALSE_DUPLICATION" != "Y" || -s "$FALSE_DUP_CONTIG_SUM" ]] \
    || { echo "ERROR: enabled false-duplication result is missing: $FALSE_DUP_CONTIG_SUM"; exit 1; }
[[ "$RUN_FALSE_DUPLICATION" != "Y" || -s "$FALSE_DUP_CLIPPING_CONTIG_SUM" ]] \
    || { echo "ERROR: false-duplication clipping contig summary is missing: $FALSE_DUP_CLIPPING_CONTIG_SUM"; exit 1; }
[[ "$RUN_FALSE_DUPLICATION" != "Y" || -s "$FALSE_DUP_CLIPPING_MAG_SUM" ]] \
    || { echo "ERROR: false-duplication clipping MAG summary is missing: $FALSE_DUP_CLIPPING_MAG_SUM"; exit 1; }
if [[ "$RUN_CIRCULARITY" == "Y" ]]; then
    [[ -s "$CIRCULARITY_ONLY_FINAL" && -s "$CIRCULARITY_ONLY_AUDIT" ]] || {
        echo "ERROR: enabled circularity-only final/audit result is missing"
        exit 1
    }
fi
if [[ "$RUN_FALSE_CLOSURE" == "Y" ]]; then
    [[ -s "$FALSE_CLOSURE_FINAL" && -s "$FALSE_CLOSURE_AUDIT" ]] || {
        echo "ERROR: enabled false-closure final/audit result is missing"
        exit 1
    }
fi
samtools faidx "$ASSEMBLY_FA"

TMP_CONTIG="${SUMMARY_CONTIG}.tmp.$$"
TMP_MAG="${SUMMARY_MAG}.tmp.$$"
trap 'rm -f "$TMP_CONTIG" "$TMP_MAG"' EXIT

echo "[10_summary] Started: $(date)"
echo "SCRIPT_ID: $SCRIPT_ID"

python3 "$SCRIPT_DIR/../scripts/summary_tables.py" build \
    "$ASSEMBLY_FA.fai" "$CLIP_FINAL" "$GTDB_CHIMERA_REPORT_CONTIG" "$ALLELE_CONTIG" \
    "$PHANTOM_FINAL" "$REPEAT_SUM" "$REPEAT_MAG_SUM" "$FALSE_DUP_CONTIG_SUM" \
    "$FALSE_DUP_CLIPPING_CONTIG_SUM" "$FALSE_DUP_CLIPPING_MAG_SUM" "$MAG_CONTIG_MAP" \
    "$CIRCULARITY_ONLY_FINAL" "$CIRCULARITY_ONLY_AUDIT" \
    "$FALSE_CLOSURE_FINAL" "$FALSE_CLOSURE_AUDIT" \
    "$TMP_CONTIG" "$TMP_MAG" \
    "$RUN_CLIPPING" "$RUN_TAXONOMIC_CHIMERA" "$RUN_ASSEMBLY_VARIANTS" "$RUN_PHANTOM" \
    "$RUN_REPEAT" "$RUN_FALSE_DUPLICATION" "$RUN_CIRCULARITY" \
    "$RUN_FALSE_CLOSURE"

python3 "$SCRIPT_DIR/../scripts/summary_tables.py" arrange "$TMP_MAG"
python3 "$SCRIPT_DIR/../scripts/summary_tables.py" validate "$TMP_CONTIG" "$TMP_MAG"
mv -f "$TMP_CONTIG" "$SUMMARY_CONTIG"
mv -f "$TMP_MAG" "$SUMMARY_MAG"
rm -f "$(dirname "$SUMMARY_MAG")/integrated_QC_summary.total_contigs.tsv" \
    "$(dirname "$SUMMARY_MAG")/integrated_QC_flags.by_MAG.tsv" \
    "$(dirname "$SUMMARY_MAG")/complete_genome_validity_QC.by_MAG.tsv"

echo "Outputs:"
echo "  $SUMMARY_MAG"
echo "  $SUMMARY_CONTIG"
echo "[10_summary] Finished: $(date)"
