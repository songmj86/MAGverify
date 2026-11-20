#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../modules" && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/MAGs" "$TMP/source" "$TMP/01_mapping"
printf '>contig1\nACGT\n' > "$TMP/combined.fa"
cp "$TMP/combined.fa" "$TMP/MAGs/MAG1.fa"
cp "$TMP/combined.fa" "$TMP/source/assembler.fa"
printf 'dummy\n' > "$TMP/01_mapping/all.bam"
mkdir -p "$TMP/01_mapping/.cache"
printf 'fingerprint\n' > "$TMP/01_mapping/.cache/all.bam.inputs.sha256"
manifest="$TMP/01_mapping/mapping_manifest.tsv"
printf 'Assembler\tAssembly_FASTA\tHiFi_all_BAM\tHiFi_primary_BAM\tHiFi_primary_supplementary_BAM\tOutput_directory\n' > "$manifest"
printf 'assembler\t%s\t%s\t%s\t%s\t%s\n' \
    "$TMP/source/assembler.fa" "$TMP/01_mapping/all.bam" \
    "$TMP/01_mapping/primary.bam" "$TMP/01_mapping/primary_supplementary.bam" \
    "$TMP/01_mapping" >> "$manifest"

export MAG_DIR="$TMP/MAGs"
export ASSEMBLER_MAPPING_MANIFEST="$manifest"
python3 "$ROOT/../scripts/clipping_check_clipping_inputs.py" --assembly "$TMP/combined.fa" --mag-dir "$MAG_DIR" --manifest "$ASSEMBLER_MAPPING_MANIFEST" > "$TMP/success.log"
grep -q '1 MAG FASTAs; 1 assembler' "$TMP/success.log"

# Missing prior mapping must be rejected before starting any downstream analysis.
rm "$TMP/01_mapping/all.bam"
if bash "$ROOT/03_clipping.sh" "$TMP/combined.fa" "$TMP/03_clipping" > "$TMP/missing.log" 2>&1; then
    echo 'ERROR: clipping accepted missing BAM' >&2
    exit 1
fi
grep -q 'HiFi_all_BAM for assembler' "$TMP/missing.log"

# Missing MAG_DIR must be reported directly, not as a samtools or anvi-o error.
if MAG_DIR="$TMP/not_mags" bash "$ROOT/03_clipping.sh" "$TMP/combined.fa" "$TMP/03_clipping" > "$TMP/no_mags.log" 2>&1; then
    echo 'ERROR: clipping accepted missing MAG_DIR' >&2
    exit 1
fi
grep -q 'MAG_DIR missing' "$TMP/no_mags.log"

# A bad manifest is rejected even if the file exists.
printf 'bad\theader\nfoo\tbar\n' > "$manifest"
if python3 "$ROOT/../scripts/clipping_check_clipping_inputs.py" --assembly "$TMP/combined.fa" --mag-dir "$MAG_DIR" --manifest "$ASSEMBLER_MAPPING_MANIFEST" > "$TMP/header.log" 2>&1; then
    echo 'ERROR: clipping accepted invalid manifest header' >&2
    exit 1
fi
grep -q 'Invalid mapping manifest header' "$TMP/header.log"

echo 'Clipping input validation tests passed.'
