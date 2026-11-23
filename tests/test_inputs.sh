#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../modules" && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/mags"

printf 'contig_a\t4\t0\t4\t5\ncontig_b\t6\t0\t6\t7\n' > "$TMP/assembly.fasta.fai"
printf '>contig_a\nACGT\n' > "$TMP/mags/MAG01.fa"
printf '>contig_b\nGATTAC\n' > "$TMP/mags/MAG02.fa"

python3 "$ROOT/../scripts/prep_build_mag_contig_map.py" \
    --mag-dir "$TMP/mags" --assembly-fai "$TMP/assembly.fasta.fai" \
    --output "$TMP/map.tsv" > /dev/null
cmp - "$TMP/map.tsv" <<'EXPECTED'
MAG_ID	Contig
MAG01	contig_a
MAG02	contig_b
EXPECTED

python3 "$ROOT/../scripts/prep_build_mag_contig_map.py" \
    --mag-dir "$TMP/mags" --assembly-fai "$TMP/assembly.fasta.fai" \
    --output "$TMP/map.tsv" > /dev/null

printf '>contig_a\nACGTA\n' > "$TMP/mags/MAG01.fa"
if python3 "$ROOT/../scripts/prep_build_mag_contig_map.py" \
    --mag-dir "$TMP/mags" --assembly-fai "$TMP/assembly.fasta.fai" \
    --output "$TMP/invalid.tsv" > /dev/null 2>&1; then
    echo 'ERROR: inconsistent sequence length unexpectedly accepted' >&2
    exit 1
fi
printf 'MAG-to-contig input tests passed.\n'
