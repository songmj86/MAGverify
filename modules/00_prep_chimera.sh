#!/usr/bin/env bash
# Build an MMseqs2 protein taxonomy DB using manually downloaded GTDB R232 files.
# Based on the MMseqs2 GTDB database workflow (tar2db + createtaxdb).
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo 'Internal helper: use MAGverify prep-chimera --gtdb-dir DIR --output-dir DIR [--database-prefix PREFIX] --threads N' >&2
    exit 2
fi

PREP=$(mkdir -p "$1" && cd "$1" && pwd)
THREADS="$2"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'threads must be a positive integer' >&2; exit 2; }
[[ -d "$3" ]] || { echo "GTDB input directory not found: $3" >&2; exit 2; }
GTDB_DIR=$(cd "$3" && pwd)
DB_DIR=$(dirname "$4")
DB_NAME=$(basename "$4")
[[ -n "$DB_NAME" && "$DB_NAME" != '.' && "$DB_NAME" != '..' ]] || { echo 'Invalid database prefix' >&2; exit 2; }
mkdir -p "$DB_DIR"
DB="$(cd "$DB_DIR" && pwd)/$DB_NAME"
for executable in mmseqs awk sed; do
    command -v "$executable" >/dev/null 2>&1 || { echo "Missing command: $executable" >&2; exit 1; }
done

TMP="$PREP/mmseqs_tmp"
mkdir -p "$TMP/taxonomy"

# These files must be downloaded by the user before running this command.
PROTEINS="$GTDB_DIR/gtdb_proteins_aa_reps_r232.tar.gz"
BAC="$GTDB_DIR/bac120_taxonomy_r232.tsv"
ARC="$GTDB_DIR/ar53_taxonomy_r232.tsv"
for file in "$PROTEINS" "$BAC" "$ARC"; do
    [[ -s "$file" ]] || { echo "GTDB R232 input missing or empty: $file" >&2; exit 1; }
done
printf 'Using GTDB R232 inputs from: %s\n' "$GTDB_DIR"
printf 'MMseqs2 taxonomy database prefix: %s\n' "$DB"

# Each nested *_protein.faa.gz is ingested and its genome accession
# becomes the identifier needed for mapping GTDB taxonomy to sequences.
mmseqs tar2db "$PROTEINS" "$TMP/tardb" --tar-include 'faa.gz$' --threads "$THREADS"
sed 's|_protein\.faa\.gz||g' "$TMP/tardb.lookup" > "$TMP/tardb.lookup.tmp"
mv "$TMP/tardb.lookup.tmp" "$TMP/tardb.lookup"
mmseqs createdb "$TMP/tardb" "$DB"

# Reproduce the MMseqs2 GTDB taxonomy workflow using the fixed R232 TSVs.
# Do not use the NCBI taxdump: GTDB rank assignments are distinct.
awk -v taxdir="$TMP/taxonomy" '
BEGIN {
    FS = "[\t;]"
    rank["c"] = "class"; rank["d"] = "superkingdom"
    rank["f"] = "family"; rank["g"] = "genus"
    rank["o"] = "order"; rank["p"] = "phylum"; rank["s"] = "species"
    taxCnt = 1
    print "1\t|\t1\t|\tno rank\t|\t-\t|" > taxdir"/nodes.dmp"
    print "1\t|\troot\t|\t-\t|\tscientific name\t|" > taxdir"/names.dmp"
}
{
    accession = $1
    prevTaxon = 1
    if (accession == "" || NF < 2) {
        print "Invalid taxonomy row: " $0 > "/dev/stderr"
        exit 2
    }
    for (i = 2; i <= NF; i++) {
        token = $i
        sub(/\r$/, "", token)
        r = substr(token, 1, 1)
        if (!(r in rank) || substr(token, 2, 2) != "__") {
            print "Invalid GTDB rank token: " token > "/dev/stderr"
            exit 2
        }
        if (token == r "__") continue
        # Taxon IDs are lineage-specific, not keyed by the rank label alone.
        # Otherwise unclassified ranks (e.g. g__) merge unrelated clades.
        key = prevTaxon SUBSEP token
        if (key in ids) {
            prevTaxon = ids[key]
        } else {
            taxCnt++
            ids[key] = taxCnt
            name = substr(token, 4)
            gsub(/_/, " ", name)
            printf("%s\t|\t%s\t|\t%s\t|\t-\t|\n", taxCnt, prevTaxon, rank[r]) > taxdir"/nodes.dmp"
            printf("%s\t|\t%s\t|\t-\t|\tscientific name\t|\n", taxCnt, name) > taxdir"/names.dmp"
            prevTaxon = taxCnt
        }
    }
    printf("%s\t%s\n", accession, prevTaxon) > taxdir"/mapping_genomes"
}' "$BAC" "$ARC"
: > "$TMP/taxonomy/merged.dmp"
: > "$TMP/taxonomy/delnodes.dmp"
mmseqs createtaxdb "$DB" "$TMP/taxdb" \
  --ncbi-tax-dump "$TMP/taxonomy" \
  --tax-mapping-file "$TMP/taxonomy/mapping_genomes" \
  --tax-mapping-mode 1 \
  --threads "$THREADS"

[[ -s "${DB}.dbtype" && -s "${DB}_mapping" ]] || {
    echo 'Incomplete MMseqs2 GTDB taxonomy database; check the command log.' >&2
    exit 1
}
printf 'GTDB release R232 (232.0)\n' > "${DB}.version"
printf 'GTDB R232 MMseqs2 database ready: %s\n' "$DB"
