#!/usr/bin/env bash
# MAGverify preparation commands.
# Run through the CLI: MAGverify prep-chimera | prep-scg | prep-circularity | prep-taxonomy
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'HELP'
Usage: bash 00_prep.sh SUBCOMMAND [ARGS ...]

Subcommands:
  chimera     TMP_DIR THREADS GTDB_DIR DB_PREFIX
  scg         ASSEMBLY_FASTA OUTPUT_DIR [THREADS]
  circularity MAG_DIR OUTPUT_DIR ASSEMBLER_FILE [ASSEMBLER_FILE ...]
  taxonomy    GTDBTK_DIR MAG_DIR OUTPUT_TSV

Run through the MAGverify CLI instead:
  MAGverify prep-chimera -d DIR -o DIR -t N
  MAGverify prep-scg -a DIR -o DIR -t N
  MAGverify prep-circularity -c FILE [-c FILE ...] -g DIR -o DIR
  MAGverify prep-taxonomy -d DIR -g DIR -o DIR
HELP
}

prep_chimera() {
    if [[ $# -ne 4 ]]; then
        echo 'Internal helper: use MAGverify prep-chimera -d DIR -o DIR -t N' >&2
        exit 2
    fi

    TMP=$(mkdir -p "$1" && cd "$1" && pwd)
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

    mkdir -p "$TMP/taxonomy"

    PROTEINS="$GTDB_DIR/gtdb_proteins_aa_reps_r232.tar.gz"
    BAC="$GTDB_DIR/bac120_taxonomy_r232.tsv"
    ARC="$GTDB_DIR/ar53_taxonomy_r232.tsv"
    for file in "$PROTEINS" "$BAC" "$ARC"; do
        [[ -s "$file" ]] || { echo "GTDB R232 input missing or empty: $file" >&2; exit 1; }
    done
    printf 'Using GTDB R232 inputs from: %s\n' "$GTDB_DIR"
    printf 'MMseqs2 taxonomy database prefix: %s\n' "$DB"
    printf 'Temporary MMseqs2 directory: %s\n' "$TMP"

    mmseqs tar2db "$PROTEINS" "$TMP/tardb" --tar-include 'faa.gz$' --threads "$THREADS"
    sed 's|_protein\.faa\.gz||g' "$TMP/tardb.lookup" > "$TMP/tardb.lookup.tmp"
    mv "$TMP/tardb.lookup.tmp" "$TMP/tardb.lookup"
    mmseqs createdb "$TMP/tardb" "$DB"

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
}

prep_scg() {
    if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
        cat <<'HELP'
Usage: bash 00_prep.sh scg ASSEMBLY_FASTA OUTPUT_DIR [THREADS]
Input: all-assemblies FASTA.
Additional software: pyrodigal; anvi'o (anvi-script-reformat-fasta, anvi-gen-contigs-database, anvi-run-hmms).
Genes come from the shared Pyrodigal prediction (lib_genes.sh).
Output: SCG_TSV (default OUTPUT_DIR/scg.contig.tsv); anvi'o files stay in OUTPUT_DIR.
Genes: GENE_CALLS_DIR (default OUTPUT_DIR).
HELP
        exit 0
    fi
    ASSEMBLY_FA="${1:?Usage: bash 00_prep.sh scg <ASSEMBLY_FA> <OUT_DIR> [THREADS]}"
    OUT_DIR="${2:?Usage: bash 00_prep.sh scg <ASSEMBLY_FA> <OUT_DIR> [THREADS]}"
    THREADS="${3:-${THREADS:-4}}"
    CONFIG_DIR="$SCRIPT_DIR"
    source "$SCRIPT_DIR/defaults.sh"
    if [[ -n "${ANVIO_ENV_PREFIX:-}" && -d "$ANVIO_ENV_PREFIX/bin" ]]; then
        export PATH="$ANVIO_ENV_PREFIX/bin:$PATH"
    fi

    SCG_SOURCE_REGEX="${SCG_SOURCE_REGEX:-^(Bacteria_71|Archaea_76)$}"
    GENE_CALLS_DIR="${GENE_CALLS_DIR:-$OUT_DIR}"

    mkdir -p "$OUT_DIR" "$(dirname "${SCG_TSV:-$OUT_DIR/x}")"

    for t in pyrodigal anvi-script-reformat-fasta anvi-gen-contigs-database anvi-run-hmms python3; do
        command -v "$t" >/dev/null 2>&1 || {
            echo "ERROR: $t not found. Activate the MAGverify Conda environment first."
            echo "Example: conda activate MAGverify"
            exit 1
        }
    done

    [[ -s "$ASSEMBLY_FA" ]] || {
        echo "ERROR: assembly FASTA not found or empty: $ASSEMBLY_FA"
        exit 1
    }

    REFORMATTED="$OUT_DIR/Contig.anvio_simple.fa"
    NAME_MAP="$OUT_DIR/Contig.anvio_simple.name_map.tsv"
    GENE_CALLS="$OUT_DIR/Contig.anvio_simple.gene_calls.tsv"
    DB="$OUT_DIR/contigs.db"
    OUT_TSV="${SCG_TSV:-$OUT_DIR/scg.contig.tsv}"
    LOG="$OUT_DIR/prep_scg.log"
    CACHE_META="$OUT_DIR/inputs.sha256"
    TMP_CACHE="${CACHE_META}.tmp.$$"
    TMP_OUT="${OUT_TSV}.tmp.$$"
    trap 'rm -f "$TMP_CACHE" "$TMP_OUT"' EXIT

    exec > >(tee -a "$LOG") 2>&1

    echo "[00_SCG] Started: $(date)"
    echo "Input assembly : $ASSEMBLY_FA"
    echo "Output dir     : $OUT_DIR"
    echo "Threads        : $THREADS"
    echo "HMM collections: $SCG_SOURCE_REGEX"
    echo

    bash "$SCRIPT_DIR/lib_genes.sh" "$ASSEMBLY_FA" "$GENE_CALLS_DIR" "$THREADS"
    GENES_FP="$(<"$GENE_CALLS_DIR/.cache/genes.inputs.sha256")"

    INPUT_FP="$(python3 - "$ASSEMBLY_FA" "$GENES_FP" <<'PY_FP'
import hashlib
import json
import os
import shutil
import sys

def stat(path):
    value = os.stat(path)
    return [os.path.realpath(path), value.st_size, value.st_mtime_ns]

payload = {
    'assembly': stat(sys.argv[1]),
    'genes': sys.argv[2],
    'reformat_binary': stat(shutil.which('anvi-script-reformat-fasta')),
    'contigs_db_binary': stat(shutil.which('anvi-gen-contigs-database')),
    'hmms_binary': stat(shutil.which('anvi-run-hmms')),
    'method': 'simplify_names|external_pyrodigal_genes|contigs_db|run_hmms_v3',
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_FP
    )"

    PREP_REUSED=0
    if [[ -s "$CACHE_META" && "$(<"$CACHE_META")" == "$INPUT_FP" \
          && -s "$REFORMATTED" && -s "$NAME_MAP" && -s "$GENE_CALLS" && -s "$DB" ]]; then
        PREP_REUSED=1
        echo "Fingerprint-matched anvi'o DB/HMM preparation found."
        echo "Reusing: $DB"
    else
        echo "Removing stale or incomplete anvi'o preparation outputs."
        rm -f "$REFORMATTED" "$NAME_MAP" "$GENE_CALLS" "$DB" "$CACHE_META" "$OUT_TSV"
    fi

    if (( PREP_REUSED == 0 )); then
        echo
        echo "[1/4] Reformatting FASTA for anvi'o simple deflines"
        anvi-script-reformat-fasta "$ASSEMBLY_FA" \
            --simplify-names \
            --report-file "$NAME_MAP" \
            -o "$REFORMATTED"

        python3 "$SCRIPT_DIR/../scripts/prep_genes.py" anvio \
            "$GENE_CALLS_DIR/genes.gff" "$GENE_CALLS_DIR/genes.faa" "$NAME_MAP" "$GENE_CALLS" \
            --version "$(bash "$SCRIPT_DIR/lib_genes.sh" --version)"

        echo
        echo "[2/4] Generating anvi'o contigs.db with the shared Pyrodigal genes"
        rm -f "$DB"
        anvi-gen-contigs-database \
            -f "$REFORMATTED" \
            -o "$DB" \
            -n "LR_assembly" \
            --external-gene-calls "$GENE_CALLS" \
            -T "$THREADS"

        echo
        echo "[3/4] Running anvi-run-hmms"
        anvi-run-hmms \
            -c "$DB" \
            -I Bacteria_71,Archaea_76 \
            -T "$THREADS" \
            --just-do-it

        printf '%s\n' "$INPUT_FP" > "$TMP_CACHE"
        mv -f "$TMP_CACHE" "$CACHE_META"
    else
        echo
        echo "[1/4] FASTA reformatting: reused"
        echo "[2/4] contigs.db         : reused"
        echo "[3/4] anvi-run-hmms      : reused"
    fi

    echo
    echo "[4/4] Exporting SCG HMM hits"

    python3 - "$DB" "$TMP_OUT" "$NAME_MAP" "$SCG_SOURCE_REGEX" <<'PY'
import csv
import os
import re
import sqlite3
import sys
from pathlib import Path

DB, OUT_TSV, NAME_MAP, SOURCE_REGEX = sys.argv[1:]
source_pattern = re.compile(SOURCE_REGEX)

def table_exists(cur, name):
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (name,))
    return cur.fetchone() is not None

def columns(cur, table):
    cur.execute(f'PRAGMA table_info("{table}")')
    return [r[1] for r in cur.fetchall()]

def first_existing(cols, candidates):
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    for cand in candidates:
        pat = cand.lower()
        for c in cols:
            if pat in c.lower():
                return c
    return None

def load_name_map(path):
    """Return mapping from anvi simplified contig name to original contig name.

    anvi-script-reformat-fasta --report-file output can vary slightly among
    versions. This parser tries common layouts and also handles simple 2-column
    maps. If uncertain, it keeps both directions only when one side looks like a
    simplified anvi name.
    """
    mapping = {}
    if not path:
        return mapping
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return mapping

    with p.open(errors='ignore') as f:
        rows = [line.rstrip('\n').split('\t') for line in f if line.strip()]
    if not rows:
        return mapping

    header = [x.lower().strip().replace(' ', '_').replace('-', '_') for x in rows[0]]
    has_header = any(x in header for x in ['old_name','new_name','original_name','new_defline','old_defline','contig','contig_name'])

    data = rows[1:] if has_header else rows
    h = {v:i for i,v in enumerate(header)} if has_header else {}

    def clean(x):
        x = x.strip()
        if x.startswith('>'):
            x = x[1:]
        return x.split()[0]

    if has_header:
        old_idx = None
        new_idx = None
        for k in ['old_name','original_name','old_defline','original_defline','original']:
            if k in h:
                old_idx = h[k]; break
        for k in ['new_name','new_defline','new','contig','contig_name']:
            if k in h:
                new_idx = h[k]; break
        if old_idx is not None and new_idx is not None:
            for r in data:
                if len(r) <= max(old_idx, new_idx):
                    continue
                old = clean(r[old_idx])
                new = clean(r[new_idx])
                if old and new:
                    mapping[new] = old
            return mapping

    for r in data:
        if len(r) < 2:
            continue
        a, b = clean(r[0]), clean(r[1])
        if not a or not b:
            continue
        a_simple = bool(re.match(r'^(contig|c|seq|scaffold)[_\-]?[0-9]+$', a, re.I))
        b_simple = bool(re.match(r'^(contig|c|seq|scaffold)[_\-]?[0-9]+$', b, re.I))
        if a_simple and not b_simple:
            mapping[a] = b
        elif b_simple and not a_simple:
            mapping[b] = a
        else:
            mapping[b] = a
    return mapping

name_map = load_name_map(NAME_MAP)

conn = sqlite3.connect(DB)
cur = conn.cursor()

if not table_exists(cur, 'hmm_hits'):
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
    available = ', '.join(r[0] for r in cur.fetchall())
    raise SystemExit(f"ERROR: table 'hmm_hits' not found in contigs.db. Available tables: {available}")

hmm_cols = columns(cur, 'hmm_hits')
gene_col_hmm = first_existing(hmm_cols, ['gene_callers_id', 'gene_id', 'gene_callers_id'])
hmm_name_col = first_existing(hmm_cols, ['hmm_name', 'hmm', 'gene_name', 'name', 'model_name', 'hmm_id'])
source_col = first_existing(hmm_cols, ['source', 'hmm_source', 'collection'])

if gene_col_hmm is None:
    raise SystemExit(f"ERROR: could not find gene_callers_id-like column in hmm_hits. Columns: {hmm_cols}")
if hmm_name_col is None:
    raise SystemExit(f"ERROR: could not find hmm_name-like column in hmm_hits. Columns: {hmm_cols}")

gene_table = None
for cand in ['genes_in_contigs', 'gene_calls', 'genes']:
    if table_exists(cur, cand):
        cols = columns(cur, cand)
        if first_existing(cols, ['gene_callers_id','gene_id']) and first_existing(cols, ['contig','contig_name','sequence','seq_id']):
            gene_table = cand
            break

if gene_table is None:
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
    available = ', '.join(r[0] for r in cur.fetchall())
    raise SystemExit(f"ERROR: could not find a gene-to-contig table. Available tables: {available}")

gene_cols = columns(cur, gene_table)
gene_col_gene = first_existing(gene_cols, ['gene_callers_id', 'gene_id'])
contig_col = first_existing(gene_cols, ['contig', 'contig_name', 'sequence', 'seq_id'])

q = f'''
SELECT g."{contig_col}", h."{hmm_name_col}", h."{gene_col_hmm}", {('h."' + source_col + '"') if source_col else "''"}
FROM hmm_hits h
JOIN "{gene_table}" g
ON h."{gene_col_hmm}" = g."{gene_col_gene}"
ORDER BY g."{contig_col}", h."{gene_col_hmm}", h."{hmm_name_col}"
'''

rows = []
skipped_source = 0
for contig, hmm_name, gene_id, source in cur.execute(q):
    contig = str(contig)
    out_contig = name_map.get(contig, contig)
    hmm_name = str(hmm_name)
    source = str(source) if source is not None else ''
    if source and not source_pattern.search(source):
        skipped_source += 1
        continue
    marker = f'{source}|{hmm_name}' if source else hmm_name
    rows.append((out_contig, marker, str(gene_id), source if source else 'NA'))

conn.close()
rows = sorted(set(rows))

if len(rows) == 0:
    raise SystemExit(
        'ERROR: no SCG HMM hits passed the HMM collection filter; verify anvi-run-hmms '
        'and SCG_SOURCE_REGEX'
    )

with open(OUT_TSV, 'w', newline='') as out:
    w = csv.writer(out, delimiter='\t')
    w.writerow(['contig', 'hmm_name', 'gene_callers_id', 'HMM_collection'])
    w.writerows(rows)

print(f'hmm_hits columns      : {", ".join(hmm_cols)}')
print(f'gene table used       : {gene_table}')
print(f'gene table columns    : {", ".join(gene_cols)}')
print(f'name_map entries      : {len(name_map)}')
print(f'SCG/HMM hits written  : {len(rows)}')
print(f'HMM hits excluded by collection filter: {skipped_source}')
print(f'Wrote                 : {OUT_TSV}')

PY

    mv -f "$TMP_OUT" "$OUT_TSV"

    echo
    head "$OUT_TSV" || true
    echo
    wc -l "$OUT_TSV" || true

    echo
    echo "[00_SCG] Finished: $(date)"
    echo "Outputs:"
    echo "  SCG TSV  : $OUT_TSV"
    echo "  DB       : $DB"
    echo "  FASTA    : $REFORMATTED"
    echo "  NAME_MAP : $NAME_MAP"
}

prep_circularity() {
    if [[ $# -lt 3 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        echo 'Usage: bash 00_prep.sh circularity MAG_DIR OUTPUT_DIR ASSEMBLER_FILE [ASSEMBLER_FILE ...]' >&2
        echo 'Run through: MAGverify prep-circularity -c FILE [-c FILE ...] -g MAG_DIR -o OUTPUT' >&2
        exit 2
    fi
    MAG_DIR="$1"
    OUTPUT_DIR="$2"
    shift 2
    [[ -d "$MAG_DIR" ]] || { echo "ERROR: MAG directory not found: $MAG_DIR" >&2; exit 2; }
    args=()
    for file in "$@"; do
        [[ -s "$file" ]] || { echo "ERROR: assembler file missing or empty: $file" >&2; exit 2; }
        args+=(-c "$file")
    done
    python3 "$SCRIPT_DIR/../scripts/prep_assembler_circularity.py" "${args[@]}" -g "$MAG_DIR" -o "$OUTPUT_DIR"
}

prep_taxonomy() {
    if [[ $# -ne 3 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        echo 'Usage: bash 00_prep.sh taxonomy GTDBTK_DIR MAG_DIR OUTPUT_TSV' >&2
        echo 'Run through: MAGverify prep-taxonomy -d GTDBTK_DIR -g MAG_DIR -o OUTPUT' >&2
        exit 2
    fi
    [[ -d "$1" ]] || { echo "ERROR: GTDB-Tk directory not found: $1" >&2; exit 2; }
    [[ -d "$2" ]] || { echo "ERROR: MAG directory not found: $2" >&2; exit 2; }
    mkdir -p "$(dirname "$3")"
    python3 "$SCRIPT_DIR/../scripts/prep_gtdbtk_summary_to_mag_taxonomy.py" \
        --gtdbtk-dir "$1" -g "$2" -o "$3"
}

subcommand="${1:-}"
[[ $# -gt 0 ]] && shift
case "$subcommand" in
    chimera) prep_chimera "$@" ;;
    scg) prep_scg "$@" ;;
    circularity) prep_circularity "$@" ;;
    taxonomy) prep_taxonomy "$@" ;;
    -h|--help) usage ;;
    "") usage >&2; exit 2 ;;
    *)
        echo "ERROR: unknown prep subcommand: $subcommand" >&2
        usage >&2
        exit 2
        ;;
esac
