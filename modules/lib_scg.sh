#!/usr/bin/env bash
set -euo pipefail

# Export anvi'o single-copy marker genes from an assembly.
# Usage: bash lib_scg.sh <ASSEMBLY_FA> <OUT_DIR> [THREADS]

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    cat <<'HELP'
Usage: bash lib_scg.sh ASSEMBLY_FASTA OUTPUT_DIR [THREADS]
Input: combined assembly FASTA.
Additional software: anvi'o (anvi-script-reformat-fasta, anvi-gen-contigs-database, anvi-run-hmms).
Output: OUTPUT_DIR/SCG_from_anvio.original_contig.tsv
HELP
    exit 0
fi
ASSEMBLY_FA="${1:?Usage: bash lib_scg.sh <ASSEMBLY_FA> <OUT_DIR> [THREADS]}"
OUT_DIR="${2:?Usage: bash lib_scg.sh <ASSEMBLY_FA> <OUT_DIR> [THREADS]}"
THREADS="${3:-${THREADS:-4}}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="$SCRIPT_DIR"
source "$SCRIPT_DIR/defaults.sh"
if [[ -n "${ANVIO_ENV_PREFIX:-}" && -d "$ANVIO_ENV_PREFIX/bin" ]]; then
    export PATH="$ANVIO_ENV_PREFIX/bin:$PATH"
fi

# Collections used by the small-genome proxy.
SCG_SOURCE_REGEX="${SCG_SOURCE_REGEX:-^(Bacteria_71|Archaea_76)$}"

mkdir -p "$OUT_DIR"

for t in anvi-script-reformat-fasta anvi-gen-contigs-database anvi-run-hmms python3; do
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
DB="$OUT_DIR/contigs.db"
OUT_TSV="$OUT_DIR/SCG_from_anvio.original_contig.tsv"
LOG="$OUT_DIR/00_SCG.log"
CACHE_META="$OUT_DIR/inputs.sha256"
TMP_CACHE="${CACHE_META}.tmp.$$"
TMP_OUT="${OUT_TSV}.tmp.$$"
trap 'rm -f "$TMP_CACHE" "$TMP_OUT"' EXIT

exec > >(tee -a "$LOG") 2>&1

echo "[00_SCG] Started: $(date)"
echo "Input assembly : $ASSEMBLY_FA"
echo "Output dir     : $OUT_DIR"
echo "Threads        : $THREADS"
echo "SCG sources    : $SCG_SOURCE_REGEX"
echo

# Fingerprint anvi'o preparation so the expensive DB/HMM steps are reused.

INPUT_FP="$(python3 - "$ASSEMBLY_FA" <<'PY_FP'
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
    'reformat_binary': stat(shutil.which('anvi-script-reformat-fasta')),
    'contigs_db_binary': stat(shutil.which('anvi-gen-contigs-database')),
    'hmms_binary': stat(shutil.which('anvi-run-hmms')),
    'method': 'simplify_names|contigs_db|run_hmms_v2',
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_FP
)"

PREP_REUSED=0
if [[ -s "$CACHE_META" && "$(<"$CACHE_META")" == "$INPUT_FP" \
      && -s "$REFORMATTED" && -s "$NAME_MAP" && -s "$DB" ]]; then
    PREP_REUSED=1
    echo "Fingerprint-matched anvi'o DB/HMM preparation found."
    echo "Reusing: $DB"
else
    echo "Removing stale or incomplete anvi'o preparation outputs."
    rm -f "$REFORMATTED" "$NAME_MAP" "$DB" "$CACHE_META" "$OUT_TSV"
fi

# Steps 1-3. Prepare anvi'o DB and HMM hits when cache cannot be reused.

if (( PREP_REUSED == 0 )); then
    echo
    echo "[1/4] Reformatting FASTA for anvi'o simple deflines"
    anvi-script-reformat-fasta "$ASSEMBLY_FA" \
        --simplify-names \
        --report-file "$NAME_MAP" \
        -o "$REFORMATTED"

    echo
    echo "[2/4] Generating anvi'o contigs.db"
    anvi-gen-contigs-database \
        -f "$REFORMATTED" \
        -o "$DB" \
        -n "LR_assembly" \
        -T "$THREADS"

    echo
    echo "[3/4] Running anvi-run-hmms"
    # Only the collections exported below are needed; skip the other installed HMMs.
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

# Step 4. Export Bacteria_71 / Archaea_76 hits with original contig names.

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
    # partial fallback
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

    # Fallback: assume two-column map. Try to identify simplified side.
    for r in data:
        if len(r) < 2:
            continue
        a, b = clean(r[0]), clean(r[1])
        if not a or not b:
            continue
        # Common anvi simplified names: contig_000000000001, c_000000000001, etc.
        a_simple = bool(re.match(r'^(contig|c|seq|scaffold)[_\-]?[0-9]+$', a, re.I))
        b_simple = bool(re.match(r'^(contig|c|seq|scaffold)[_\-]?[0-9]+$', b, re.I))
        if a_simple and not b_simple:
            mapping[a] = b
        elif b_simple and not a_simple:
            mapping[b] = a
        else:
            # Last resort: many anvi reports are old\tnew, so map new->old.
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

# Find a gene-to-contig table.
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
        'ERROR: no SCG HMM hits passed the source filter; verify anvi-run-hmms '
        'and SCG_SOURCE_REGEX'
    )

with open(OUT_TSV, 'w', newline='') as out:
    w = csv.writer(out, delimiter='\t')
    w.writerow(['contig', 'hmm_name', 'gene_callers_id', 'source'])
    w.writerows(rows)

print(f'hmm_hits columns      : {", ".join(hmm_cols)}')
print(f'gene table used       : {gene_table}')
print(f'gene table columns    : {", ".join(gene_cols)}')
print(f'name_map entries      : {len(name_map)}')
print(f'SCG/HMM hits written  : {len(rows)}')
print(f'HMM hits excluded by source filter: {skipped_source}')
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
