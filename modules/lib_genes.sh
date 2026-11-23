#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat >&2 <<'HELP'
Usage: bash lib_genes.sh ASSEMBLY_FASTA OUTPUT_DIR [THREADS]
       bash lib_genes.sh --fingerprint ASSEMBLY_FASTA
       bash lib_genes.sh --version
Output: OUTPUT_DIR/genes.faa, genes.gff, genes.coordinates.tsv
HELP
    exit "${1:-2}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
LINKER="$SCRIPT_DIR/../scripts/prep_genes.py"
METHOD="pyrodigal_meta_v1"

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage 0

PYRODIGAL_BIN="${PYRODIGAL_BIN:-$(command -v pyrodigal 2>/dev/null || true)}"
[[ -n "$PYRODIGAL_BIN" && -x "$PYRODIGAL_BIN" ]] || {
    echo "ERROR: pyrodigal not found. Activate the MAGverify Conda environment first." >&2
    exit 1
}

fingerprint() {
    python3 - "$1" "$PYRODIGAL_BIN" "$LINKER" "$METHOD" <<'PY_FP'
import hashlib, json, os, sys
assembly, binary, linker, method = sys.argv[1:]
def stat(path):
    value = os.stat(path)
    return [os.path.realpath(path), value.st_size, value.st_mtime_ns]
payload = {'assembly': stat(assembly), 'predictor': stat(binary),
           'linker': stat(linker), 'method': method}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_FP
}

if [[ "${1:-}" == "--version" ]]; then
    version="$("$PYRODIGAL_BIN" --version 2>/dev/null | grep -oE '[0-9]+(\.[0-9]+)+' | head -n 1 || true)"
    printf '%s\n' "${version:-NA}"
    exit 0
fi

if [[ "${1:-}" == "--fingerprint" ]]; then
    [[ -s "${2:-}" ]] || usage
    fingerprint "$2"
    exit 0
fi

ASSEMBLY_FA="${1:-}"
OUT_DIR="${2:-}"
THREADS="${3:-${THREADS:-4}}"
[[ -n "$ASSEMBLY_FA" && -n "$OUT_DIR" ]] || usage
[[ -s "$ASSEMBLY_FA" ]] || { echo "ERROR: assembly FASTA not found or empty: $ASSEMBLY_FA" >&2; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: THREADS must be a positive integer: $THREADS" >&2; exit 2; }

FAA="$OUT_DIR/genes.faa"
GFF="$OUT_DIR/genes.gff"
COORD="$OUT_DIR/genes.coordinates.tsv"
CACHE="$OUT_DIR/.cache/genes.inputs.sha256"
mkdir -p "$OUT_DIR/.cache"

if command -v flock >/dev/null 2>&1; then
    exec 9>"$OUT_DIR/.cache/genes.lock"
    flock 9
fi

FP="$(fingerprint "$ASSEMBLY_FA")"
if [[ -s "$CACHE" && "$(<"$CACHE")" == "$FP" && -s "$FAA" && -s "$GFF" && -s "$COORD" ]]; then
    echo "Gene prediction: $OUT_DIR (reused)"
    exit 0
fi

echo "Gene prediction: pyrodigal -p meta, $THREADS threads -> $OUT_DIR"
rm -f "$CACHE"
TMP="$(mktemp -d "$OUT_DIR/.genes.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
"$PYRODIGAL_BIN" -p meta -j "$THREADS" -i "$ASSEMBLY_FA" \
    -a "$TMP/genes.faa" -f gff -o "$TMP/genes.gff"
python3 "$LINKER" link "$TMP/genes.faa" "$TMP/genes.gff" "$TMP/genes.coordinates.tsv"
[[ -s "$TMP/genes.coordinates.tsv" ]] || { echo "ERROR: pyrodigal produced no genes" >&2; exit 1; }
mv -f "$TMP/genes.faa" "$FAA"
mv -f "$TMP/genes.gff" "$GFF"
mv -f "$TMP/genes.coordinates.tsv" "$COORD"
printf '%s\n' "$FP" > "$CACHE.tmp.$$"
mv -f "$CACHE.tmp.$$" "$CACHE"
echo "Predicted genes: $(wc -l < "$COORD")"
