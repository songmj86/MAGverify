#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    echo 'Usage: bash 00_prep_circularity.sh MAG_DIR OUTPUT_DIR ASSEMBLER_FILE [ASSEMBLER_FILE ...]' >&2
    echo 'Run through: MAGverify prep-circularity -c FILE [-c FILE ...] -g MAG_DIR -o OUTPUT' >&2
    exit 2
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
