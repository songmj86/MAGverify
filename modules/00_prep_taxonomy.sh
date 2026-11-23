#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    echo 'Usage: bash 00_prep_taxonomy.sh GTDBTK_DIR MAG_DIR OUTPUT_TSV' >&2
    echo 'Run through: MAGverify prep-taxonomy -d GTDBTK_DIR -g MAG_DIR -o OUTPUT' >&2
    exit 2
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -d "$1" ]] || { echo "ERROR: GTDB-Tk directory not found: $1" >&2; exit 2; }
[[ -d "$2" ]] || { echo "ERROR: MAG directory not found: $2" >&2; exit 2; }
mkdir -p "$(dirname "$3")"
python3 "$SCRIPT_DIR/../scripts/prep_gtdbtk_summary_to_mag_taxonomy.py" \
    --gtdbtk-dir "$1" -g "$2" -o "$3"
