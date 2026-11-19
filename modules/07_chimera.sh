#!/usr/bin/env bash
# Gene taxonomy preprocessing is embedded here to keep chimera analysis self-contained.
if [[ "${1:-}" == "__taxonomy_internal" || "${1:-}" == "__taxonomy_scope" ]]; then
    _internal_mode="$1"
    shift
    if [[ "$_internal_mode" == "__taxonomy_scope" ]]; then
        set -- __scope "$@"
    fi
set -Eeuo pipefail

SCRIPT_ID="00_gene_taxonomy"
GENE_PREDICTION_METHOD_VERSION="full_assembly_meta_gene_prediction_v1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SELF="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"

detect_default_threads() {
    local value
    if command -v nproc >/dev/null 2>&1; then
        value="$(nproc --all)"
    else
        value="$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')"
    fi
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || value=1
    printf '%s\n' "$value"
}

public_usage() {
    cat >&2 <<'EOF'
Usage:
  MMSEQS_TAX_DB=<GTDB_R232_DB_PREFIX> MAG_CONTIG_MAP=<MAP_TSV> \
    MAGverify chimera -a ASSEMBLIES --mag-dir MAGs --database MMSEQS_PREFIX ...

  Internal: bash 07_chimera.sh __taxonomy_internal <genome|contig> \
    <ASSEMBLY_FA> <GTDB_R232_DB_PREFIX> <MAG_CONTIG_MAP> <OUT_DIR> [THREADS]

  Internal taxonomy preprocessing is called by the chimera analysis.
  THREADS comes from the MAGverify chimera command.

Modes:
  genome  Produce the MAG-associated-contig GTDB R232 taxonomy table only.
  contig  Reuse a compatible MAG result when available, classify non-MAG contigs,
          and produce the merged all-assembly-contig GTDB R232 taxonomy table.
          Non-MAG contigs with <8 predicted genes are skipped before MMseqs2
          by default because they cannot satisfy >=4 classified genes per side.

Outputs:
  gene_taxonomy.mmseqs2.GTDB_R232.{MAG,nonMAG,all}_contigs.tsv

The first run predicts genes once across the full assembly. Runs that use the
same OUT_DIR, assembly, and predictor settings reuse this validated prediction
across genome/contig modes.

The MMseqs2 database argument is the GTDB R232 database prefix, not a directory.
EOF
    exit "${1:-2}"
}

scope_usage() {
    echo "ERROR: invalid internal scope invocation" >&2
    public_usage
}

# Public genome/contig controller. The same file is invoked with the private
# __scope mode for MAG and non-MAG work, so no second taxonomy shell is needed.
if [[ "${1:-}" != "__scope" ]]; then
    DEFAULT_THREADS="$(detect_default_threads)"
    if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
        public_usage 0
    fi
    if (( $# == 2 )); then
        CONFIG_DIR="$SCRIPT_DIR"
        source "$SCRIPT_DIR/defaults.sh"
        MAP_DEFAULT="$OUT/00_prep/MAG_contig_map.tsv"
        set -- "${TAXONOMY_MODE:-contig}" "$1" "${MMSEQS_TAX_DB:-}" \
            "${MAG_CONTIG_MAP:-$MAP_DEFAULT}" "$2" "${THREADS:-$DEFAULT_THREADS}"
    fi
    [[ $# -ge 5 && $# -le 6 ]] || public_usage

    MODE="$1"
    ASSEMBLY_FA="$2"
    MMSEQS_TAX_DB="$3"
    MAG_CONTIG_MAP="$4"
    OUT_DIR="$5"
    THREADS="${6:-$DEFAULT_THREADS}"

    [[ "$MODE" == "genome" || "$MODE" == "contig" ]] || public_usage
    [[ -n "$MMSEQS_TAX_DB" && -n "$MAG_CONTIG_MAP" ]] || {
        echo 'ERROR: MMSEQS_TAX_DB and MAG_CONTIG_MAP must be set for taxonomy.' >&2
        exit 2
    }
    [[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: THREADS must be a positive integer: $THREADS" >&2
        exit 2
    }

    MERGER="$SCRIPT_DIR/../scripts/chimera_merge_gene_taxonomy_tables.py"
    [[ -s "$MERGER" ]] || {
        echo "ERROR: missing required helper: $MERGER" >&2
        exit 1
    }
    for tool in python3 awk; do
        command -v "$tool" >/dev/null 2>&1 || {
            echo "ERROR: required tool not found in PATH: $tool" >&2
            exit 1
        }
    done

    MMSEQS_TAX_DB="${MMSEQS_TAX_DB%.dbtype}"
    OUTPUT_INFIX=".GTDB_R232"

    mkdir -p "$OUT_DIR"
    OUT_DIR="$(cd "$OUT_DIR" && pwd -P)"
    MAG_OUT="$OUT_DIR/gene_taxonomy.mmseqs2${OUTPUT_INFIX}.MAG_contigs.tsv"
    NONMAG_OUT="$OUT_DIR/gene_taxonomy.mmseqs2${OUTPUT_INFIX}.nonMAG_contigs.tsv"
    ALL_OUT="$OUT_DIR/gene_taxonomy.mmseqs2${OUTPUT_INFIX}.all_contigs.tsv"

    run_genome() {
        echo "[genome] MAG-associated contigs only (complete gene taxonomy; no chimera prefilter)"
        TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG=1 bash "$SELF" __taxonomy_scope mag \
            "$ASSEMBLY_FA" "$MMSEQS_TAX_DB" "$MAG_OUT" \
            "$MAG_CONTIG_MAP" "$THREADS"
        echo "MAG result: $MAG_OUT"
    }

    run_contig() {
        MAG_META="${MAG_OUT}.run_metadata.tsv"
        if [[ -s "$MAG_OUT" && -s "$MAG_META" ]] && python3 - "$MAG_OUT" "$MAG_META" "$MMSEQS_TAX_DB" "$ASSEMBLY_FA" "$MAG_CONTIG_MAP" <<'PY_VALIDATE_GENOME'
import csv
import os
import sys

out_path, meta_path, database, assembly_path, map_path = sys.argv[1:]
expected_header = ['contig','start','end','domain','phylum','class','order','family','gene_id']
try:
    with open(out_path, encoding='utf-8-sig') as handle:
        header = handle.readline().rstrip('\n\r').split('\t')
    with open(meta_path, encoding='utf-8-sig', newline='') as handle:
        meta = {row['key']: row['value'] for row in csv.DictReader(handle, delimiter='\t')}
except (OSError, csv.Error):
    raise SystemExit(1)
if header != expected_header:
    raise SystemExit(1)
expected = ('GTDB', 'GTDB_R232', 'disabled_empty', '<empty>')
for key, value in zip(
    ('taxonomy_source', 'taxonomy_release', 'blacklist_policy', 'blacklist_value'),
    expected,
):
    if meta.get(key) != value:
        raise SystemExit(1)
if meta.get('output_schema') != 'nine_column_rank_names_v1':
    raise SystemExit(1)
if os.path.realpath(meta.get('mmseqs_database', '')) != os.path.realpath(database):
    raise SystemExit(1)
# Cheap staleness guard for the compatible-result fast path. The full internal
# fingerprint remains authoritative whenever this check fails and genome mode
# is rebuilt.
try:
    output_mtime = os.path.getmtime(out_path)
    if output_mtime < os.path.getmtime(assembly_path) or output_mtime < os.path.getmtime(map_path):
        raise SystemExit(1)
    cache_dir = meta.get('gene_prediction_cache', '')
    fp = meta.get('gene_prediction_fingerprint', '')
    fp_file = os.path.join(cache_dir, 'inputs.sha256')
    if not cache_dir or not fp or not os.path.isfile(fp_file):
        raise SystemExit(1)
    with open(fp_file, encoding='utf-8') as handle:
        if handle.read().strip() != fp:
            raise SystemExit(1)
except OSError:
    raise SystemExit(1)
PY_VALIDATE_GENOME
        then
            echo "[contig] Reusing compatible existing MAG taxonomy: $MAG_OUT"
        else
            echo "[contig] Compatible MAG taxonomy not found; rebuilding genome scope"
            run_genome
        fi
        [[ -s "$MAG_OUT" && -s "$MAG_META" ]] || {
            echo "ERROR: genome-mode MAG taxonomy or metadata is missing" >&2
            exit 1
        }

        GENOME_PREDICTOR="$(awk -F'\t' '$1=="gene_predictor"{print $2}' "$MAG_META")"
        [[ -n "$GENOME_PREDICTOR" ]] || {
            echo "ERROR: genome-mode gene_predictor metadata is blank" >&2
            exit 1
        }
        if [[ "${GENE_PREDICTOR:-auto}" == "auto" ]]; then
            export GENE_PREDICTOR="$GENOME_PREDICTOR"
            echo "[contig] Reusing genome-mode predictor backend: $GENOME_PREDICTOR"
        elif [[ "$GENE_PREDICTOR" != "$GENOME_PREDICTOR" ]]; then
            echo "ERROR: contig-mode GENE_PREDICTOR=$GENE_PREDICTOR differs from genome-mode=$GENOME_PREDICTOR" >&2
            exit 1
        fi

        CHIMERA_MIN_CLASSIFIED_GENES="${TAX_MIN_CLASSIFIED_GENES:-4}"
        [[ "$CHIMERA_MIN_CLASSIFIED_GENES" =~ ^[1-9][0-9]*$ ]] || {
            echo "ERROR: TAX_MIN_CLASSIFIED_GENES must be a positive integer" >&2
            exit 2
        }
        SAFE_PREFILTER_MAX=$(( 2 * CHIMERA_MIN_CLASSIFIED_GENES ))
        TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG="${TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG:-$SAFE_PREFILTER_MAX}"
        [[ "$TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG" =~ ^[1-9][0-9]*$ ]] || {
            echo "ERROR: TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG must be a positive integer" >&2
            exit 2
        }
        if (( TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG > SAFE_PREFILTER_MAX )); then
            echo "ERROR: taxonomy prefilter is too stringent for TAX_MIN_CLASSIFIED_GENES=$CHIMERA_MIN_CLASSIFIED_GENES" >&2
            echo "Maximum lossless threshold is $SAFE_PREFILTER_MAX predicted genes/contig." >&2
            exit 2
        fi
        export TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG
        echo "[contig] Classifying non-MAG contigs (lossless MMseqs2 prefilter: >=${TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG} predicted genes/contig; chimera minimum classified genes/side=${CHIMERA_MIN_CLASSIFIED_GENES})"
        bash "$SELF" __taxonomy_scope nonmag \
            "$ASSEMBLY_FA" "$MMSEQS_TAX_DB" "$NONMAG_OUT" \
            "$MAG_CONTIG_MAP" "$THREADS"

        NONMAG_META="${NONMAG_OUT}.run_metadata.tsv"
        python3 - "$MAG_META" "$NONMAG_META" <<'PY_COMPARE_META'
import csv
import os
import sys

def load(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        return {row['key']: row['value'] for row in csv.DictReader(handle, delimiter='\t')}

mag = load(sys.argv[1])
nonmag = load(sys.argv[2])
for key in (
    'gene_predictor', 'predictor_realpath', 'predictor_signature',
    'gene_prediction_fingerprint', 'mmseqs_database', 'taxonomy_source',
    'taxonomy_release', 'taxonomy_db_label', 'gtdb_db_label',
    'blacklist_policy', 'blacklist_value', 'output_schema',
):
    if key == 'mmseqs_database':
        same = os.path.realpath(mag.get(key, '')) == os.path.realpath(nonmag.get(key, ''))
    else:
        same = mag.get(key) == nonmag.get(key)
    if not same:
        raise SystemExit(
            f'ERROR: genome/contig metadata differ for {key}: '
            f'{mag.get(key)!r} versus {nonmag.get(key)!r}'
        )
PY_COMPARE_META

        echo "[contig] Merging MAG and non-MAG results"
        python3 "$MERGER" "$ALL_OUT" "$MAG_OUT" "$NONMAG_OUT"
        echo "Non-MAG result : $NONMAG_OUT"
        echo "All-contig result: $ALL_OUT"
    }

    case "$MODE" in
        genome) run_genome ;;
        contig) run_contig ;;
    esac
    exit 0
fi

shift  # private __scope marker

# Internal scope engine: build one nine-column GTDB R232 gene-taxonomy table
# with Prodigal/Pyrodigal and MMseqs2. The GTDB MMseqs2 database uses the
# empty blacklist policy so NCBI-taxid blacklist entries do not suppress GTDB
# taxonomy assignments.
#
# Gene-prediction backends:
#   GENE_PREDICTOR=auto               prefer parallel Prodigal, otherwise Pyrodigal
#   GENE_PREDICTOR=pyrodigal          Pyrodigal CLI with -j threading
#   GENE_PREDICTOR=prodigal_parallel  intact-contig chunks, parallel Prodigal jobs
#   GENE_PREDICTOR=prodigal_serial    one Prodigal process

TAXONOMY_SOURCE="GTDB"
TAXONOMY_RELEASE="GTDB_R232"
BLACKLIST_POLICY="disabled_empty"
BLACKLIST_VALUE="<empty>"
APPLY_EMPTY_BLACKLIST="Y"

[[ "${1:-}" == "mag" || "${1:-}" == "nonmag" || "${1:-}" == "all" ]] \
    || scope_usage
SCOPE="$1"
shift

if command -v nproc >/dev/null 2>&1; then
    DEFAULT_THREADS="$(nproc --all)"
else
    DEFAULT_THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')"
fi
[[ "$DEFAULT_THREADS" =~ ^[1-9][0-9]*$ ]] || DEFAULT_THREADS=1

case "$SCOPE" in
    mag|nonmag)
        [[ $# -ge 4 && $# -le 5 ]] || scope_usage
        ASSEMBLY_FA="$1"
        MMSEQS_TAX_DB="$2"
        OUT_TSV="$3"
        MAG_CONTIG_MAP="$4"
        THREADS="${5:-$DEFAULT_THREADS}"
        ;;
    all)
        [[ $# -ge 3 && $# -le 4 ]] || scope_usage
        ASSEMBLY_FA="$1"
        MMSEQS_TAX_DB="$2"
        OUT_TSV="$3"
        MAG_CONTIG_MAP="NA"
        THREADS="${4:-$DEFAULT_THREADS}"
        ;;
    *) scope_usage ;;
esac
PRODIGAL_LINKER="$SCRIPT_DIR/../scripts/chimera_link_prodigal_accessions.py"
FASTA_SELECTOR="$SCRIPT_DIR/../scripts/chimera_select_sequences.py"
GENE_SCOPE_SELECTOR="$SCRIPT_DIR/../scripts/chimera_select_sequences.py"
MMSEQS_TAXONOMY_MERGER="$SCRIPT_DIR/../scripts/chimera_merge_mmseqs_taxonomy_with_coordinates.py"
TAXONOMY_PREFILTER="$SCRIPT_DIR/../scripts/chimera_filter_gene_predictions_by_min_contig_gene_count.py"

for file in \
    "$PRODIGAL_LINKER" "$FASTA_SELECTOR" "$GENE_SCOPE_SELECTOR" \
    "$MMSEQS_TAXONOMY_MERGER" "$TAXONOMY_PREFILTER"; do
    [[ -s "$file" ]] || {
        echo "ERROR: missing helper script: $file" >&2
        exit 1
    }
done

MMSEQS_ENV="${MMSEQS_ENV:-${CONDA_PREFIX:-}}"
MMSEQS_BIN="${MMSEQS_BIN:-$(command -v mmseqs 2>/dev/null || true)}"
if [[ -n "$MMSEQS_ENV" && -x "$MMSEQS_ENV/bin/mmseqs" ]]; then
    MMSEQS_BIN="$MMSEQS_ENV/bin/mmseqs"
fi
[[ -x "$MMSEQS_BIN" ]] || {
    echo "ERROR: MMSEQS_BIN is not executable: $MMSEQS_BIN" >&2
    exit 1
}
if [[ -n "$MMSEQS_ENV" ]]; then export PATH="${MMSEQS_ENV}/bin:${PATH}"; fi

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: THREADS must be a positive integer: $THREADS" >&2
    exit 2
}

GENE_PREDICTOR="${GENE_PREDICTOR:-auto}"
GENE_PRED_THREADS="${GENE_PRED_THREADS:-$THREADS}"
PYRODIGAL_POOL="${PYRODIGAL_POOL:-thread}"
PRODIGAL_JOBS="${PRODIGAL_JOBS:-$THREADS}"
PRODIGAL_CHUNK_BP="${PRODIGAL_CHUNK_BP:-100000000}"
MMSEQS_THREADS="${MMSEQS_THREADS:-$THREADS}"
TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG="${TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG:-1}"
KEEP_FAILED_WORK="${KEEP_FAILED_WORK:-Y}"
PYRODIGAL_BIN="${PYRODIGAL_BIN:-$(command -v pyrodigal 2>/dev/null || true)}"
PRODIGAL_BIN="${PRODIGAL_BIN:-$(command -v prodigal 2>/dev/null || true)}"

for value_name in GENE_PRED_THREADS PRODIGAL_JOBS PRODIGAL_CHUNK_BP MMSEQS_THREADS TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG; do
    value="${!value_name}"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: $value_name must be a positive integer: $value" >&2
        exit 2
    }
done
[[ "$PYRODIGAL_POOL" == "thread" || "$PYRODIGAL_POOL" == "process" ]] || {
    echo "ERROR: PYRODIGAL_POOL must be thread or process" >&2
    exit 2
}
[[ "$KEEP_FAILED_WORK" == "Y" || "$KEEP_FAILED_WORK" == "N" ]] || {
    echo "ERROR: KEEP_FAILED_WORK must be Y or N" >&2
    exit 2
}
if (( PRODIGAL_JOBS > THREADS )); then
    PRODIGAL_JOBS="$THREADS"
fi

export OMP_NUM_THREADS="$MMSEQS_THREADS"

case "$GENE_PREDICTOR" in
    auto)
        if [[ -n "$PRODIGAL_BIN" && -x "$PRODIGAL_BIN" ]]; then
            GENE_PREDICTOR="prodigal_parallel"
        elif [[ -n "$PYRODIGAL_BIN" && -x "$PYRODIGAL_BIN" ]]; then
            GENE_PREDICTOR="pyrodigal"
        else
            echo "ERROR: neither pyrodigal nor prodigal is available" >&2
            exit 1
        fi
        ;;
    pyrodigal)
        [[ -n "$PYRODIGAL_BIN" && -x "$PYRODIGAL_BIN" ]] || {
            echo "ERROR: GENE_PREDICTOR=pyrodigal but pyrodigal is unavailable" >&2
            exit 1
        }
        ;;
    prodigal_parallel|prodigal_serial)
        [[ -n "$PRODIGAL_BIN" && -x "$PRODIGAL_BIN" ]] || {
            echo "ERROR: $GENE_PREDICTOR selected but prodigal is unavailable" >&2
            exit 1
        }
        ;;
    *)
        echo "ERROR: unsupported GENE_PREDICTOR: $GENE_PREDICTOR" >&2
        exit 2
        ;;
esac

for tool in python3 mktemp grep sort; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: required tool not found in PATH: $tool" >&2
        exit 1
    }
done
[[ "$GENE_PREDICTOR" != "prodigal_parallel" ]] || command -v xargs >/dev/null 2>&1 || {
    echo "ERROR: xargs is needed for prodigal_parallel" >&2
    exit 1
}

[[ -s "$ASSEMBLY_FA" ]] || {
    echo "ERROR: missing or empty ASSEMBLY_FA: $ASSEMBLY_FA" >&2
    exit 1
}
if [[ "$SCOPE" != "all" ]]; then
    [[ -s "$MAG_CONTIG_MAP" ]] || {
        echo "ERROR: missing or empty MAG_CONTIG_MAP: $MAG_CONTIG_MAP" >&2
        exit 1
    }
fi

MMSEQS_TAX_DB="${MMSEQS_TAX_DB%.dbtype}"
TAXONOMY_DB_LABEL="$(basename "$MMSEQS_TAX_DB")"
GTDB_DB_LABEL="$TAXONOMY_DB_LABEL"
ALLOW_TAXONOMY_DB_MISMATCH="${ALLOW_TAXONOMY_DB_MISMATCH:-${ALLOW_NON_GTDB_DB:-N}}"
[[ "$ALLOW_TAXONOMY_DB_MISMATCH" == "Y" || "$ALLOW_TAXONOMY_DB_MISMATCH" == "N" ]] || {
    echo "ERROR: ALLOW_TAXONOMY_DB_MISMATCH must be Y or N" >&2
    exit 2
}
if [[ "$ALLOW_TAXONOMY_DB_MISMATCH" != "Y" && ! "$TAXONOMY_DB_LABEL" =~ [Gg][Tt][Dd][Bb] ]]; then
    echo "ERROR: GTDB R232 mode requires a GTDB database prefix: $MMSEQS_TAX_DB" >&2
    echo "Set ALLOW_TAXONOMY_DB_MISMATCH=Y only for an intentionally renamed GTDB database." >&2
    exit 1
fi
[[ ! -d "$MMSEQS_TAX_DB" ]] || {
    echo "ERROR: supply the MMseqs2 database prefix, not a directory: $MMSEQS_TAX_DB" >&2
    exit 1
}
[[ -s "${MMSEQS_TAX_DB}.dbtype" ]] || {
    echo "ERROR: MMseqs2 database marker is missing: ${MMSEQS_TAX_DB}.dbtype" >&2
    exit 1
}

[[ -s "$MMSEQS_TAX_DB" && -s "${MMSEQS_TAX_DB}.index" ]] || {
    echo "ERROR: incomplete MMseqs2 sequence database prefix: $MMSEQS_TAX_DB" >&2
    exit 1
}

# MMseqs2 18 can print "Unknown" for a valid database through the optional
# `dbtype` helper. Read the canonical four-byte .dbtype marker first instead.
# Sequence DB type 0 is amino acid and type 1 is nucleotide.
DB_TYPE_CODE="$(python3 - "${MMSEQS_TAX_DB}.dbtype" <<'PY_DBTYPE'
import struct
import sys
from pathlib import Path
path = Path(sys.argv[1])
data = path.read_bytes()
if len(data) < 4:
    print('NA')
else:
    print(struct.unpack('<I', data[:4])[0] & 0xFFFF)
PY_DBTYPE
)"
DB_TYPE_RAW="$("$MMSEQS_BIN" dbtype "$MMSEQS_TAX_DB" 2>&1 || true)"
case "$DB_TYPE_CODE" in
    0) DB_TYPE="Aminoacid" ;;
    1)
        echo "ERROR: MMseqs2 target database is nucleotide, but protein genes require an amino-acid database: $MMSEQS_TAX_DB" >&2
        exit 1
        ;;
    *)
        if grep -qi 'amino' <<<"$DB_TYPE_RAW"; then
            DB_TYPE="Aminoacid"
        else
            DB_TYPE="Unresolved_by_dbtype_preflight"
            echo "WARNING: target DB type could not be confirmed from ${MMSEQS_TAX_DB}.dbtype (code=$DB_TYPE_CODE) or mmseqs dbtype: $DB_TYPE_RAW" >&2
        fi
        ;;
esac
if [[ ! -s "${MMSEQS_TAX_DB}_taxonomy" ]]; then
    for suffix in _mapping _nodes.dmp _names.dmp _merged.dmp; do
        [[ -s "${MMSEQS_TAX_DB}${suffix}" ]] || {
            echo "ERROR: taxonomy database is incomplete: ${MMSEQS_TAX_DB}${suffix}" >&2
            exit 1
        }
    done
fi

mkdir -p "$(dirname "$OUT_TSV")"
OUT_DIR="$(cd "$(dirname "$OUT_TSV")" && pwd -P)"
OUT_TSV="$OUT_DIR/$(basename "$OUT_TSV")"
CACHE_META="${OUT_TSV}.inputs.sha256"
RUN_META="${OUT_TSV}.run_metadata.tsv"
TMP_OUT="${OUT_TSV}.tmp.$$"
TMP_CACHE="${CACHE_META}.tmp.$$"
TMP_RUN_META="${RUN_META}.tmp.$$"
WORK=""
PRED_CACHE_TMP_FAA=""
PRED_CACHE_TMP_COORD=""
PRED_CACHE_TMP_META=""
PRED_CACHE_TMP_FP=""
SCOPE_CACHE_TMP_FAA=""
SCOPE_CACHE_TMP_COORD=""
SCOPE_CACHE_TMP_META=""
SCOPE_CACHE_TMP_FP=""

cleanup() {
    status=$?
    rm -f -- "$TMP_OUT" "$TMP_CACHE" "$TMP_RUN_META"
    for cache_tmp in \
        "$PRED_CACHE_TMP_FAA" "$PRED_CACHE_TMP_COORD" \
        "$PRED_CACHE_TMP_META" "$PRED_CACHE_TMP_FP" \
        "$SCOPE_CACHE_TMP_FAA" "$SCOPE_CACHE_TMP_COORD" \
        "$SCOPE_CACHE_TMP_META" "$SCOPE_CACHE_TMP_FP"; do
        [[ -z "$cache_tmp" ]] || rm -f -- "$cache_tmp"
    done
    if [[ -n "$WORK" && -d "$WORK" ]]; then
        if (( status != 0 )) && [[ "$KEEP_FAILED_WORK" == "Y" ]]; then
            echo "WARNING: retaining failed-run work directory for recovery: $WORK" >&2
            if [[ -s "${GENE_COORD:-}" && -s "${TAX_RAW:-}" ]]; then
                echo "Recovery command:" >&2
                printf 'python3 %q %q %q %q\n' \
                    "$MMSEQS_TAXONOMY_MERGER" "$GENE_COORD" "$TAX_RAW" "$OUT_TSV" >&2
            fi
        else
            case "$WORK" in
                "$OUT_DIR"/.mmseqs2_gene_taxonomy.*) rm -rf -- "$WORK" ;;
                *) echo "WARNING: refusing to remove unexpected work path: $WORK" >&2 ;;
            esac
        fi
    fi
    return "$status"
}
trap cleanup EXIT

PREDICTOR_BIN="$PRODIGAL_BIN"
[[ "$GENE_PREDICTOR" != "pyrodigal" ]] || PREDICTOR_BIN="$PYRODIGAL_BIN"
PREDICTOR_REALPATH="$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$PREDICTOR_BIN")"
PREDICTOR_SIGNATURE="$(python3 - "$PREDICTOR_BIN" <<'PY_PREDICTOR_SIG'
import hashlib
import os
import sys
path = sys.argv[1]
stat = os.stat(path)
payload = f'{os.path.realpath(path)}\0{stat.st_size}\0{stat.st_mtime_ns}'
print(hashlib.sha256(payload.encode()).hexdigest())
PY_PREDICTOR_SIG
)"

# This fingerprint deliberately excludes MAG/non-MAG scope, the MAG map, and
# the MMseqs2 taxonomy database. Prodigal runs once on the complete
# assembly and its predictions are partitioned afterwards.
PREDICTION_FP="$(python3 - \
    "$ASSEMBLY_FA" "$PRODIGAL_LINKER" "$PREDICTOR_BIN" \
    "$GENE_PREDICTOR" "$GENE_PRED_THREADS" "$PYRODIGAL_POOL" \
    "$PRODIGAL_JOBS" "$PRODIGAL_CHUNK_BP" \
    "$GENE_PREDICTION_METHOD_VERSION" <<'PY_PREDICTION_FP'
import hashlib
import json
import os
import sys

(
    assembly, linker, predictor_bin, predictor, predictor_threads,
    pyrodigal_pool, prodigal_jobs, prodigal_chunk_bp, method_version,
) = sys.argv[1:]

def stat(path):
    value = os.stat(path)
    return [os.path.realpath(path), value.st_size, value.st_mtime_ns]

payload = {
    'assembly': stat(assembly),
    'linker': stat(linker),
    'predictor_binary': stat(predictor_bin),
    'predictor': predictor,
    'predictor_threads': int(predictor_threads),
    'pyrodigal_pool': pyrodigal_pool,
    'prodigal_jobs': int(prodigal_jobs),
    'prodigal_chunk_bp': int(prodigal_chunk_bp),
    'method': method_version,
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_PREDICTION_FP
)"

PRED_CACHE_ROOT="${GENE_PREDICTION_CACHE_DIR:-$OUT_DIR/gene_prediction_cache}"
PRED_CACHE_ENTRY="$PRED_CACHE_ROOT/all/$PREDICTION_FP"
PRED_CACHE_FAA="$PRED_CACHE_ENTRY/genes.all.faa"
PRED_CACHE_COORD="$PRED_CACHE_ENTRY/gene_coordinates.all.tsv"
PRED_CACHE_META="$PRED_CACHE_ENTRY/run_metadata.tsv"
PRED_CACHE_FP="$PRED_CACHE_ENTRY/inputs.sha256"

SCOPE_SELECTION_FP="all"
SCOPE_CACHE_ENTRY="$PRED_CACHE_ENTRY"
SCOPE_CACHE_FAA="$PRED_CACHE_FAA"
SCOPE_CACHE_COORD="$PRED_CACHE_COORD"
SCOPE_CACHE_META="$PRED_CACHE_META"
SCOPE_CACHE_FP="$PRED_CACHE_FP"
if [[ "$SCOPE" != "all" ]]; then
    SCOPE_SELECTION_FP="$(python3 - \
        "$SCOPE" "$ASSEMBLY_FA" "$MAG_CONTIG_MAP" "$PREDICTION_FP" \
        "$GENE_SCOPE_SELECTOR" "$FASTA_SELECTOR" <<'PY_SCOPE_FP'
import hashlib
import json
import os
import sys

scope, assembly, map_file, prediction_fp, scope_selector, map_reader = sys.argv[1:]

def stat(path):
    value = os.stat(path)
    return [os.path.realpath(path), value.st_size, value.st_mtime_ns]

payload = {
    'scope': scope,
    'assembly': stat(assembly),
    'map': stat(map_file),
    'prediction_fingerprint': prediction_fp,
    'scope_selector': stat(scope_selector),
    'map_reader': stat(map_reader),
    'method': 'partition_full_gene_prediction_by_mag_map_v1',
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_SCOPE_FP
)"
    SCOPE_CACHE_ENTRY="$PRED_CACHE_ROOT/scopes/$SCOPE/$SCOPE_SELECTION_FP"
    SCOPE_CACHE_FAA="$SCOPE_CACHE_ENTRY/genes.${SCOPE}.faa"
    SCOPE_CACHE_COORD="$SCOPE_CACHE_ENTRY/gene_coordinates.${SCOPE}.tsv"
    SCOPE_CACHE_META="$SCOPE_CACHE_ENTRY/run_metadata.tsv"
    SCOPE_CACHE_FP="$SCOPE_CACHE_ENTRY/inputs.sha256"
fi

INPUT_FP="$(python3 - \
    "$SCOPE" "$ASSEMBLY_FA" "$MAG_CONTIG_MAP" "$MMSEQS_TAX_DB" \
    "$PRODIGAL_LINKER" "$FASTA_SELECTOR" "$GENE_SCOPE_SELECTOR" \
    "$MMSEQS_TAXONOMY_MERGER" "$TAXONOMY_PREFILTER" \
    "$PREDICTOR_BIN" "$MMSEQS_BIN" "$GENE_PREDICTOR" "$0" \
    "$TAXONOMY_SOURCE" "$TAXONOMY_RELEASE" "$BLACKLIST_POLICY" \
    "$PREDICTION_FP" "$SCOPE_SELECTION_FP" "$TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG" <<'PY_FP'
import glob
import hashlib
import json
import os
import sys

(
    scope, assembly, map_file, db, linker, selector, scope_selector,
    taxonomy_merger, taxonomy_prefilter,
    predictor_bin, mmseqs_bin, predictor, core_script,
    taxonomy_source, taxonomy_release, blacklist_policy,
    prediction_fingerprint, scope_selection_fingerprint, min_predicted_genes,
) = sys.argv[1:]

def stat(path):
    value = os.stat(path)
    return [os.path.realpath(path), value.st_size, value.st_mtime_ns]

db_files = sorted(path for path in glob.glob(db + '*') if os.path.isfile(path))
if not db_files:
    raise SystemExit(f'ERROR: no MMseqs2 database files found for prefix: {db}')

payload = {
    'scope': scope,
    'assembly': stat(assembly),
    'map': None if map_file in {'NA', 'None', ''} else stat(map_file),
    'database_files': [stat(path) for path in db_files],
    'linker': stat(linker),
    'selector': stat(selector),
    'scope_selector': stat(scope_selector),
    'taxonomy_merger': stat(taxonomy_merger),
    'taxonomy_prefilter': stat(taxonomy_prefilter),
    'predictor_binary': stat(predictor_bin),
    'mmseqs_binary': stat(mmseqs_bin),
    'predictor': predictor,
    'core_script': stat(core_script),
    'taxonomy_source': taxonomy_source,
    'taxonomy_release': taxonomy_release,
    'blacklist_policy': blacklist_policy,
    'gene_prediction_fingerprint': prediction_fingerprint,
    'scope_selection_fingerprint': scope_selection_fingerprint,
    'taxonomy_min_predicted_genes_per_contig': int(min_predicted_genes),
    'method': (
        'full_assembly_meta_gene_prediction|post_prediction_scope_partition|'
        'mmseqs_taxonomy_lca_mode_3|'
        'lca_ranks_superkingdom_phylum_class_order_family|tax_lineage_2_domain_recovery|'
        'gtdb_empty_blacklist|blank_rank_unclassified|disk_backed_gene_id_merge|schema_rank_names_v2|'
        'lossless_min_predicted_gene_prefilter'
    ),
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_FP
)"

cached_output_is_valid() {
    python3 - "$OUT_TSV" "$RUN_META" <<'PY_VALIDATE_CACHE'
import csv
import sys
from pathlib import Path

out_path, meta_path = map(Path, sys.argv[1:])
expected = ['contig', 'start', 'end', 'domain', 'phylum', 'class', 'order', 'family', 'gene_id']
try:
    with meta_path.open(encoding='utf-8-sig', newline='') as handle:
        metadata = {row['key']: row['value'] for row in csv.DictReader(handle, delimiter='\t')}
    selected = int(metadata['selected_contigs'])
    genes = int(metadata['predicted_genes'])
    with out_path.open(encoding='utf-8-sig') as handle:
        header = handle.readline().rstrip('\n\r').split('\t')
        first = handle.readline().rstrip('\n\r')
except (OSError, KeyError, ValueError, csv.Error):
    raise SystemExit(1)
if header != expected:
    raise SystemExit(1)
if selected < 0 or genes < 0:
    raise SystemExit(1)
if genes == 0:
    if first:
        raise SystemExit(1)
else:
    if selected < 1 or not first or len(first.split('\t')) != len(expected):
        raise SystemExit(1)
PY_VALIDATE_CACHE
}

prediction_cache_is_valid() {
    [[ -s "$PRED_CACHE_FAA" && -s "$PRED_CACHE_COORD" \
       && -s "$PRED_CACHE_META" && -s "$PRED_CACHE_FP" \
       && "$(<"$PRED_CACHE_FP")" == "$PREDICTION_FP" ]] || return 1
    python3 - \
        "$PRED_CACHE_FAA" "$PRED_CACHE_COORD" "$PRED_CACHE_META" \
        "$PREDICTION_FP" <<'PY_VALIDATE_PREDICTION_CACHE'
import csv
import sys
from pathlib import Path

faa_path, coord_path, meta_path = map(Path, sys.argv[1:4])
expected_fingerprint = sys.argv[4]
try:
    with meta_path.open(encoding='utf-8-sig', newline='') as handle:
        metadata = {
            row['key']: row['value']
            for row in csv.DictReader(handle, delimiter='\t')
        }
    selected = int(metadata['selected_contigs'])
    expected_genes = int(metadata['predicted_genes'])
    if metadata['input_fingerprint'] != expected_fingerprint:
        raise ValueError('fingerprint mismatch')
    if metadata['cache_schema'] != 'full_assembly_gene_prediction_cache_v1':
        raise ValueError('cache schema mismatch')
    if selected < 1 or expected_genes < 1:
        raise ValueError('invalid cached counts')

    faa_genes = 0
    with faa_path.open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                faa_genes += 1

    coordinate_genes = 0
    with coord_path.open(encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) != 4:
                raise ValueError('invalid coordinate row')
            int(fields[2])
            int(fields[3])
            coordinate_genes += 1
except (OSError, KeyError, ValueError, csv.Error):
    raise SystemExit(1)

if faa_genes != expected_genes or coordinate_genes != expected_genes:
    raise SystemExit(1)
PY_VALIDATE_PREDICTION_CACHE
}

scope_cache_is_valid() {
    [[ "$SCOPE" != "all" ]] || return 1
    [[ -e "$SCOPE_CACHE_FAA" && -e "$SCOPE_CACHE_COORD" \
       && -s "$SCOPE_CACHE_META" && -s "$SCOPE_CACHE_FP" \
       && "$(<"$SCOPE_CACHE_FP")" == "$SCOPE_SELECTION_FP" ]] || return 1
    python3 - \
        "$SCOPE_CACHE_FAA" "$SCOPE_CACHE_COORD" "$SCOPE_CACHE_META" \
        "$SCOPE_SELECTION_FP" "$SCOPE" <<'PY_VALIDATE_SCOPE_CACHE'
import csv
import sys
from pathlib import Path

faa_path, coord_path, meta_path = map(Path, sys.argv[1:4])
expected_fingerprint, expected_scope = sys.argv[4:]
try:
    with meta_path.open(encoding='utf-8-sig', newline='') as handle:
        metadata = {
            row['key']: row['value']
            for row in csv.DictReader(handle, delimiter='\t')
        }
    selected = int(metadata['selected_contigs'])
    expected_genes = int(metadata['predicted_genes'])
    if metadata['input_fingerprint'] != expected_fingerprint:
        raise ValueError('fingerprint mismatch')
    if metadata['cache_schema'] != 'scope_gene_prediction_cache_v1':
        raise ValueError('cache schema mismatch')
    if metadata['scope'] != expected_scope or selected < 0 or expected_genes < 0:
        raise ValueError('scope/count mismatch')

    faa_genes = sum(
        1 for line in faa_path.open(encoding='utf-8') if line.startswith('>')
    )
    coordinate_genes = 0
    with coord_path.open(encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) != 4:
                raise ValueError('invalid coordinate row')
            int(fields[2])
            int(fields[3])
            coordinate_genes += 1
except (OSError, KeyError, ValueError, csv.Error):
    raise SystemExit(1)

if faa_genes != expected_genes or coordinate_genes != expected_genes:
    raise SystemExit(1)
PY_VALIDATE_SCOPE_CACHE
}

write_run_metadata() {
    {
        printf 'key\tvalue\n'
        printf 'script_id\t%s\n' "$SCRIPT_ID"
        printf 'scope\t%s\n' "$SCOPE"
        printf 'gene_predictor\t%s\n' "$GENE_PREDICTOR"
        printf 'full_assembly_contigs\t%s\n' "$N_ALL_CONTIGS"
        printf 'selected_contigs\t%s\n' "$N_SCOPE_CONTIGS"
        printf 'full_assembly_predicted_genes\t%s\n' "$N_ALL_PREDICTED_GENES"
        printf 'predicted_genes\t%s\n' "$N_PREDICTED_GENES"
        printf 'predictor_realpath\t%s\n' "$PREDICTOR_REALPATH"
        printf 'predictor_signature\t%s\n' "$PREDICTOR_SIGNATURE"
        printf 'gene_prediction_threads\t%s\n' "$GENE_PRED_THREADS"
        printf 'pyrodigal_pool\t%s\n' "$PYRODIGAL_POOL"
        printf 'prodigal_jobs\t%s\n' "$PRODIGAL_JOBS"
        printf 'prodigal_chunk_bp\t%s\n' "$PRODIGAL_CHUNK_BP"
        printf 'prediction_method_version\t%s\n' "$GENE_PREDICTION_METHOD_VERSION"
        printf 'gene_prediction_fingerprint\t%s\n' "$PREDICTION_FP"
        printf 'gene_prediction_cache\t%s\n' "$PRED_CACHE_ENTRY"
        printf 'gene_prediction_cache_status\t%s\n' "$PREDICTION_CACHE_STATUS"
        printf 'scope_selection_fingerprint\t%s\n' "$SCOPE_SELECTION_FP"
        printf 'scope_gene_cache\t%s\n' "$SCOPE_CACHE_ENTRY"
        printf 'scope_gene_cache_status\t%s\n' "$SCOPE_CACHE_STATUS"
        printf 'mmseqs_threads\t%s\n' "$MMSEQS_THREADS"
        printf 'taxonomy_min_predicted_genes_per_contig\t%s\n' "$TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG"
        printf 'taxonomy_query_contigs\t%s\n' "${N_TAX_QUERY_CONTIGS:-$N_SCOPE_CONTIGS}"
        printf 'taxonomy_query_genes\t%s\n' "${N_TAX_QUERY_GENES:-$N_PREDICTED_GENES}"
        printf 'taxonomy_prefilter_skipped_contigs\t%s\n' "${N_TAX_SKIPPED_CONTIGS:-0}"
        printf 'taxonomy_prefilter_skipped_genes\t%s\n' "${N_TAX_SKIPPED_GENES:-0}"
        printf 'mmseqs_database\t%s\n' "$MMSEQS_TAX_DB"
        printf 'taxonomy_source\t%s\n' "$TAXONOMY_SOURCE"
        printf 'taxonomy_release\t%s\n' "$TAXONOMY_RELEASE"
        printf 'taxonomy_db_label\t%s\n' "$TAXONOMY_DB_LABEL"
        printf 'gtdb_db_label\t%s\n' "$GTDB_DB_LABEL"
        printf 'blacklist_policy\t%s\n' "$BLACKLIST_POLICY"
        printf 'blacklist_value\t%s\n' "$BLACKLIST_VALUE"
        printf 'output_schema\tnine_column_rank_names_v1\n'
        printf 'lca_mode\t3\n'
        printf 'lca_ranks\tsuperkingdom,phylum,class,order,family\n'
        printf 'tax_lineage_mode\t2\n'
        printf 'input_fingerprint\t%s\n' "$INPUT_FP"
    } > "$TMP_RUN_META"
    mv -f "$TMP_RUN_META" "$RUN_META"
}

if [[ -s "$OUT_TSV" && -s "$CACHE_META" && -s "$RUN_META" \
      && "$(<"$CACHE_META")" == "$INPUT_FP" ]]; then
    if cached_output_is_valid; then
        echo "Fingerprint-matched taxonomy table already exists: $OUT_TSV"
        exit 0
    fi
    echo "WARNING: cached taxonomy output or metadata is incomplete; rebuilding: $OUT_TSV" >&2
fi

WORK="$(mktemp -d "$OUT_DIR/.mmseqs2_gene_taxonomy.XXXXXX")"
PREDICTION_CACHE_STATUS="created"

if prediction_cache_is_valid; then
    PRED_CACHE_COUNTS="$(python3 - "$PRED_CACHE_META" <<'PY_READ_PREDICTION_CACHE'
import csv
import sys
with open(sys.argv[1], encoding='utf-8-sig', newline='') as handle:
    metadata = {
        row['key']: row['value']
        for row in csv.DictReader(handle, delimiter='\t')
    }
print(metadata['selected_contigs'], metadata['predicted_genes'])
PY_READ_PREDICTION_CACHE
)"
    read -r N_ALL_CONTIGS N_ALL_PREDICTED_GENES <<< "$PRED_CACHE_COUNTS"
    N_PREDICTED_GENES="$N_ALL_PREDICTED_GENES"
    GENES_FAA="$PRED_CACHE_FAA"
    GENE_COORD="$PRED_CACHE_COORD"
    PREDICTION_CACHE_STATUS="reused"
else
    SCOPE_FASTA="$ASSEMBLY_FA"
    ALL_STATS="$WORK/all.stats.json"
    python3 - "$ASSEMBLY_FA" "$ALL_STATS" <<'PY_ALL_STATS'
import json
import sys
fasta, out = sys.argv[1:]
n = bp = 0
with open(fasta) as handle:
    for line in handle:
        if line.startswith('>'):
            n += 1
        else:
            bp += len(line.strip())
with open(out, 'w') as dst:
    json.dump({'mode': 'all', 'output_contigs': n, 'output_bp': bp}, dst, indent=2)
    dst.write('\n')
print(f'All-contig scope: {n} contigs; {bp} bp')
PY_ALL_STATS
    N_ALL_CONTIGS="$(python3 - "$ALL_STATS" <<'PY_COUNT'
import json, sys
print(json.load(open(sys.argv[1]))['output_contigs'])
PY_COUNT
)"
    (( N_ALL_CONTIGS > 0 )) || {
        echo "ERROR: assembly contains no FASTA records: $ASSEMBLY_FA" >&2
        exit 1
    }
    GENES_FAA="$WORK/genes.all.faa"
    GENES_GFF="$WORK/genes.all.gff"
    GENE_COORD="$WORK/gene_coordinates.all.tsv"
fi

QUERY_DB="$WORK/queryDB"
TAX_RESULT="$WORK/taxonomyResult"
TAX_RAW="$WORK/taxonomyResult.tsv"
MMSEQS_TMP="$WORK/tmp"
mkdir -p "$MMSEQS_TMP"

MMSEQS_VERSION="$("$MMSEQS_BIN" version 2>&1)"
echo "[00_gene_taxonomy_mmseqs2] Started: $(date)"
echo "SCRIPT_ID        : $SCRIPT_ID"
echo "SCOPE                 : $SCOPE"
echo "FULL_ASSEMBLY_CONTIGS : $N_ALL_CONTIGS"
echo "ASSEMBLY_FA           : $ASSEMBLY_FA"
echo "MAG_CONTIG_MAP        : $MAG_CONTIG_MAP"
echo "GENE_PREDICTOR        : $GENE_PREDICTOR"
echo "GENE_PRED_THREADS     : $GENE_PRED_THREADS"
echo "PYRODIGAL_POOL        : $PYRODIGAL_POOL"
echo "PRODIGAL_JOBS         : $PRODIGAL_JOBS"
echo "PRODIGAL_CHUNK_BP     : $PRODIGAL_CHUNK_BP"
echo "GENE_PRED_CACHE       : $PRED_CACHE_ENTRY"
echo "GENE_PRED_CACHE_STATUS: $PREDICTION_CACHE_STATUS"
echo "MMSEQS_BIN            : $MMSEQS_BIN"
echo "MMSEQS_TAX_DB         : $MMSEQS_TAX_DB"
echo "TAXONOMY_SOURCE       : $TAXONOMY_SOURCE"
echo "TAXONOMY_RELEASE      : $TAXONOMY_RELEASE"
echo "TAXONOMY_DB_LABEL     : $TAXONOMY_DB_LABEL"
echo "GTDB_DB_LABEL         : $GTDB_DB_LABEL"
echo "MMSEQS_BLACKLIST      : $BLACKLIST_POLICY ($BLACKLIST_VALUE)"
echo "MMSEQS_DB_TYPE        : $DB_TYPE"
echo "MMSEQS_VERSION        : ${MMSEQS_VERSION%%$'\n'*}"
echo "MMSEQS_THREADS        : $MMSEQS_THREADS"
echo "TAX_MIN_GENES/CONTIG  : $TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG"
echo "KEEP_FAILED_WORK      : $KEEP_FAILED_WORK"
echo "OUT_TSV               : $OUT_TSV"

if [[ "$PREDICTION_CACHE_STATUS" == "reused" ]]; then
    echo "[1/5] Reusing validated full-assembly gene predictions"
else
echo "[1/5] Predicting genes across the full assembly and creating the shared cache"
if [[ "$GENE_PREDICTOR" == "pyrodigal" ]]; then
    "$PYRODIGAL_BIN" \
        -p meta \
        -j "$GENE_PRED_THREADS" \
        --pool "$PYRODIGAL_POOL" \
        -i "$SCOPE_FASTA" \
        -a "$GENES_FAA" \
        -f gff \
        -o "$GENES_GFF"
elif [[ "$GENE_PREDICTOR" == "prodigal_serial" ]]; then
    "$PRODIGAL_BIN" -q -p meta -i "$SCOPE_FASTA" \
        -a "$GENES_FAA" -f gff -o "$GENES_GFF"
else
    CHUNK_DIR="$WORK/prodigal_chunks"
    RESULT_DIR="$WORK/prodigal_results"
    MANIFEST="$WORK/prodigal_chunks.list"
    mkdir -p "$CHUNK_DIR" "$RESULT_DIR"

    python3 - "$SCOPE_FASTA" "$CHUNK_DIR" "$MANIFEST" "$PRODIGAL_CHUNK_BP" <<'PY_SPLIT'
import sys
from pathlib import Path

fasta, out_dir, manifest, target_bp = sys.argv[1:]
out_dir = Path(out_dir)
target_bp = int(target_bp)
chunk_index = 0
chunk_bp = 0
chunk_handle = None
paths = []

def new_chunk():
    global chunk_index, chunk_bp, chunk_handle
    if chunk_handle is not None:
        chunk_handle.close()
    chunk_index += 1
    path = out_dir / f'chunk_{chunk_index:06d}.fa'
    paths.append(path)
    chunk_handle = path.open('w')
    chunk_bp = 0

def emit(header, sequence_lines):
    global chunk_bp
    length = sum(len(line.strip()) for line in sequence_lines)
    if chunk_handle is None or (chunk_bp > 0 and chunk_bp + length > target_bp):
        new_chunk()
    chunk_handle.write(header)
    if not header.endswith('\n'):
        chunk_handle.write('\n')
    for line in sequence_lines:
        chunk_handle.write(line.rstrip('\n\r') + '\n')
    chunk_bp += length

header = None
sequence = []
with open(fasta) as handle:
    for line in handle:
        if line.startswith('>'):
            if header is not None:
                emit(header, sequence)
            header = line
            sequence = []
        else:
            sequence.append(line)
if header is not None:
    emit(header, sequence)
if chunk_handle is not None:
    chunk_handle.close()
with open(manifest, 'w') as out:
    for path in paths:
        out.write(str(path) + '\n')
print(f'Prodigal chunks: {len(paths)}')
PY_SPLIT

    WORKER="$WORK/run_prodigal_chunk.sh"
    cat > "$WORKER" <<'WORKER_SH'
#!/usr/bin/env bash
set -Eeuo pipefail
prodigal_bin="$1"
result_dir="$2"
chunk="$3"
base="$(basename "$chunk" .fa)"
"$prodigal_bin" -q -p meta -i "$chunk" \
    -a "$result_dir/$base.faa" -f gff -o "$result_dir/$base.gff"
WORKER_SH
    chmod +x "$WORKER"
    xargs -r -n 1 -P "$PRODIGAL_JOBS" \
        "$WORKER" "$PRODIGAL_BIN" "$RESULT_DIR" < "$MANIFEST"

    : > "$GENES_FAA"
    printf '%s\n' '##gff-version 3' > "$GENES_GFF"
    while IFS= read -r chunk; do
        base="$(basename "$chunk" .fa)"
        [[ -e "$RESULT_DIR/$base.faa" && -e "$RESULT_DIR/$base.gff" ]] || {
            echo "ERROR: missing Prodigal chunk output: $base" >&2
            exit 1
        }
        cat "$RESULT_DIR/$base.faa" >> "$GENES_FAA"
        grep -v '^#' "$RESULT_DIR/$base.gff" >> "$GENES_GFF" || true
    done < "$MANIFEST"
fi

[[ -s "$GENES_FAA" && -s "$GENES_GFF" ]] || {
    echo "ERROR: gene predictor produced no protein/GFF output" >&2
    exit 1
}

echo "[1/5] Linking full-assembly protein accessions to exact gene coordinates"
python3 "$PRODIGAL_LINKER" "$GENES_FAA" "$GENES_GFF" "$GENE_COORD"
[[ -s "$GENE_COORD" ]] || {
    echo "ERROR: failed to link protein accessions to GFF coordinates" >&2
    exit 1
}
N_PREDICTED_GENES="$(wc -l < "$GENE_COORD")"
[[ "$N_PREDICTED_GENES" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: invalid predicted-gene count from $GENE_COORD: $N_PREDICTED_GENES" >&2
    exit 1
}
N_ALL_PREDICTED_GENES="$N_PREDICTED_GENES"

mkdir -p "$PRED_CACHE_ENTRY"
PRED_CACHE_TMP_FAA="$PRED_CACHE_ENTRY/genes.all.faa.tmp.$$"
PRED_CACHE_TMP_COORD="$PRED_CACHE_ENTRY/gene_coordinates.all.tsv.tmp.$$"
PRED_CACHE_TMP_META="$PRED_CACHE_ENTRY/run_metadata.tsv.tmp.$$"
PRED_CACHE_TMP_FP="$PRED_CACHE_ENTRY/inputs.sha256.tmp.$$"

mv -f "$GENES_FAA" "$PRED_CACHE_TMP_FAA"
mv -f "$PRED_CACHE_TMP_FAA" "$PRED_CACHE_FAA"
mv -f "$GENE_COORD" "$PRED_CACHE_TMP_COORD"
mv -f "$PRED_CACHE_TMP_COORD" "$PRED_CACHE_COORD"
{
    printf 'key\tvalue\n'
    printf 'cache_schema\tfull_assembly_gene_prediction_cache_v1\n'
    printf 'scope\tall\n'
    printf 'gene_predictor\t%s\n' "$GENE_PREDICTOR"
    printf 'selected_contigs\t%s\n' "$N_ALL_CONTIGS"
    printf 'predicted_genes\t%s\n' "$N_ALL_PREDICTED_GENES"
    printf 'predictor_realpath\t%s\n' "$PREDICTOR_REALPATH"
    printf 'predictor_signature\t%s\n' "$PREDICTOR_SIGNATURE"
    printf 'gene_prediction_threads\t%s\n' "$GENE_PRED_THREADS"
    printf 'pyrodigal_pool\t%s\n' "$PYRODIGAL_POOL"
    printf 'prodigal_jobs\t%s\n' "$PRODIGAL_JOBS"
    printf 'prodigal_chunk_bp\t%s\n' "$PRODIGAL_CHUNK_BP"
    printf 'prediction_method_version\t%s\n' "$GENE_PREDICTION_METHOD_VERSION"
    printf 'input_fingerprint\t%s\n' "$PREDICTION_FP"
} > "$PRED_CACHE_TMP_META"
mv -f "$PRED_CACHE_TMP_META" "$PRED_CACHE_META"
printf '%s\n' "$PREDICTION_FP" > "$PRED_CACHE_TMP_FP"
mv -f "$PRED_CACHE_TMP_FP" "$PRED_CACHE_FP"

GENES_FAA="$PRED_CACHE_FAA"
GENE_COORD="$PRED_CACHE_COORD"
prediction_cache_is_valid || {
    echo "ERROR: newly written gene-prediction cache failed validation: $PRED_CACHE_ENTRY" >&2
    exit 1
}
fi

if [[ "$SCOPE" == "all" ]]; then
    N_SCOPE_CONTIGS="$N_ALL_CONTIGS"
    N_PREDICTED_GENES="$N_ALL_PREDICTED_GENES"
    GENES_FAA="$PRED_CACHE_FAA"
    GENE_COORD="$PRED_CACHE_COORD"
    SCOPE_CACHE_STATUS="full_assembly_cache"
    echo "[2/5] Using all full-assembly gene predictions"
elif scope_cache_is_valid; then
    SCOPE_CACHE_COUNTS="$(python3 - "$SCOPE_CACHE_META" <<'PY_READ_SCOPE_CACHE'
import csv
import sys
with open(sys.argv[1], encoding='utf-8-sig', newline='') as handle:
    metadata = {
        row['key']: row['value']
        for row in csv.DictReader(handle, delimiter='\t')
    }
print(metadata['selected_contigs'], metadata['predicted_genes'])
PY_READ_SCOPE_CACHE
)"
    read -r N_SCOPE_CONTIGS N_PREDICTED_GENES <<< "$SCOPE_CACHE_COUNTS"
    GENES_FAA="$SCOPE_CACHE_FAA"
    GENE_COORD="$SCOPE_CACHE_COORD"
    SCOPE_CACHE_STATUS="reused"
    echo "[2/5] Reusing validated $SCOPE gene-selection cache"
else
    SCOPE_MODE="include"
    [[ "$SCOPE" != "nonmag" ]] || SCOPE_MODE="exclude"
    SCOPE_WORK_FAA="$WORK/genes.${SCOPE}.faa"
    SCOPE_WORK_COORD="$WORK/gene_coordinates.${SCOPE}.tsv"
    SCOPE_WORK_STATS="$WORK/${SCOPE}.stats.json"
    echo "[2/5] Selecting $SCOPE genes from the full-assembly prediction cache"
    python3 "$GENE_SCOPE_SELECTOR" genes "$SCOPE_MODE" \
        "$ASSEMBLY_FA" "$MAG_CONTIG_MAP" \
        "$PRED_CACHE_FAA" "$PRED_CACHE_COORD" \
        "$SCOPE_WORK_FAA" "$SCOPE_WORK_COORD" "$SCOPE_WORK_STATS"
    SCOPE_COUNTS="$(python3 - "$SCOPE_WORK_STATS" <<'PY_SCOPE_COUNTS'
import json
import sys
stats = json.load(open(sys.argv[1]))
print(stats['output_contigs'], stats['output_genes'])
PY_SCOPE_COUNTS
)"
    read -r N_SCOPE_CONTIGS N_PREDICTED_GENES <<< "$SCOPE_COUNTS"

    mkdir -p "$SCOPE_CACHE_ENTRY"
    SCOPE_CACHE_TMP_FAA="$SCOPE_CACHE_ENTRY/genes.${SCOPE}.faa.tmp.$$"
    SCOPE_CACHE_TMP_COORD="$SCOPE_CACHE_ENTRY/gene_coordinates.${SCOPE}.tsv.tmp.$$"
    SCOPE_CACHE_TMP_META="$SCOPE_CACHE_ENTRY/run_metadata.tsv.tmp.$$"
    SCOPE_CACHE_TMP_FP="$SCOPE_CACHE_ENTRY/inputs.sha256.tmp.$$"
    mv -f "$SCOPE_WORK_FAA" "$SCOPE_CACHE_TMP_FAA"
    mv -f "$SCOPE_CACHE_TMP_FAA" "$SCOPE_CACHE_FAA"
    mv -f "$SCOPE_WORK_COORD" "$SCOPE_CACHE_TMP_COORD"
    mv -f "$SCOPE_CACHE_TMP_COORD" "$SCOPE_CACHE_COORD"
    {
        printf 'key\tvalue\n'
        printf 'cache_schema\tscope_gene_prediction_cache_v1\n'
        printf 'scope\t%s\n' "$SCOPE"
        printf 'full_prediction_fingerprint\t%s\n' "$PREDICTION_FP"
        printf 'selected_contigs\t%s\n' "$N_SCOPE_CONTIGS"
        printf 'predicted_genes\t%s\n' "$N_PREDICTED_GENES"
        printf 'input_fingerprint\t%s\n' "$SCOPE_SELECTION_FP"
    } > "$SCOPE_CACHE_TMP_META"
    mv -f "$SCOPE_CACHE_TMP_META" "$SCOPE_CACHE_META"
    printf '%s\n' "$SCOPE_SELECTION_FP" > "$SCOPE_CACHE_TMP_FP"
    mv -f "$SCOPE_CACHE_TMP_FP" "$SCOPE_CACHE_FP"
    GENES_FAA="$SCOPE_CACHE_FAA"
    GENE_COORD="$SCOPE_CACHE_COORD"
    SCOPE_CACHE_STATUS="created"
    scope_cache_is_valid || {
        echo "ERROR: newly written scope cache failed validation: $SCOPE_CACHE_ENTRY" >&2
        exit 1
    }
fi

echo "SCOPE_CONTIGS         : $N_SCOPE_CONTIGS"
echo "FULL_PREDICTED_GENES  : $N_ALL_PREDICTED_GENES"
echo "PREDICTED_GENES       : $N_PREDICTED_GENES"
echo "SCOPE_GENE_CACHE      : $SCOPE_CACHE_ENTRY"
echo "SCOPE_CACHE_STATUS    : $SCOPE_CACHE_STATUS"

if (( N_PREDICTED_GENES == 0 )); then
    printf '%s\n' $'contig\tstart\tend\tdomain\tphylum\tclass\torder\tfamily\tgene_id' > "$TMP_OUT"
    mv -f "$TMP_OUT" "$OUT_TSV"
    printf '%s\n' "$INPUT_FP" > "$TMP_CACHE"
    mv -f "$TMP_CACHE" "$CACHE_META"
    write_run_metadata
    echo "No genes were selected for scope=$SCOPE; wrote a header-only table: $OUT_TSV"
    exit 0
fi

TAX_QUERY_FAA="$GENES_FAA"
N_TAX_QUERY_CONTIGS="$N_SCOPE_CONTIGS"
N_TAX_QUERY_GENES="$N_PREDICTED_GENES"
N_TAX_SKIPPED_CONTIGS=0
N_TAX_SKIPPED_GENES=0
ALLOW_MISSING_TAXONOMY=N
if (( TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG > 1 )); then
    echo "[3/6] Prefiltering genes from contigs that cannot meet the taxonomic-transition gene-count criterion"
    TAX_QUERY_FAA="$WORK/genes.taxonomy_eligible.faa"
    TAX_FILTER_STATS="$WORK/taxonomy_prefilter.stats.json"
    python3 "$TAXONOMY_PREFILTER" \
        "$GENES_FAA" "$GENE_COORD" "$TAX_QUERY_FAA" "$TAX_FILTER_STATS" \
        "$TAXONOMY_MIN_PREDICTED_GENES_PER_CONTIG"
    read -r N_TAX_QUERY_CONTIGS N_TAX_QUERY_GENES N_TAX_SKIPPED_CONTIGS N_TAX_SKIPPED_GENES < <(
        python3 - "$TAX_FILTER_STATS" <<'PY_PREFILTER_STATS'
import json, sys
s=json.load(open(sys.argv[1]))
print(s['selected_contigs'], s['selected_genes'], s['skipped_contigs'], s['skipped_genes'])
PY_PREFILTER_STATS
    )
    ALLOW_MISSING_TAXONOMY=Y
else
    echo "[3/6] Taxonomy prefilter disabled (minimum predicted genes per contig = 1)"
fi

echo "TAXONOMY_QUERY_CONTIGS : $N_TAX_QUERY_CONTIGS / $N_SCOPE_CONTIGS"
echo "TAXONOMY_QUERY_GENES   : $N_TAX_QUERY_GENES / $N_PREDICTED_GENES"

: > "$TAX_RAW"
if (( N_TAX_QUERY_GENES > 0 )); then
    echo "[4/6] Creating MMseqs2 query database"
    "$MMSEQS_BIN" createdb "$TAX_QUERY_FAA" "$QUERY_DB"

    echo "[5/6] Assigning gene taxonomy with approximate 2bLCA"
    MMSEQS_TAXONOMY_ARGS=(
        --threads "$MMSEQS_THREADS"
        --lca-mode 3
        --lca-ranks superkingdom,phylum,class,order,family
        --tax-lineage 2
    )
    if [[ "$APPLY_EMPTY_BLACKLIST" == "Y" ]]; then
        MMSEQS_TAXONOMY_ARGS+=(--blacklist "")
    fi
    "$MMSEQS_BIN" taxonomy \
        "$QUERY_DB" "$MMSEQS_TAX_DB" "$TAX_RESULT" "$MMSEQS_TMP" \
        "${MMSEQS_TAXONOMY_ARGS[@]}"
    "$MMSEQS_BIN" createtsv "$QUERY_DB" "$TAX_RESULT" "$TAX_RAW" \
        --threads "$MMSEQS_THREADS"
    [[ -s "$TAX_RAW" ]] || {
        echo "ERROR: MMseqs2 produced an empty taxonomy TSV" >&2
        exit 1
    }
else
    echo "[4/6] No contig meets the taxonomy prefilter; MMseqs2 taxonomy is skipped"
fi

echo "[6/6] Producing the nine-column downstream table"
MERGE_ARGS=("$GENE_COORD" "$TAX_RAW" "$TMP_OUT")
if [[ "$ALLOW_MISSING_TAXONOMY" == "Y" ]]; then
    MERGE_ARGS+=(--allow-missing-as-unclassified)
fi
python3 "$MMSEQS_TAXONOMY_MERGER" "${MERGE_ARGS[@]}"

mv -f "$TMP_OUT" "$OUT_TSV"
printf '%s\n' "$INPUT_FP" > "$TMP_CACHE"
mv -f "$TMP_CACHE" "$CACHE_META"
write_run_metadata

echo "Output: $OUT_TSV"
echo "Metadata: $RUN_META"
echo "[gene_taxonomy] Finished: $(date)"
exit 0
fi

set -Eeuo pipefail

SCRIPT_ID="07_chimera"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

usage() {
    cat >&2 <<'USAGE'
Usage:
  bash 07_chimera.sh [GENE_TAX_TSV OUT_DIR]

Gene taxonomy is generated automatically if missing and MMSEQS_TAX_DB is provided.
SCG annotations are generated when anvi'o HMM tools are available.

Structural-chimera policy:
  - Scope is all contigs in the three long-read assemblies. Inter-contig MAG taxonomy is not evaluated.
  - GTDB domain/phylum/class/order transitions are candidate junctions, not final chimera calls.
  - Local taxonomy uses a bounded adaptive gene search: <=10 predicted genes searched per side,
    <=5 classified genes used, >=4 classified genes required and >=80% dominance by default.
  - Existing 03_clipping.sh results are reused. Clipping is never recomputed here.
  - Final taxonomic chimera requires a taxonomic-shift region containing a validated
    Clipping=Y site from 03_clipping.sh.
  - High-confidence taxonomic chimera additionally requires Stringent_clipping=Y
    (validated Clipping=Y plus >=10x anvi'o coverage and 100% clipping).
  - Primary + supplementary HiFi structural validation is retained as an audit annotation only.
  - Secondary alignments are excluded from that structural audit.
  - HiFi continuity requires >=2 unique spanning QNAMEs with >=500 bp matched anchors on both sides.
  - The structural audit records breakpoint support and unresolved regions but does not determine
    the final Taxonomic_chimera call.
  - GC and HiFi depth are descriptive annotations for all assembly contigs.
  - SCG redundancy is annotated only for MAG-associated contigs using the MAG GTDB domain-matched anvi'o collection (Bacteria_71 or Archaea_76) and never creates a structural call.

The final contig table is assembly-wide. Set CHIMERA_SUMMARY_ONLY=Y to regenerate the MAG-associated subset and MAG-level tables from it.
USAGE
    exit 2
}

[[ $# -le 2 ]] || usage
GENE_TAX_OVERRIDE="${1:-}"

source "$SCRIPT_DIR/lib_config.sh"
source "$SCRIPT_DIR/lib_qc.sh"

TAX_SEARCH_MAX_GENES="${TAX_SEARCH_MAX_GENES:-10}"
TAX_TARGET_CLASSIFIED_GENES="${TAX_TARGET_CLASSIFIED_GENES:-5}"
TAX_MIN_CLASSIFIED_GENES="${TAX_MIN_CLASSIFIED_GENES:-4}"
TAX_MIN_BLOCK_DOMINANCE="${TAX_MIN_BLOCK_DOMINANCE:-0.80}"
CHIMERA_HIFI_ANCHOR_BP="${CHIMERA_HIFI_ANCHOR_BP:-500}"
CHIMERA_MIN_SUPPORT_QNAMES="${CHIMERA_MIN_SUPPORT_QNAMES:-2}"
CHIMERA_BAM_JOBS="${CHIMERA_BAM_JOBS:-3}"

python3 - "$TAX_SEARCH_MAX_GENES" "$TAX_TARGET_CLASSIFIED_GENES" \
    "$TAX_MIN_CLASSIFIED_GENES" "$TAX_MIN_BLOCK_DOMINANCE" \
    "$CHIMERA_HIFI_ANCHOR_BP" "$CHIMERA_MIN_SUPPORT_QNAMES" "$CHIMERA_BAM_JOBS" <<'PY'
import sys
search_max, target, minimum, dominance, anchor, min_q, bam_jobs = sys.argv[1:]
search_max, target, minimum, anchor, min_q, bam_jobs = map(int, (search_max, target, minimum, anchor, min_q, bam_jobs))
dominance = float(dominance)
if search_max < 1: raise SystemExit('ERROR: TAX_SEARCH_MAX_GENES must be >=1')
if target < 1 or target > search_max: raise SystemExit('ERROR: TAX_TARGET_CLASSIFIED_GENES invalid')
if minimum < 1 or minimum > target: raise SystemExit('ERROR: TAX_MIN_CLASSIFIED_GENES invalid')
if not 0.5 <= dominance <= 1.0: raise SystemExit('ERROR: TAX_MIN_BLOCK_DOMINANCE must be in [0.5,1]')
if anchor < 1 or min_q < 1 or bam_jobs < 1: raise SystemExit('ERROR: HiFi anchor/support/BAM-job thresholds must be positive')
PY

GTDB_RELEASE="${GTDB_RELEASE:-GTDB_R232}"
CONFIGURED_GTDB_GENE_TAX_TSV="${GTDB_GENE_TAX_TSV:-}"
GENE_TAX_SELECTED="${GENE_TAX_OVERRIDE:-$CONFIGURED_GTDB_GENE_TAX_TSV}"
[[ -n "$GENE_TAX_SELECTED" ]] || { echo "ERROR: GTDB gene-taxonomy TSV is not configured" >&2; exit 1; }
if [[ ! -s "$GENE_TAX_SELECTED" ]]; then
    [[ -n "${MMSEQS_TAX_DB:-}" ]] || {
        echo "ERROR: gene taxonomy is missing. Specify --database MMSEQS_PREFIX for MAGverify chimera, or --taxonomy TSV." >&2
        exit 1
    }
    [[ -s "${MMSEQS_TAX_DB%.dbtype}.dbtype" ]] || {
        echo "ERROR: MMseqs2 taxonomy database prefix not found: $MMSEQS_TAX_DB" >&2
        exit 1
    }
    echo "Generating gene taxonomy for chimera analysis"
    bash "$SCRIPT_DIR/07_chimera.sh" __taxonomy_internal contig \
        "$ASSEMBLY_FA" "$MMSEQS_TAX_DB" "$MAG_CONTIG_MAP" \
        "$GTDB_GENE_TAX_DIR" "${THREADS:-100}"
    [[ -s "$GENE_TAX_SELECTED" ]] || {
        echo "ERROR: gene taxonomy was not generated: $GENE_TAX_SELECTED" >&2
        exit 1
    }
fi

# If the all-contig table was built with the lossless MMseqs2 prefilter, make
# sure its threshold remains compatible with the current chimera gene-count
# criterion. This prevents a later lowering of TAX_MIN_CLASSIFIED_GENES from
# silently making previously skipped short contigs theoretically evaluable.
TAXONOMY_PREFILTER_THRESHOLD="NA"
if [[ "$GENE_TAX_SELECTED" == *.all_contigs.tsv ]]; then
    NONMAG_TAX_META="${GENE_TAX_SELECTED%.all_contigs.tsv}.nonMAG_contigs.tsv.run_metadata.tsv"
    if [[ -s "$NONMAG_TAX_META" ]]; then
        TAXONOMY_PREFILTER_THRESHOLD="$(awk -F'\t' '$1=="taxonomy_min_predicted_genes_per_contig"{print $2}' "$NONMAG_TAX_META")"
        if [[ -n "$TAXONOMY_PREFILTER_THRESHOLD" && "$TAXONOMY_PREFILTER_THRESHOLD" != "NA" ]]; then
            [[ "$TAXONOMY_PREFILTER_THRESHOLD" =~ ^[1-9][0-9]*$ ]] || {
                echo "ERROR: invalid taxonomy prefilter threshold in $NONMAG_TAX_META" >&2
                exit 1
            }
            SAFE_PREFILTER_MAX=$(( 2 * TAX_MIN_CLASSIFIED_GENES ))
            if (( TAXONOMY_PREFILTER_THRESHOLD > SAFE_PREFILTER_MAX )); then
                echo "ERROR: all-contig taxonomy was prefiltered at >=${TAXONOMY_PREFILTER_THRESHOLD} genes/contig," >&2
                echo "but current TAX_MIN_CLASSIFIED_GENES=$TAX_MIN_CLASSIFIED_GENES permits a maximum lossless threshold of ${SAFE_PREFILTER_MAX}." >&2
                echo "Rerun MAGverify chimera with --database and a compatible taxonomy threshold." >&2
                exit 1
            fi
        fi
    fi
fi

CHIMERA_DIR="${GTDB_CHIMERA_DIR:-${CHIMERA_ROOT_DIR}/$GTDB_RELEASE}"
WORK_DIR="$CHIMERA_DIR/.candidate_work"
EVENTS_OUT="$CHIMERA_DIR/taxonomic_transition_structural_validation.${GTDB_RELEASE}.events.tsv"
STRUCTURAL_EVENTS_OUT="$CHIMERA_DIR/structural_chimera.${GTDB_RELEASE}.events.tsv"
HIGH_CONF_EVENTS_OUT="$CHIMERA_DIR/high_confidence_structural_chimera.${GTDB_RELEASE}.events.tsv"
UNRESOLVED_EVENTS_OUT="$CHIMERA_DIR/unresolved_taxonomic_transition.${GTDB_RELEASE}.events.tsv"
CONTIG_OUT="${GTDB_CHIMERA_CONTIG:-$CHIMERA_DIR/contig_structural_chimera.${GTDB_RELEASE}.final.tsv}"
STRUCTURAL_CONTIG_OUT="$CHIMERA_DIR/structural_chimera.${GTDB_RELEASE}.contig.tsv"
HIGH_CONF_CONTIG_OUT="$CHIMERA_DIR/high_confidence_structural_chimera.${GTDB_RELEASE}.contig.tsv"
UNRESOLVED_CONTIG_OUT="$CHIMERA_DIR/unresolved_taxonomic_transition.${GTDB_RELEASE}.contig.tsv"
SCG_AUDIT_OUT="$CHIMERA_DIR/SCG_redundancy.${GTDB_RELEASE}.contig.tsv"
MAG_CONTIG_OUT="${GTDB_CHIMERA_MAG_CONTIG:-$CHIMERA_DIR/MAG_contig_structural_chimera.${GTDB_RELEASE}.final.tsv}"
MAG_OUT="${GTDB_CHIMERA_MAG:-$CHIMERA_DIR/MAG_structural_chimera.${GTDB_RELEASE}.final.tsv}"
RUN_META_OUT="$CHIMERA_DIR/run_metadata.tsv"

CANDIDATE_EVALUATOR="$SCRIPT_DIR/../scripts/chimera_evaluate_taxonomic_chimera_whole_contig_GTDB.py"
STRUCTURAL_VALIDATOR="$SCRIPT_DIR/../scripts/chimera_validate_taxonomic_chimera_hifi.py"
MAG_SUMMARIZER="$SCRIPT_DIR/../scripts/chimera_summarize_taxonomic_chimera_by_mag_GTDB_whole.py"
FINAL_RECLASSIFIER="$SCRIPT_DIR/../scripts/chimera_reclassify_taxonomic_chimera_by_clipping.py"
FINAL_CHIMERA_DIR="${TAXONOMIC_CHIMERA_CLIPPING_OUT_DIR:-$CHIMERA_DIR/clipping_defined}"
SUMMARY_ONLY="${CHIMERA_SUMMARY_ONLY:-${GTDB_CHIMERA_SUMMARY_ONLY:-N}}"
case "$SUMMARY_ONLY" in Y|N) ;; *) echo "ERROR: CHIMERA_SUMMARY_ONLY must be Y or N" >&2; exit 2;; esac

LOG="$OUT/00_log/07_chimera_GTDB_R232.log"
mkdir -p "$(dirname "$LOG")" "$CHIMERA_DIR"
exec > >(tee -a "$LOG") 2>&1

write_mag_outputs() {
    [[ -s "$CONTIG_OUT" ]] || { echo "ERROR: final structural-chimera contig table missing: $CONTIG_OUT" >&2; return 1; }
    python3 "$MAG_SUMMARIZER" \
        --contig-table "$CONTIG_OUT" \
        --mag-map "$MAG_CONTIG_MAP" \
        --mag-contig-out "$MAG_CONTIG_OUT" \
        --mag-summary-out "$MAG_OUT"
}

write_final_taxonomic_outputs() {
    for f in "$EVENTS_OUT" "$CONTIG_OUT" "$MAG_CONTIG_OUT" "$MAG_OUT" "$FINAL_RECLASSIFIER"; do
        [[ -s "$f" ]] || { echo "ERROR: final taxonomic-chimera reclassification input missing: $f" >&2; return 1; }
    done
    python3 "$FINAL_RECLASSIFIER" \
        --regions-in "$EVENTS_OUT" \
        --clipping "$CLIP_FINAL" \
        --contigs-in "$CONTIG_OUT" \
        --mag-contigs-in "$MAG_CONTIG_OUT" \
        --mags-in "$MAG_OUT" \
        --out-dir "$FINAL_CHIMERA_DIR" \
        --release "$GTDB_RELEASE"
}

echo "[07_chimera] Started: $(date)"
echo "SCRIPT_ID                  : $SCRIPT_ID"
echo "GENE_TAX_TSV                   : $GENE_TAX_SELECTED"
echo "TAXONOMY_PREFILTER_THRESHOLD   : $TAXONOMY_PREFILTER_THRESHOLD"
echo "SCG_TSV                        : $SCG_TSV"
echo "CLIP_FINAL                     : $CLIP_FINAL"
echo "MAPPING_MANIFEST               : $ASSEMBLER_MAPPING_MANIFEST"
echo "HIFI_MAIN_POLICY               : primary+supplementary only, QNAME-deduplicated; secondary excluded"
echo "HIFI_ANCHOR_BP                 : $CHIMERA_HIFI_ANCHOR_BP"
echo "MIN_SUPPORT_QNAMES             : $CHIMERA_MIN_SUPPORT_QNAMES"
echo "PARALLEL_BAM_JOBS              : $CHIMERA_BAM_JOBS"
echo "ANALYSIS_SCOPE                 : all long-read assembly contigs"
echo "SCG_POLICY                     : MAG-associated contigs only; domain-matched Bacteria_71/Archaea_76 annotation; no direct structural call"
echo "INTER_CONTIG_TAXONOMY          : disabled"

if [[ "$SUMMARY_ONLY" == "Y" ]]; then
    write_mag_outputs
    write_final_taxonomic_outputs
    echo "[07_chimera] Finished: $(date)"
    exit 0
fi

for spec in \
    "ASSEMBLY_FA|$ASSEMBLY_FA" \
    "MAG_CONTIG_MAP|$MAG_CONTIG_MAP" \
    "GTDB_GENE_TAX_TSV|$GENE_TAX_SELECTED" \
    "ASSEMBLER_MAPPING_MANIFEST|$ASSEMBLER_MAPPING_MANIFEST" \
    "CANDIDATE_EVALUATOR|$CANDIDATE_EVALUATOR" \
    "STRUCTURAL_VALIDATOR|$STRUCTURAL_VALIDATOR" \
    "MAG_SUMMARIZER|$MAG_SUMMARIZER" \
    "FINAL_RECLASSIFIER|$FINAL_RECLASSIFIER"; do
    name="${spec%%|*}"; path="${spec#*|}"
    [[ -s "$path" ]] || { echo "ERROR: required input missing or empty: $name=$path" >&2; exit 1; }
done
for tool in samtools python3 sha256sum; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing tool: $tool" >&2; exit 1; }
done

if [[ ! -s "${ASSEMBLY_FA}.fai" || "$ASSEMBLY_FA" -nt "${ASSEMBLY_FA}.fai" ]]; then
    echo "Building/updating assembly FASTA index"
    samtools faidx "$ASSEMBLY_FA"
else
    echo "Reusing existing assembly FASTA index: ${ASSEMBLY_FA}.fai"
fi

# SCG is supplementary only. Generate the shared table when anvi'o is available,
# but never make the structural call depend on SCG availability.
SCG_AVAILABLE=N
if [[ ! -s "$SCG_TSV" ]]; then
    if command -v anvi-gen-contigs-database >/dev/null 2>&1 && command -v anvi-run-hmms >/dev/null 2>&1; then
        echo "SCG table absent; generating shared SCG/HMM table once for chimera annotation and false-closure reuse."
        bash "$SCRIPT_DIR/lib_scg.sh" "$ASSEMBLY_FA" "$SCG_OUT_DIR" "${THREADS:-1}"
    else
        echo "WARNING: SCG table is absent and anvi'o HMM tools are unavailable; structural chimera analysis will continue without SCG annotation." >&2
    fi
fi
if [[ -s "$SCG_TSV" ]]; then
    if [[ -s "$MAG_GTDB_TAXONOMY_TSV" ]]; then
        SCG_AVAILABLE=Y
    else
        echo "WARNING: SCG table is present but MAG GTDB-Tk taxonomy is unavailable; SCG redundancy annotation is disabled." >&2
    fi
fi

# Final taxonomic-chimera classification depends on validated clipping.
[[ -s "$CLIP_FINAL" ]] || {
    echo "ERROR: final taxonomic-chimera calls require current 03_clipping.sh output: $CLIP_FINAL" >&2
    exit 1
}
qc_validate_current_clipping_result || exit 1
CLIPPING_AVAILABLE=Y

# Validate all structural BAMs once before candidate generation.
while IFS=$'\t' read -r assembler assembly hifi_all hifi_primary hifi_secondary hifi_supplementary \
    hifi_primary_supplementary output_dir; do
    [[ "$assembler" == "Assembler" ]] && continue
    for bam in "$hifi_primary_supplementary" "$hifi_primary"; do
        [[ -s "$bam" ]] || { echo "ERROR: missing HiFi BAM for $assembler: $bam" >&2; exit 1; }
        samtools quickcheck -v "$bam" >/dev/null 2>&1 || { echo "ERROR: BAM quickcheck failed: $bam" >&2; exit 1; }
        if [[ ! -s "${bam}.bai" && ! -s "${bam}.csi" ]]; then samtools index -@ "${THREADS:-1}" -c "$bam"; fi
    done
done < "$ASSEMBLER_MAPPING_MANIFEST"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
CANDIDATE_EVENTS="$WORK_DIR/taxonomic_transition_candidates.events.tsv"
CANDIDATE_CONTIGS="$WORK_DIR/taxonomic_transition_candidates.contig.tsv"
DUMMY_SUPPORTED_EVENTS="$WORK_DIR/legacy_supported.events.tsv"
DUMMY_INDEPENDENT_EVENTS="$WORK_DIR/legacy_independent.events.tsv"
DUMMY_SUPPORTED_CONTIG="$WORK_DIR/legacy_supported.contig.tsv"
DUMMY_INDEPENDENT_CONTIG="$WORK_DIR/legacy_independent.contig.tsv"

# Step 1: candidate discovery only. Legacy candidate-class labels remain confined to WORK_DIR.
CANDIDATE_ARGS=(
    --assembly "$ASSEMBLY_FA"
    --fai "${ASSEMBLY_FA}.fai"
    --depth-manifest "$ASSEMBLER_MAPPING_MANIFEST"
    --gene-tax-tsv "$GENE_TAX_SELECTED"
    --mag-map "$MAG_CONTIG_MAP"
    --events-output "$CANDIDATE_EVENTS"
    --supported-events-output "$DUMMY_SUPPORTED_EVENTS"
    --independent-events-output "$DUMMY_INDEPENDENT_EVENTS"
    --contig-output "$CANDIDATE_CONTIGS"
    --supported-contig-output "$DUMMY_SUPPORTED_CONTIG"
    --independent-contig-output "$DUMMY_INDEPENDENT_CONTIG"
    --max-search-genes "$TAX_SEARCH_MAX_GENES"
    --target-classified-genes "$TAX_TARGET_CLASSIFIED_GENES"
    --min-classified-genes "$TAX_MIN_CLASSIFIED_GENES"
    --min-dominance "$TAX_MIN_BLOCK_DOMINANCE"
    --taxonomy-source GTDB
    --taxonomy-release "$GTDB_RELEASE"
    --parallel-bam-jobs "$CHIMERA_BAM_JOBS"
)
if [[ -s "$MAG_GTDB_TAXONOMY_TSV" ]]; then CANDIDATE_ARGS+=(--mag-taxonomy-tsv "$MAG_GTDB_TAXONOMY_TSV"); fi
if [[ "$CLIPPING_AVAILABLE" == "Y" ]]; then CANDIDATE_ARGS+=(--clipping "$CLIP_FINAL"); fi
python3 "$CANDIDATE_EVALUATOR" "${CANDIDATE_ARGS[@]}"

# Step 2: direct structural validation retained as an audit annotation.
VALIDATE_ARGS=(
    --candidate-events "$CANDIDATE_EVENTS"
    --candidate-contigs "$CANDIDATE_CONTIGS"
    --mapping-manifest "$ASSEMBLER_MAPPING_MANIFEST"
    --mag-map "$MAG_CONTIG_MAP"
    --events-output "$EVENTS_OUT"
    --structural-events-output "$STRUCTURAL_EVENTS_OUT"
    --high-confidence-events-output "$HIGH_CONF_EVENTS_OUT"
    --unresolved-events-output "$UNRESOLVED_EVENTS_OUT"
    --contig-output "$CONTIG_OUT"
    --structural-contig-output "$STRUCTURAL_CONTIG_OUT"
    --high-confidence-contig-output "$HIGH_CONF_CONTIG_OUT"
    --unresolved-contig-output "$UNRESOLVED_CONTIG_OUT"
    --scg-audit-output "$SCG_AUDIT_OUT"
    --anchor-bp "$CHIMERA_HIFI_ANCHOR_BP"
    --min-support-qnames "$CHIMERA_MIN_SUPPORT_QNAMES"
    --threads "${THREADS:-1}"
    --parallel-bam-jobs "$CHIMERA_BAM_JOBS"
)
if [[ "$SCG_AVAILABLE" == "Y" ]]; then VALIDATE_ARGS+=(--scg-tsv "$SCG_TSV" --mag-taxonomy-tsv "$MAG_GTDB_TAXONOMY_TSV"); fi
if [[ "$CLIPPING_AVAILABLE" == "Y" ]]; then VALIDATE_ARGS+=(--clipping "$CLIP_FINAL"); fi
python3 "$STRUCTURAL_VALIDATOR" "${VALIDATE_ARGS[@]}"

write_mag_outputs
write_final_taxonomic_outputs

# Candidate work products deliberately removed so no taxonomy-only table can be mistaken
# for a final chimera call.
rm -rf "$WORK_DIR"

{
    printf 'key\tvalue\n'
    printf 'script_id\t%s\n' "$SCRIPT_ID"
    printf 'taxonomy_source\tGTDB\n'
    printf 'taxonomy_release\t%s\n' "$GTDB_RELEASE"
    printf 'gene_taxonomy_tsv\t%s\n' "$GENE_TAX_SELECTED"
    printf 'gene_taxonomy_table_sha256\t%s\n' "$(sha256sum "$GENE_TAX_SELECTED" | awk '{print $1}')"
    printf 'taxonomy_prefilter_threshold\t%s\n' "$TAXONOMY_PREFILTER_THRESHOLD"
    printf 'analysis_scope\tall_long_read_assembly_contigs\n'
    printf 'candidate_policy\tbounded_local_gene_taxonomy_transition_only\n'
    printf 'final_call_policy\ttaxonomic_shift_plus_validated_Clipping_Y\n'
    printf 'hifi_main_alignment_policy\tprimary_plus_supplementary_only_secondary_excluded_QNAME_deduplicated\n'
    printf 'hifi_anchor_bp\t%s\n' "$CHIMERA_HIFI_ANCHOR_BP"
    printf 'minimum_support_qnames\t%s\n' "$CHIMERA_MIN_SUPPORT_QNAMES"
    printf 'parallel_bam_jobs\t%s\n' "$CHIMERA_BAM_JOBS"
    printf 'validated_clipping_policy\treuse_03_clipping_no_recomputation\n'
    printf 'high_confidence_policy\ttaxonomic_shift_plus_Stringent_clipping_Y\n'
    printf 'primary_supp_structural_validation_policy\taudit_only_not_final_classifier\n'
    printf 'scg_policy\tMAG_associated_contigs_only_domain_matched_Bacteria_71_or_Archaea_76_redundancy_annotation_no_structural_call\n'
    printf 'scg_available\t%s\n' "$SCG_AVAILABLE"
    printf 'gc_depth_policy\tdescriptive_only\n'
    printf 'inter_contig_taxonomy_policy\tdisabled\n'
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$RUN_META_OUT"

echo "Outputs:"
for f in "$EVENTS_OUT" "$STRUCTURAL_EVENTS_OUT" "$HIGH_CONF_EVENTS_OUT" "$UNRESOLVED_EVENTS_OUT" \
    "$CONTIG_OUT" "$STRUCTURAL_CONTIG_OUT" "$HIGH_CONF_CONTIG_OUT" "$UNRESOLVED_CONTIG_OUT" \
    "$SCG_AUDIT_OUT" "$MAG_CONTIG_OUT" "$MAG_OUT" "$RUN_META_OUT"; do echo "  $f"; done
echo "  Final clipping-defined taxonomic-chimera outputs: $FINAL_CHIMERA_DIR"
echo "[07_chimera] Finished: $(date)"
