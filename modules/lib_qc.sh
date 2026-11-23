#!/usr/bin/env bash

qc_sq_count() {
    local bam="$1"
    samtools view -H "$bam" 2>/dev/null \
        | awk 'BEGIN{n=0} /^@SQ/{n++} END{print n+0}'
}

qc_bam_is_valid() {
    local bam="$1"
    local n_sq
    [[ -s "$bam" ]] || return 1
    samtools quickcheck "$bam" >/dev/null 2>&1 || return 1
    n_sq="$(qc_sq_count "$bam")"
    [[ "$n_sq" =~ ^[0-9]+$ ]] && (( n_sq > 0 ))
}

qc_bam_matches_fai() {
    local bam="$1"
    local fai="$2"

    qc_bam_is_valid "$bam" || return 1
    python3 - "$bam" "$fai" <<'PY_QC_REFS'
import subprocess
import sys

bam, fai = sys.argv[1:]
expected = {}
with open(fai) as handle:
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        if len(fields) >= 2:
            name, length = fields[0], int(fields[1])
            if name in expected:
                raise SystemExit(f'ERROR: duplicate FASTA sequence ID: {name}')
            expected[name] = length

observed = {}
header = subprocess.check_output(['samtools', 'view', '-H', bam], text=True)
for line in header.splitlines():
    if not line.startswith('@SQ\t'):
        continue
    tags = {}
    for field in line.split('\t')[1:]:
        if ':' in field:
            key, value = field.split(':', 1)
            tags[key] = value
    if 'SN' in tags and 'LN' in tags:
        name = tags['SN']
        if name in observed:
            raise SystemExit(f'ERROR: duplicate BAM @SQ SN: {name}')
        observed[name] = int(tags['LN'])

if observed != expected:
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    wrong = sorted(
        name for name in set(expected) & set(observed)
        if expected[name] != observed[name]
    )
    print(
        f'BAM/reference mismatch: missing={len(missing)}, '
        f'unexpected={len(unexpected)}, wrong_length={len(wrong)}',
        file=sys.stderr,
    )
    raise SystemExit(1)
PY_QC_REFS
}

qc_compute_mapping_fingerprint() {
    local assembly="$1"
    local reads="$2"
    local preset="$3"
    local index_chunk="$4"
    local max_secondary="${MINIMAP2_MAX_SECONDARY:-5}"

    python3 - "$assembly" "$reads" "$preset" "$index_chunk" "$max_secondary" <<'PY_QC_FP'
import hashlib
import json
import os
import subprocess
import sys

assembly, reads, preset, index_chunk, max_secondary = sys.argv[1:]

def stat_payload(path):
    st = os.stat(path)
    return {
        'realpath': os.path.realpath(path),
        'size': st.st_size,
        'mtime_ns': st.st_mtime_ns,
    }

fai = assembly + '.fai'
digest = hashlib.sha256()
with open(fai, 'rb') as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        digest.update(chunk)

def command_first_line(command):
    raw = subprocess.check_output(command, stderr=subprocess.STDOUT)
    lines = raw.splitlines()
    if not lines:
        raise SystemExit(f'ERROR: command returned no version output: {command[0]}')
    return lines[0].decode('utf-8', errors='replace').strip()

version = command_first_line(['minimap2', '--version'])
samtools_version = command_first_line(['samtools', '--version'])

payload = {
    'assembly': stat_payload(assembly),
    'assembly_fai_sha256': digest.hexdigest(),
    'reads': stat_payload(reads),
    'preset': preset,
    'index_chunk': index_chunk,
    'mapping_arguments': [
        '-ax', preset, '-p', '1', '--secondary-seq', '-I', index_chunk,
        'samtools view -F 4'
    ],
    'minimap2_defaults_assumed': {
        'secondary': 'yes',
        'max_secondary': int(max_secondary),
    },
    'minimap2_version': version,
    'samtools_version': samtools_version,
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_QC_FP
}

qc_hash_text() {
    python3 - "$@" <<'PY_QC_TEXT_HASH'
import hashlib
import sys
print(hashlib.sha256('\0'.join(sys.argv[1:]).encode()).hexdigest())
PY_QC_TEXT_HASH
}

qc_compute_mapping_stage_fingerprints() {
    local assembly="$1"
    local reads="$2"
    local preset="$3"
    local index_chunk="$4"
    local mapping_fp samtools_version sorted_fp primary_fp

    mapping_fp="$(qc_compute_mapping_fingerprint \
        "$assembly" "$reads" "$preset" "$index_chunk")" || return 1
    samtools_version="$(samtools --version | awk 'NR==1{print}')" || return 1
    sorted_fp="$(qc_hash_text \
        "$mapping_fp" "sort_mapped_only" "$samtools_version")" || return 1
    primary_fp="$(qc_hash_text \
        "$sorted_fp" "samtools view -F 2308")" || return 1
    printf '%s\t%s\t%s\n' "$mapping_fp" "$sorted_fp" "$primary_fp"
}

qc_cache_matches() {
    local path="$1"
    local expected="$2"
    [[ -s "$path" ]] && [[ "$(<"$path")" == "$expected" ]]
}

qc_hifi_bam_mapping_semantically_matches() {
    local bam="$1"
    local assembly="$2"
    local reads="$3"
    local preset="$4"
    local index_chunk="$5"
    local max_secondary="${6:-${MINIMAP2_MAX_SECONDARY:-5}}"

    [[ -s "$bam" && -s "$assembly" && -s "$reads" && -s "${assembly}.fai" ]] || return 1
    samtools quickcheck "$bam" >/dev/null 2>&1 || return 1
    qc_bam_matches_fai "$bam" "${assembly}.fai" || return 1

    python3 - "$bam" "$assembly" "$reads" "$preset" "$index_chunk" "$max_secondary" <<'PY_QC_HIFI_SEMANTIC'
import os
import re
import shlex
import subprocess
import sys

bam, assembly, reads, preset, expected_chunk, expected_n = sys.argv[1:]
expected_n = int(expected_n)

try:
    header = subprocess.check_output(
        ["samtools", "view", "-H", bam], stderr=subprocess.STDOUT, text=True
    )
except Exception as exc:
    print(f"HiFi BAM semantic audit failed: cannot read header: {exc}", file=sys.stderr)
    raise SystemExit(1)

hd = next((line for line in header.splitlines() if line.startswith("@HD\t")), "")
if "\tSO:coordinate" not in hd:
    print("HiFi BAM semantic audit failed: BAM is not coordinate-sorted", file=sys.stderr)
    raise SystemExit(1)

pg_lines = [line for line in header.splitlines() if line.startswith("@PG\t")]
commands = []
for line in pg_lines:
    fields = line.split("\t")[1:]
    rec = {}
    for field in fields:
        if ":" in field:
            key, value = field.split(":", 1)
            rec[key] = value
    if rec.get("PN") == "minimap2" or "minimap2" in rec.get("CL", ""):
        commands.append(rec.get("CL", ""))

if not commands:
    print("HiFi BAM semantic audit failed: minimap2 @PG command not found", file=sys.stderr)
    raise SystemExit(1)

def option_value(tokens, opt):
    for i, tok in enumerate(tokens):
        if tok == opt and i + 1 < len(tokens):
            return tokens[i + 1]
        if tok.startswith(opt + "="):
            return tok.split("=", 1)[1]
        if len(opt) == 2 and tok.startswith(opt) and tok != opt:
            return tok[len(opt):]
    return None

def parse_size(value):
    m = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([kKmMgGtT]?)", value or "")
    if not m:
        return None
    number = float(m.group(1))
    suffix = m.group(2).upper()
    scale = {"": 1, "K": 1_000, "M": 1_000_000, "G": 1_000_000_000, "T": 1_000_000_000_000}[suffix]
    return int(number * scale)

reference_bases = 0
with open(assembly + ".fai", "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            reference_bases += int(parts[1])

expected_chunk_size = parse_size(expected_chunk)
assembly_real = os.path.realpath(assembly)
reads_real = os.path.realpath(reads)

failures = []
for command in commands:
    try:
        tokens = shlex.split(command)
    except ValueError:
        continue

    try:
        mm_i = next(i for i, tok in enumerate(tokens) if os.path.basename(tok) == "minimap2")
        tokens = tokens[mm_i + 1:]
    except StopIteration:
        pass

    preset_value = option_value(tokens, "-ax")
    if preset_value is None:
        x_value = option_value(tokens, "-x")
        has_a = "-a" in tokens
        preset_value = x_value if has_a else None
    p_value = option_value(tokens, "-p")
    n_value = option_value(tokens, "-N")
    chunk_value = option_value(tokens, "-I")

    if preset_value != preset:
        failures.append(f"preset={preset_value!r}")
        continue
    try:
        if float(p_value) != 1.0:
            failures.append(f"-p={p_value!r}")
            continue
    except (TypeError, ValueError):
        failures.append(f"-p={p_value!r}")
        continue
    try:
        if int(n_value) != expected_n:
            failures.append(f"-N={n_value!r}")
            continue
    except (TypeError, ValueError):
        failures.append(f"-N={n_value!r}")
        continue
    if "--secondary-seq" not in tokens:
        failures.append("--secondary-seq absent")
        continue
    if "--secondary=no" in tokens or "--secondary=no" in command:
        failures.append("secondary alignments disabled")
        continue

    observed_chunk_size = parse_size(chunk_value)
    chunk_ok = chunk_value == expected_chunk
    if not chunk_ok and observed_chunk_size is not None and expected_chunk_size is not None:
        chunk_ok = observed_chunk_size >= reference_bases and expected_chunk_size >= reference_bases
    if not chunk_ok:
        failures.append(f"-I={chunk_value!r}")
        continue

    existing_paths = []
    for tok in tokens:
        if tok == "-":
            continue
        try:
            if os.path.exists(tok):
                existing_paths.append(os.path.realpath(tok))
        except OSError:
            pass
    if assembly_real not in existing_paths:
        failures.append("assembly path absent from minimap2 command")
        continue
    if reads_real not in existing_paths:
        failures.append("read path absent from minimap2 command")
        continue

    raise SystemExit(0)

reason = failures[-1] if failures else "no parseable compatible minimap2 command"
print(f"HiFi BAM semantic audit failed: {reason}", file=sys.stderr)
raise SystemExit(1)
PY_QC_HIFI_SEMANTIC
}

qc_write_cache_atomic() {
    local path="$1"
    local value="$2"
    local tmp="${path}.tmp.$$"
    mkdir -p "$(dirname "$path")"
    printf '%s\n' "$value" > "$tmp"
    mv -f "$tmp" "$path"
}

qc_file_sha256() {
    local path="$1"
    python3 - "$path" <<'PY_QC_FILE_HASH'
import hashlib
import sys

digest = hashlib.sha256()
with open(sys.argv[1], 'rb') as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        digest.update(chunk)
print(digest.hexdigest())
PY_QC_FILE_HASH
}

qc_result_cache_matches() {
    local meta="$1"
    local input_fingerprint="$2"
    local result="$3"
    local cached_fingerprint cached_hash extra actual_hash

    [[ -s "$meta" && -s "$result" ]] || return 1
    IFS=$'\t' read -r cached_fingerprint cached_hash extra < "$meta"
    [[ -z "${extra:-}" ]] || return 1
    [[ "$cached_fingerprint" == "$input_fingerprint" ]] || return 1
    [[ "$cached_hash" =~ ^[0-9a-fA-F]{64}$ ]] || return 1
    actual_hash="$(qc_file_sha256 "$result")" || return 1
    [[ "${actual_hash,,}" == "${cached_hash,,}" ]]
}

qc_write_result_cache_atomic() {
    local meta="$1"
    local input_fingerprint="$2"
    local result="$3"
    local result_hash tmp

    result_hash="$(qc_file_sha256 "$result")"
    tmp="${meta}.tmp.$$"
    mkdir -p "$(dirname "$meta")"
    printf '%s\t%s\n' "$input_fingerprint" "$result_hash" > "$tmp"
    mv -f "$tmp" "$meta"
}

qc_compute_anvio_clipping_raw_fingerprint() {
    local clip_bam_meta="$1"
    local default_ratio="$2"
    local min_dist_to_end="$3"
    local algorithm_id="$4"
    local executable="$5"
    local anvio_version="$6"

    [[ -s "$clip_bam_meta" ]] || return 1
    [[ -x "$executable" ]] || return 1
    qc_hash_text \
        "$(<"$clip_bam_meta")" \
        "$algorithm_id" \
        "$default_ratio" \
        "$min_dist_to_end" \
        "$(readlink -f "$executable")" \
        "$anvio_version" \
        "--just-do-it"
}

qc_compute_anvio_clipping_result_fingerprint() {
    local raw_cache_meta="$1"
    local raw_clipping="$2"
    local stringent_depth="$3"
    local stringent_ratio="$4"
    local normalization_algorithm="$5"
    local validation_algorithm="$6"
    local normalizer_script="$7"
    local validator_script="$8"

    [[ -s "$raw_cache_meta" ]] || return 1
    [[ -s "$raw_clipping" ]] || return 1
    [[ -s "$normalizer_script" ]] || return 1
    [[ -s "$validator_script" ]] || return 1
    qc_hash_text \
        "$(<"$raw_cache_meta")" \
        "$(qc_file_sha256 "$raw_clipping")" \
        "$normalization_algorithm" \
        "$stringent_depth" \
        "$stringent_ratio" \
        "$(qc_file_sha256 "$normalizer_script")" \
        "$validation_algorithm" \
        "$(qc_file_sha256 "$validator_script")"
}

qc_compute_source_specific_clipping_fingerprint() {
    local manifest="$1"
    local default_ratio="$2"
    local min_dist_to_end="$3"
    local terminal_window="$4"
    local stringent_depth="$5"
    local stringent_ratio="$6"
    local discovery_algorithm="$7"
    local normalization_algorithm="$8"
    local validation_algorithm="$9"
    shift 9
    local normalizer_script="$1"
    local validator_script="$2"
    local classifier_script="$3"
    local executable="$4"
    local anvio_version="$5"

    python3 - \
        "$manifest" "$default_ratio" "$min_dist_to_end" "$terminal_window" \
        "$stringent_depth" "$stringent_ratio" \
        "$discovery_algorithm" "$normalization_algorithm" "$validation_algorithm" \
        "$normalizer_script" "$validator_script" "$classifier_script" \
        "$executable" "$anvio_version" <<'PY_QC_SOURCE_CLIP_FP'
import csv
import hashlib
import json
import os
import sys

(
    manifest, default_ratio, min_dist, terminal_window, stringent_depth,
    stringent_ratio, discovery_algorithm, normalization_algorithm,
    validation_algorithm, normalizer, validator, classifier, executable,
    anvio_version,
) = sys.argv[1:]

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def manifest_sha256(path):
    # The first seven manifest columns only: the flagstat counts appended by map do not change results.
    with open(path, 'rb') as handle:
        lines = handle.read().split(b'\n')
    core = b'\n'.join(b'\t'.join(line.split(b'\t')[:7]) for line in lines)
    return hashlib.sha256(core).hexdigest()

def stat_payload(path):
    st = os.stat(path)
    return {
        'realpath': os.path.realpath(path),
        'size': st.st_size,
        'mtime_ns': st.st_mtime_ns,
    }

for path in (manifest, normalizer, validator, classifier, executable):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f'ERROR: clipping fingerprint input missing: {path}')

sources = []
with open(manifest, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'Assembler', 'Assembly_FASTA', 'HiFi_all_BAM'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'ERROR: mapping manifest lacks columns: {sorted(missing)}')
    for row in reader:
        assembler = (row.get('Assembler') or '').strip()
        assembly = (row.get('Assembly_FASTA') or '').strip()
        bam = (row.get('HiFi_all_BAM') or '').strip()
        if not assembler or not assembly or not bam:
            raise SystemExit('ERROR: blank mapping manifest field')
        fai = assembly + '.fai'
        bam_meta = os.path.join(os.path.dirname(bam), '.cache', os.path.basename(bam) + '.inputs.sha256')
        for path in (assembly, fai, bam, bam_meta):
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise SystemExit(f'ERROR: clipping input missing: {path}')
        sources.append({
            'assembler': assembler,
            'assembly': stat_payload(assembly),
            'assembly_fai_sha256': sha256(fai),
            'hifi_bam': stat_payload(bam),
            'hifi_bam_mapping_fingerprint': open(bam_meta, encoding='utf-8').read().strip(),
        })

sources.sort(key=lambda x: x['assembler'])
payload = {
    'manifest_sha256': manifest_sha256(manifest),
    'sources': sources,
    'candidate_discovery': {
        'clipping_ratio': float(default_ratio),
        'min_dist_to_end': int(min_dist),
        'terminal_window_bp': int(terminal_window),
        'algorithm': discovery_algorithm,
        'anvio_executable': os.path.realpath(executable),
        'anvio_version': anvio_version,
    },
    'high_confidence': {
        'min_depth': int(stringent_depth),
        'clipping_ratio': float(stringent_ratio),
    },
    'normalization_algorithm': normalization_algorithm,
    'validation_algorithm': validation_algorithm,
    'normalizer_sha256': sha256(normalizer),
    'validator_sha256': sha256(validator),
    'classifier_sha256': sha256(classifier),
}
print(hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
PY_QC_SOURCE_CLIP_FP
}

qc_validate_current_clipping_result() {
    local anvio_version expected_final_fp

    [[ -s "$ASSEMBLER_MAPPING_MANIFEST" ]] || {
        echo "ERROR: assembler mapping manifest missing; run 01_mapping.sh" >&2
        return 1
    }
    [[ -s "$CLIP_FINAL" && -s "$CLIP_FINAL_CACHE_META" ]] || {
        echo "ERROR: clipping result/cache missing; run 03_clipping.sh" >&2
        return 1
    }
    [[ -x "$ANVIO_FIND_MISASSEMBLIES" ]] || return 1
    for path in "$CLIP_VALIDATE_SCRIPT" "$CLIP_SUMMARIZE_SCRIPT"; do
        [[ -s "$path" ]] || return 1
    done

    while IFS=$'\t' read -r assembler assembly hifi_bam _rest; do
        [[ "$assembler" == "Assembler" ]] && continue
        [[ -n "$assembler" && -s "$assembly" && -s "${assembly}.fai" && -s "$hifi_bam" ]] || {
            echo "ERROR: incomplete clipping input for $assembler" >&2
            return 1
        }
        samtools quickcheck -v "$hifi_bam" >/dev/null 2>&1 || {
            echo "ERROR: HiFi BAM failed quickcheck: $hifi_bam" >&2
            return 1
        }
        qc_bam_matches_fai "$hifi_bam" "${assembly}.fai" || {
            echo "ERROR: HiFi BAM/reference mismatch: $assembler" >&2
            return 1
        }
    done < "$ASSEMBLER_MAPPING_MANIFEST"

    anvio_version="$($ANVIO_ENV_PREFIX/bin/python - <<'PY_QC_ANVIO_VERSION'
import anvio
print(anvio.__version__)
PY_QC_ANVIO_VERSION
)" || return 1

    expected_final_fp="$(qc_compute_source_specific_clipping_fingerprint \
        "$ASSEMBLER_MAPPING_MANIFEST" "$ANVIO_DEFAULT_CLIPPING_RATIO" \
        "$MIN_DIST_TO_END" "$TERMINAL_CLIPPING_WINDOW_BP" \
        "$STRINGENT_CLIP_MIN_DEPTH" "$STRINGENT_CLIP_RATIO" \
        "$ANVIO_CLIPPING_ALGORITHM_ID" "$CLIPPING_NORMALIZATION_ALGORITHM_ID" \
        "$CLIPPING_VALIDATION_ALGORITHM_ID" "$CLIP_VALIDATE_SCRIPT" \
        "$CLIP_VALIDATE_SCRIPT" "$CLIP_SUMMARIZE_SCRIPT" \
        "$ANVIO_FIND_MISASSEMBLIES" "$anvio_version")" || return 1

    qc_result_cache_matches \
        "$CLIP_FINAL_CACHE_META" "$expected_final_fp" "$CLIP_FINAL" || {
        echo "ERROR: final clipping cache is stale; run 03_clipping.sh" >&2
        return 1
    }

    python3 - "$CLIP_FINAL" <<'PY_QC_INTERNAL_CLIP'
import csv
import sys
path = sys.argv[1]
with open(path, encoding='utf-8-sig', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'Assembler', 'region_type', 'Clipping'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'ERROR: clipping table lacks new schema columns: {sorted(missing)}')
    bad = [row for row in reader if (row.get('region_type') or '').strip() != 'INTERNAL']
if bad:
    raise SystemExit(f'ERROR: internal clipping table must contain INTERNAL rows only: {path}')
PY_QC_INTERNAL_CLIP
}

qc_parallel_jobs() {
    local items="$1" jobs="${PARALLEL_JOBS:-$1}"
    (( jobs > items )) && jobs="$items"
    (( jobs < 1 )) && jobs=1
    echo "$jobs"
}

qc_duration() {
    local s="$1"
    if (( s >= 3600 )); then printf '%dh%02dm' $(( s / 3600 )) $(( s % 3600 / 60 ))
    elif (( s >= 60 )); then printf '%dm' $(( s / 60 ))
    else printf '%ds' "$s"; fi
}

qc_progress_watch() {
    local label="$1" total="$2" dir="$3" pattern="$4" parent="$$"
    local interval="${PROGRESS_INTERVAL:-300}" start="$SECONDS" finished elapsed remaining
    (( total > 0 )) || return 0
    while sleep "$interval"; do
        kill -0 "$parent" 2>/dev/null || break
        finished="$(find "$dir" -maxdepth 1 -name "$pattern" 2>/dev/null | wc -l)"
        elapsed=$(( SECONDS - start ))
        remaining="estimating"
        (( finished > 0 )) && remaining="$(qc_duration $(( elapsed * (total - finished) / finished )))"
        echo "[$label] $finished/$total ($(( 100 * finished / total ))%) | elapsed $(qc_duration "$elapsed") | remaining $remaining"
    done
}

qc_progress_stop() {
    kill "$1" 2>/dev/null || true
    wait "$1" 2>/dev/null || true
}

qc_split_threads() {
    local threads="$1" jobs="$2" each
    each=$(( threads / jobs ))
    (( each < 1 )) && each=1
    echo "$each"
}

qc_wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do
        sleep 1
    done
}

qc_wait_all() {
    local pid
    for pid in $(jobs -rp); do
        wait "$pid" 2>/dev/null || true
    done
}

qc_manifest_columns() {
    python3 - "$@" <<'PY_MANIFEST_COLUMNS'
import csv, sys
path, *columns = sys.argv[1:]
with open(path, encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    missing = [c for c in columns if c not in (reader.fieldnames or [])]
    if missing:
        sys.exit(f"ERROR: mapping manifest lacks columns {missing}: {path}")
    for row in reader:
        if (row.get("Assembler") or "").strip():
            print("\t".join((row.get(c) or "").strip() for c in columns))
PY_MANIFEST_COLUMNS
}

qc_bam_cache() {
    printf '%s/.cache/%s.inputs.sha256\n' "$(dirname "$1")" "$(basename "$1")"
}
