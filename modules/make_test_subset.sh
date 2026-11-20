#!/usr/bin/env bash
# Build a small MAGverify test dataset from a full one.
#
# Keeps a few MAGs (whole FASTAs), the matching contigs from each assembler FASTA
# (same IDs and sequences), some unbinned contigs, and only the HiFi reads that map
# to this reduced reference. FASTAs are streamed, so memory use stays small.
#
# Usage:
#   bash make_test_subset.sh \
#     --contigs-dir DIR --mag-dir DIR --reads FASTQ[.gz] --out DIR \
#     [--mags Bin1L,Bin2L | --n-mags 6] [--unbinned-per-assembler 200] [--threads 32]
set -euo pipefail

CONTIGS_DIR="" MAG_DIR="" READS="" OUT=""
MAGS="" N_MAGS=6 UNBINNED=200 THREADS=32
while [[ $# -gt 0 ]]; do
    case "$1" in
        --contigs-dir) CONTIGS_DIR="$2"; shift 2 ;;
        --mag-dir) MAG_DIR="$2"; shift 2 ;;
        --reads) READS="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --mags) MAGS="$2"; shift 2 ;;
        --n-mags) N_MAGS="$2"; shift 2 ;;
        --unbinned-per-assembler) UNBINNED="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) sed -n 2,11p "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
for v in CONTIGS_DIR MAG_DIR READS OUT; do
    [[ -n "${!v}" ]] || { echo "Missing --${v,,}" | tr _ - >&2; exit 2; }
done
for exe in python3 minimap2 samtools gzip; do
    command -v "$exe" >/dev/null || { echo "Missing command: $exe" >&2; exit 1; }
done

mkdir -p "$OUT/contigs" "$OUT/MAGs"

# 1) Choose MAGs and write the reduced assembler FASTAs and MAG FASTAs.
python3 - "$CONTIGS_DIR" "$MAG_DIR" "$OUT" "$MAGS" "$N_MAGS" "$UNBINNED" <<'PY'
import random
import shutil
import sys
from pathlib import Path

contigs_dir, mag_dir, out, mags_arg, n_mags, unbinned = sys.argv[1:]
contigs_dir, mag_dir, out = Path(contigs_dir), Path(mag_dir), Path(out)
n_mags, unbinned = int(n_mags), int(unbinned)
EXT = ('.fa', '.fna', '.fasta')
MAX_UNBINNED_BP = 200_000


def read_fasta(path):
    name, seq = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def write_record(handle, name, seq):
    handle.write(f'>{name}\n')
    for i in range(0, len(seq), 80):
        handle.write(seq[i:i + 80] + '\n')


mag_info = {}
for path in sorted(mag_dir.iterdir()):
    if path.suffix.lower() in EXT:
        names, size = [], 0
        for name, seq in read_fasta(path):
            names.append(name)
            size += len(seq)
        mag_info[path.stem] = (path, names, size)

if mags_arg:
    chosen = [m.strip() for m in mags_arg.split(',') if m.strip()]
    missing = [m for m in chosen if m not in mag_info]
    if missing:
        sys.exit(f'MAGs not found in --mag-dir: {missing}')
else:
    # Smallest MAGs keep the run short; include single-contig MAGs for circularity.
    by_size = sorted(mag_info, key=lambda m: mag_info[m][2])
    single = [m for m in by_size if len(mag_info[m][1]) == 1]
    multi = [m for m in by_size if len(mag_info[m][1]) > 1]
    half = max(1, n_mags // 2)
    chosen = (single[:half] + multi)[:n_mags]
    if len(chosen) < n_mags:
        chosen += [m for m in single[half:] if m not in chosen][:n_mags - len(chosen)]

keep = {c for m in chosen for c in mag_info[m][1]}
found, total = set(), 0
rng = random.Random(1)
for asm in sorted(p for p in contigs_dir.iterdir() if p.suffix.lower() in EXT):
    binned, reservoir, seen = 0, [], 0
    with open(out / 'contigs' / asm.name, 'w') as handle:
        for name, seq in read_fasta(asm):
            if name in keep:
                write_record(handle, name, seq)
                found.add(name)
                binned += 1
            elif len(seq) <= MAX_UNBINNED_BP and unbinned:
                seen += 1
                if len(reservoir) < unbinned:
                    reservoir.append((name, seq))
                else:
                    j = rng.randrange(seen)
                    if j < unbinned:
                        reservoir[j] = (name, seq)
        for name, seq in reservoir:
            write_record(handle, name, seq)
    total += binned + len(reservoir)
    print(f'{asm.name}: {binned} MAG contigs + {len(reservoir)} unbinned contigs', flush=True)
absent = keep - found
if absent:
    sys.exit(f'MAG contigs not found in any assembler FASTA: {sorted(absent)[:5]}')

for m in chosen:
    path = mag_info[m][0]
    shutil.copyfile(path, out / 'MAGs' / path.name)
    print(f'MAG {m}: {len(mag_info[m][1])} contigs, {mag_info[m][2]:,} bp')
print(f'Total contigs in reduced reference: {total}')
PY

# 2) Keep only the HiFi reads whose primary alignment falls on the reduced reference.
cat "$OUT"/contigs/*.fa* > "$OUT/.subset_reference.fa"
echo "Selecting HiFi reads that map to the reduced reference..."
minimap2 -t "$THREADS" -ax map-hifi "$OUT/.subset_reference.fa" "$READS" \
    | samtools fastq -@ 4 -F 0x904 - \
    | gzip -c > "$OUT/read.fastq.gz"
rm -f "$OUT/.subset_reference.fa"

n_reads=$(( $(gzip -cd "$OUT/read.fastq.gz" | wc -l) / 4 ))
echo
echo "Test dataset: $OUT"
echo "  contigs/       $(ls "$OUT/contigs" | wc -l) assembler FASTAs ($(du -sh "$OUT/contigs" | cut -f1))"
echo "  MAGs/          $(ls "$OUT/MAGs" | wc -l) MAG FASTAs"
echo "  read.fastq.gz  $n_reads reads ($(du -sh "$OUT/read.fastq.gz" | cut -f1))"
