#!/usr/bin/env python3
"""Build the MAG-to-contig map (MAG_ID, Contig) from MAG FASTAs and the all-assemblies FAI."""
import argparse
import csv
import os
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description='Build/validate MAG-to-contig ownership map from MAG FASTAs and assembly FAI.')
    p.add_argument('--mag-dir', required=True)
    p.add_argument('--extensions', default='fa fasta fna')
    p.add_argument('--assembly-fai', required=True)
    p.add_argument('--output', required=True)
    return p.parse_args()


def canonical(path: Path) -> str:
    name = path.name
    for suffix in ('.fasta', '.fna', '.fa'):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return name


def fasta_records(path: Path):
    name = None
    length = 0
    seen = set()
    with path.open(encoding='utf-8') as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, length
                name = line[1:].split()[0]
                if not name or name in seen:
                    raise SystemExit(f'ERROR: empty/duplicate FASTA ID in {path}: {name}')
                seen.add(name)
                length = 0
            else:
                if name is None:
                    raise SystemExit(f'ERROR: sequence before FASTA header in {path}')
                length += len(line)
    if name is not None:
        yield name, length


def main():
    args = parse_args()
    mag_dir = Path(args.mag_dir)
    fai = Path(args.assembly_fai)
    output = Path(args.output)

    if not mag_dir.is_dir():
        raise SystemExit(f'ERROR: MAG directory not found: {mag_dir}')
    if not fai.is_file() or fai.stat().st_size == 0:
        raise SystemExit(f'ERROR: assembly FAI missing or empty: {fai}')

    assembly = {}
    with fai.open(encoding='utf-8') as handle:
        for raw in handle:
            fields = raw.rstrip('\n').split('\t')
            if len(fields) < 2:
                continue
            name = fields[0]
            if name in assembly:
                raise SystemExit(f'ERROR: duplicate contig in assembly FAI: {name}')
            assembly[name] = int(fields[1])
    if not assembly:
        raise SystemExit(f'ERROR: assembly FAI has no contigs: {fai}')

    paths = set()
    for raw_ext in args.extensions.split():
        ext = raw_ext.strip().lstrip('.')
        if ext:
            paths.update(mag_dir.glob(f'*.{ext}'))
    if not paths:
        raise SystemExit(f'ERROR: no MAG FASTA files found in {mag_dir}')

    target_rows = []
    target_owner = {}
    target_by_mag = {}
    for path in sorted(paths):
        mag = canonical(path)
        target_by_mag.setdefault(mag, set())
        found = False
        for contig, length in fasta_records(path):
            found = True
            if contig not in assembly or assembly[contig] != length:
                raise SystemExit(
                    f'ERROR: MAG/assembly contig name or length mismatch: {mag}/{contig}'
                )
            if contig in target_owner and target_owner[contig] != mag:
                raise SystemExit(f'ERROR: contig occurs in multiple MAG files: {contig}')
            target_owner[contig] = mag
            target_by_mag[mag].add(contig)
            target_rows.append((mag, contig))
        if not found:
            raise SystemExit(f'ERROR: empty MAG FASTA: {path}')

    output.parent.mkdir(parents=True, exist_ok=True)
    if not output.is_file() or output.stat().st_size == 0:
        tmp = Path(f'{output}.tmp.{os.getpid()}')
        with tmp.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['MAG_ID', 'Contig'])
            writer.writerows(sorted(target_rows))
        os.replace(tmp, output)
        print(f'[PREP] Generated MAG_CONTIG_MAP: {output}')

    mapped_owner = {}
    mapped_by_mag = {}
    seen_rows = set()
    with output.open(encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = {name.lower(): name for name in (reader.fieldnames or [])}
        mag_col = fields.get('mag_id') or fields.get('id')
        contig_col = fields.get('contig') or fields.get('original_contig')
        if not mag_col or not contig_col:
            raise SystemExit('ERROR: MAG_CONTIG_MAP needs MAG_ID and Contig columns')
        for row in reader:
            mag = row[mag_col].strip()
            contig = row[contig_col].strip()
            if not mag or not contig or (mag, contig) in seen_rows:
                raise SystemExit(f'ERROR: blank or duplicate MAG map row: {mag}/{contig}')
            seen_rows.add((mag, contig))
            if contig not in assembly:
                raise SystemExit(f'ERROR: MAG map contig absent from assembly: {contig}')
            if contig in mapped_owner and mapped_owner[contig] != mag:
                raise SystemExit(f'ERROR: contig assigned to multiple MAGs in map: {contig}')
            mapped_owner[contig] = mag
            mapped_by_mag.setdefault(mag, set()).add(contig)

    for contig, expected_mag in target_owner.items():
        if mapped_owner.get(contig) != expected_mag:
            raise SystemExit(
                f'ERROR: target MAG map mismatch for {contig}: '
                f'{mapped_owner.get(contig, "missing")} versus {expected_mag}'
            )
    for mag, expected_contigs in target_by_mag.items():
        observed_contigs = mapped_by_mag.get(mag, set())
        if observed_contigs != expected_contigs:
            missing = sorted(expected_contigs - observed_contigs)
            extra = sorted(observed_contigs - expected_contigs)
            raise SystemExit(
                f'ERROR: target MAG FASTA/map membership differs for {mag}: '
                f'missing={missing[:10]}, extra={extra[:10]}'
            )

    print(f'[PREP] Validated MAG_CONTIG_MAP rows: {len(seen_rows)}')


if __name__ == '__main__':
    main()
