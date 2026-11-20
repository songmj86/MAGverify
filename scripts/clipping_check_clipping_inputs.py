#!/usr/bin/env python3
"""Check the required data files before running source-specific clipping."""

import argparse
import csv
from pathlib import Path


FIELDS = (
    "Assembler",
    "Assembly_FASTA",
    "HiFi_all_BAM",
    "HiFi_primary_BAM",
    "HiFi_primary_supplementary_BAM",
    "Output_directory",
)


def check_file(path, label):
    if not path or not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"{label} missing or empty: {path or '<unset>'}")


def check_inputs(assembly, mag_dir, manifest):
    check_file(assembly, "All-assemblies FASTA (In)")
    if assembly.suffix.lower() == ".gz":
        raise ValueError("All-assemblies FASTA must be uncompressed for samtools faidx")
    if not mag_dir or not mag_dir.is_dir():
        raise ValueError(
            f"MAG_DIR missing: {mag_dir or '<unset>'}. "
            "Use MAGverify clipping -g /path/to/MAGs"
        )
    mags = [p for p in mag_dir.iterdir() if p.is_file() and p.suffix.lower() in (".fa", ".fna", ".fasta")]
    if not mags:
        raise ValueError(f"No MAG FASTA files (.fa/.fna/.fasta) in MAG_DIR: {mag_dir}")
    check_file(manifest, "ASSEMBLER_MAPPING_MANIFEST (run 01_mapping.sh first)")
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(FIELDS) - set(reader.fieldnames or []):
            raise ValueError(
                "Invalid mapping manifest header. Expected tab-delimited columns: "
                + ", ".join(FIELDS)
            )
        seen = set()
        n = 0
        for line_no, row in enumerate(reader, start=2):
            n += 1
            assembler = (row.get("Assembler") or "").strip()
            if not assembler or not all(c.isalnum() or c in "._-" for c in assembler):
                raise ValueError(f"Invalid assembler name at manifest line {line_no}: {assembler!r}")
            if assembler in seen:
                raise ValueError(f"Duplicate assembler at manifest line {line_no}: {assembler}")
            seen.add(assembler)
            if None in row or any(not (row.get(f) or "").strip() for f in FIELDS):
                raise ValueError(f"Missing manifest field at line {line_no}. Re-run 01_mapping.sh")
            source = Path(row["Assembly_FASTA"])
            bam = Path(row["HiFi_all_BAM"])
            check_file(source, f"Assembly FASTA for {assembler}")
            check_file(bam, f"HiFi_all_BAM for {assembler} (run 01_mapping.sh first)")
            check_file(bam.parent / ".cache" / (bam.name + ".inputs.sha256"), f"HiFi_all_BAM cache for {assembler} (run 01_mapping.sh first)")
    if n == 0:
        raise ValueError("Mapping manifest has no assembler rows. Run 01_mapping.sh first")
    return len(mags), n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--mag-dir", default="")
    parser.add_argument("--manifest", default="")
    args = parser.parse_args()
    try:
        mag_count, assembler_count = check_inputs(
            Path(args.assembly), Path(args.mag_dir) if args.mag_dir else None,
            Path(args.manifest) if args.manifest else None,
        )
    except ValueError as error:
        parser.exit(1, f"ERROR: {error}\n")
    print(f"Clipping inputs found: {mag_count} MAG FASTAs; {assembler_count} assembler(s).")
    print("File presence checked. BAM integrity and contig identities are verified during analysis.")


if __name__ == "__main__":
    main()
