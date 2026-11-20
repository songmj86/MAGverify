#!/usr/bin/env python3
"""Convert GTDB-Tk summaries to the MAGverify MAG taxonomy TSV."""
import argparse
import csv
from pathlib import Path
import sys

HEADERS = ("ID", "GTDB_Domain", "GTDB_Phylum", "GTDB_Class", "GTDB_Order")
RANKS = ("d__", "p__", "c__", "o__")
MAG_EXTENSIONS = (".fa", ".fna", ".fasta")


def mag_ids(mag_dir: Path) -> set[str]:
    if not mag_dir.is_dir():
        raise ValueError(f"MAG directory not found: {mag_dir}")
    ids = {p.stem for p in mag_dir.iterdir() if p.is_file() and p.suffix.lower() in MAG_EXTENSIONS}
    if not ids:
        raise ValueError(f"No MAG FASTAs with {MAG_EXTENSIONS} found in {mag_dir}")
    return ids


def convert(gtdbtk_dir: Path, output: Path, mag_dir: Path | None = None) -> int:
    """Write one row per MAG. With mag_dir, only the MAGs in that directory are written;
    GTDB-Tk rows for other genomes are ignored, and every MAG must be classified."""
    wanted = mag_ids(mag_dir) if mag_dir is not None else None
    paths = [gtdbtk_dir / f"gtdbtk.{domain}.summary.tsv" for domain in ("bac120", "ar53")]
    existing = [p for p in paths if p.is_file() and p.stat().st_size > 0]
    if not existing:
        raise ValueError(f"No GTDB-Tk summary TSV found in {gtdbtk_dir}")
    data = {}
    for path in existing:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = {"user_genome", "classification"} - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{path}: missing columns {sorted(missing)}")
            for n, row in enumerate(reader, start=2):
                identifier = (row.get("user_genome") or "").strip()
                if wanted is not None and identifier not in wanted:
                    continue
                classification = (row.get("classification") or "").strip()
                if not identifier or not classification:
                    raise ValueError(f"{path}:{n}: empty genome ID or classification")
                if identifier in data:
                    raise ValueError(f"Duplicate MAG ID in GTDB-Tk output: {identifier}")
                ranks = {}
                for raw in classification.split(";"):
                    token = raw.strip()
                    if len(token) >= 3 and token[:3] in RANKS:
                        if token[:3] in ranks:
                            raise ValueError(f"{path}:{n}: duplicate taxonomy rank {token[:3]}")
                        ranks[token[:3]] = token
                if ranks.get("d__") not in ("d__Bacteria", "d__Archaea"):
                    raise ValueError(f"{path}:{n}: unresolved or invalid domain for {identifier}")
                data[identifier] = [identifier] + [ranks.get(rank, rank) for rank in RANKS]
    if wanted is not None:
        absent = sorted(wanted - data.keys())
        if absent:
            raise ValueError(f"MAGs in -g without a GTDB-Tk classification: {absent[:10]} "
                             f"(total {len(absent)}). Check that MAG file names match user_genome.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADERS)
        for identifier in sorted(data):
            writer.writerow(data[identifier])
    return len(data)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("-h", action="help", help=argparse.SUPPRESS)
    parser.add_argument("--gtdbtk-dir", required=True, type=Path, metavar="DIR",
                        help="GTDB-Tk classify_wf output directory "
                             "(gtdbtk.bac120.summary.tsv, gtdbtk.ar53.summary.tsv)")
    parser.add_argument("-g", dest="mag_dir", required=True, type=Path, metavar="DIR",
                        help="Folder of MAGs (.fa)")
    parser.add_argument("-o", dest="output_dir", required=True, type=Path, metavar="DIR",
                        help="Output folder (table: DIR/MAG_taxonomy.tsv)")
    args = parser.parse_args()
    output = args.output_dir / "MAG_taxonomy.tsv"
    try:
        count = convert(args.gtdbtk_dir, output, args.mag_dir)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(f"Wrote {count} MAG taxonomic classifications to {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
