#!/usr/bin/env python3
"""Gene taxonomy for the chimera module: gene selection, prefilter and MMseqs2 result tables.

Subcommands:
  select             Select contigs or genes for a MAG/non-MAG scope
  filter             Keep genes on contigs with enough predicted genes
  mmseqs             Join MMseqs2 taxonomy with gene coordinates
  merge              Merge the per-scope gene taxonomy tables
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, TextIO, Tuple


ENTITY_COLUMNS = (
    "mag_id", "entity_id", "genome_id", "bin_id", "genome", "bin", "id"
)

CONTIG_COLUMNS = ("contig", "original_contig", "contig_id", "sequence", "seqid")


def choose_column(fieldnames: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {name.strip().lower(): name for name in fieldnames}
    return next((lookup[name] for name in candidates if name in lookup), None)


def load_map(path: Path) -> Dict[str, str]:
    ownership: Dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: empty MAG contig map: {path}")
        entity_col = choose_column(reader.fieldnames, ENTITY_COLUMNS)
        contig_col = choose_column(reader.fieldnames, CONTIG_COLUMNS)
        if entity_col is None or contig_col is None:
            raise SystemExit(
                "ERROR: MAG contig map needs a MAG/entity column and a contig column"
            )
        for line_no, row in enumerate(reader, 2):
            entity = row[entity_col].strip()
            contig = row[contig_col].strip()
            if not entity or not contig:
                raise SystemExit(
                    f"ERROR: blank MAG/contig value in map line {line_no}: {path}"
                )
            previous = ownership.get(contig)
            if previous is not None and previous != entity:
                raise SystemExit(
                    f"ERROR: contig assigned to multiple MAGs: {contig} ({previous}, {entity})"
                )
            ownership[contig] = entity
    if not ownership:
        raise SystemExit(f"ERROR: MAG contig map contains no contigs: {path}")
    return ownership


def fasta_records(path: Path) -> Iterator[Tuple[str, str, list[str]]]:
    header: str | None = None
    name: str | None = None
    sequence: list[str] = []
    observed: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if header is not None and name is not None:
                    yield name, header, sequence
                header = raw.rstrip("\n\r")
                text = header[1:].strip()
                name = text.split(None, 1)[0] if text else ""
                if not name:
                    raise SystemExit(f"ERROR: blank FASTA identifier at line {line_no}")
                if name in observed:
                    raise SystemExit(f"ERROR: duplicate FASTA identifier: {name}")
                observed.add(name)
                sequence = []
            else:
                if header is None:
                    if raw.strip():
                        raise SystemExit("ERROR: sequence data found before first FASTA header")
                    continue
                sequence.append(raw.rstrip("\n\r"))
    if header is not None and name is not None:
        yield name, header, sequence


def fasta_main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("include", "exclude"))
    parser.add_argument("assembly_fasta", type=Path)
    parser.add_argument("mag_contig_map", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("stats_json", type=Path)
    args = parser.parse_args(argv)

    ownership = load_map(args.mag_contig_map)
    found_mapped: set[str] = set()
    n_input = 0
    n_output = 0
    bp_input = 0
    bp_output = 0

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.output_fasta.open("w", encoding="utf-8", newline="") as out:
        for name, header, sequence_lines in fasta_records(args.assembly_fasta):
            n_input += 1
            length = sum(len(line.strip()) for line in sequence_lines)
            bp_input += length
            is_mag = name in ownership
            if is_mag:
                found_mapped.add(name)
            keep = is_mag if args.mode == "include" else not is_mag
            if not keep:
                continue
            n_output += 1
            bp_output += length
            out.write(header + "\n")
            for line in sequence_lines:
                out.write(line + "\n")

    missing = sorted(set(ownership) - found_mapped)
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise SystemExit(
            f"ERROR: {len(missing)} mapped contigs are absent from assembly: {preview}{suffix}"
        )

    stats = {
        "mode": args.mode,
        "input_contigs": n_input,
        "input_bp": bp_input,
        "mapped_contigs": len(ownership),
        "output_contigs": n_output,
        "output_bp": bp_output,
    }
    args.stats_json.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(stats, sort_keys=True))


def protein_records(path: Path) -> Iterator[Tuple[str, str, list[str]]]:
    header: str | None = None
    gene_id: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if header is not None and gene_id is not None:
                    yield gene_id, header, sequence
                header = raw.rstrip("\n\r")
                text = header[1:].strip()
                gene_id = text.split(None, 1)[0] if text else ""
                if not gene_id:
                    raise SystemExit(
                        f"ERROR: blank protein identifier in {path}:{line_no}"
                    )
                sequence = []
            else:
                if header is None:
                    if raw.strip():
                        raise SystemExit(
                            f"ERROR: protein sequence before first header in {path}:{line_no}"
                        )
                    continue
                sequence.append(raw.rstrip("\n\r"))
    if header is not None and gene_id is not None:
        yield gene_id, header, sequence


def write_protein(out: TextIO, header: str, sequence: list[str]) -> None:
    out.write(header + "\n")
    for line in sequence:
        out.write(line + "\n")


def genes_main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("include", "exclude"))
    parser.add_argument("assembly_fasta", type=Path)
    parser.add_argument("mag_contig_map", type=Path)
    parser.add_argument("all_genes_faa", type=Path)
    parser.add_argument("all_gene_coordinates", type=Path)
    parser.add_argument("output_genes_faa", type=Path)
    parser.add_argument("output_gene_coordinates", type=Path)
    parser.add_argument("stats_json", type=Path)
    args = parser.parse_args(argv)

    ownership = load_map(args.mag_contig_map)
    assembly_contigs: set[str] = set()
    n_input_contigs = n_output_contigs = 0
    bp_input = bp_output = 0
    for name, _header, sequence_lines in fasta_records(args.assembly_fasta):
        assembly_contigs.add(name)
        n_input_contigs += 1
        length = sum(len(line.strip()) for line in sequence_lines)
        bp_input += length
        is_mag = name in ownership
        keep = is_mag if args.mode == "include" else not is_mag
        if keep:
            n_output_contigs += 1
            bp_output += length

    missing = sorted(set(ownership) - assembly_contigs)
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise SystemExit(
            f"ERROR: {len(missing)} mapped contigs are absent from assembly: "
            f"{preview}{suffix}"
        )

    args.output_genes_faa.parent.mkdir(parents=True, exist_ok=True)
    args.output_gene_coordinates.parent.mkdir(parents=True, exist_ok=True)
    proteins = protein_records(args.all_genes_faa)
    n_input_genes = n_output_genes = 0
    with (
        args.all_gene_coordinates.open(encoding="utf-8") as coordinates,
        args.output_genes_faa.open("w", encoding="utf-8", newline="") as faa_out,
        args.output_gene_coordinates.open(
            "w", encoding="utf-8", newline=""
        ) as coordinate_out,
    ):
        for line_no, raw in enumerate(coordinates, 1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) != 4:
                raise SystemExit(
                    f"ERROR: invalid gene-coordinate row in "
                    f"{args.all_gene_coordinates}:{line_no}"
                )
            gene_id, contig, start_text, end_text = fields
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as error:
                raise SystemExit(
                    f"ERROR: invalid gene coordinates in "
                    f"{args.all_gene_coordinates}:{line_no}"
                ) from error
            if not gene_id or not contig or start < 1 or end < start:
                raise SystemExit(
                    f"ERROR: invalid gene-coordinate values in "
                    f"{args.all_gene_coordinates}:{line_no}"
                )

            try:
                faa_gene_id, header, sequence = next(proteins)
            except StopIteration as error:
                raise SystemExit(
                    "ERROR: gene-coordinate table contains more records than protein FASTA"
                ) from error
            if faa_gene_id != gene_id:
                raise SystemExit(
                    "ERROR: protein FASTA and coordinate order disagree at record "
                    f"{line_no}: {faa_gene_id!r} versus {gene_id!r}"
                )

            n_input_genes += 1
            is_mag = contig in ownership
            keep = is_mag if args.mode == "include" else not is_mag
            if keep:
                write_protein(faa_out, header, sequence)
                coordinate_out.write(raw if raw.endswith("\n") else raw + "\n")
                n_output_genes += 1

        try:
            extra_gene_id, _header, _sequence = next(proteins)
        except StopIteration:
            pass
        else:
            raise SystemExit(
                "ERROR: protein FASTA contains records absent from the coordinate table: "
                f"first extra gene={extra_gene_id}"
            )

    stats = {
        "mode": args.mode,
        "input_contigs": n_input_contigs,
        "input_bp": bp_input,
        "mapped_contigs": len(ownership),
        "output_contigs": n_output_contigs,
        "output_bp": bp_output,
        "input_genes": n_input_genes,
        "output_genes": n_output_genes,
    }
    args.stats_json.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(stats, sort_keys=True))


def _select_main():
    import sys
    if len(sys.argv) < 2 or sys.argv[1] not in ("fasta", "genes"):
        raise SystemExit("Usage: chimera_gene_taxonomy.py select {fasta|genes} ARGS...")
    (fasta_main if sys.argv[1] == "fasta" else genes_main)(sys.argv[2:])


def cmd_select(argv):
    """Select contigs or genes for a MAG/non-MAG scope"""
    sys.argv = [sys.argv[0] + " select", *argv]
    return _select_main()


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def count_genes_by_contig(coord_path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with coord_path.open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) != 4:
                fail(f"malformed coordinate row {line_no}: expected 4 columns")
            gene, contig, start, end = fields
            if not gene or not contig:
                fail(f"blank gene/contig in coordinate row {line_no}")
            try:
                start_i, end_i = int(start), int(end)
            except ValueError as exc:
                fail(f"invalid coordinates in row {line_no}: {start}-{end}")
            if start_i < 1 or end_i < start_i:
                fail(f"invalid coordinates in row {line_no}: {start}-{end}")
            counts[contig] += 1
    if not counts:
        fail(f"coordinate table is empty: {coord_path}")
    return counts


def next_coordinate(handle, line_no_ref: list[int]):
    for raw in handle:
        line_no_ref[0] += 1
        if not raw.strip():
            continue
        fields = raw.rstrip("\n\r").split("\t")
        if len(fields) != 4:
            fail(f"malformed coordinate row {line_no_ref[0]} during streaming")
        return fields
    return None


def filter_faa(faa: Path, coords: Path, output: Path, min_genes: int) -> dict[str, int]:
    if min_genes < 1:
        fail("min_genes must be >=1")
    counts = count_genes_by_contig(coords)
    eligible = {contig for contig, n in counts.items() if n >= min_genes}

    output.parent.mkdir(parents=True, exist_ok=True)
    total_genes = selected_genes = 0
    total_contigs = len(counts)
    selected_contigs = len(eligible)

    with coords.open(encoding="utf-8") as coord_handle, faa.open(encoding="utf-8") as faa_handle, output.open("w", encoding="utf-8", newline="") as out:
        coord_line_no = [0]
        write_record = False
        for faa_line_no, raw in enumerate(faa_handle, 1):
            if raw.startswith(">"):
                fields = next_coordinate(coord_handle, coord_line_no)
                if fields is None:
                    fail(f"FASTA contains more records than coordinate table; first extra FASTA record at line {faa_line_no}")
                gene, contig, _start, _end = fields
                accession = raw[1:].split(None, 1)[0].strip()
                if accession != gene:
                    fail(
                        "FASTA/coordinate accession order mismatch: "
                        f"FASTA={accession!r}, coordinates={gene!r} at coordinate row {coord_line_no[0]}"
                    )
                total_genes += 1
                write_record = contig in eligible
                if write_record:
                    selected_genes += 1
                    out.write(raw)
            elif write_record:
                out.write(raw)

        extra = next_coordinate(coord_handle, coord_line_no)
        if extra is not None:
            fail(f"coordinate table contains more genes than FASTA; first extra gene={extra[0]!r}")

    if total_genes != sum(counts.values()):
        fail(f"FASTA/coordinate gene counts disagree: {total_genes} versus {sum(counts.values())}")

    return {
        "min_genes_per_contig": min_genes,
        "total_contigs": total_contigs,
        "selected_contigs": selected_contigs,
        "skipped_contigs": total_contigs - selected_contigs,
        "total_genes": total_genes,
        "selected_genes": selected_genes,
        "skipped_genes": total_genes - selected_genes,
    }


def _filter_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("faa", type=Path)
    parser.add_argument("coordinates", type=Path)
    parser.add_argument("output_faa", type=Path)
    parser.add_argument("stats_json", type=Path)
    parser.add_argument("min_genes", type=int)
    args = parser.parse_args()

    stats = filter_faa(args.faa, args.coordinates, args.output_faa, args.min_genes)
    args.stats_json.parent.mkdir(parents=True, exist_ok=True)
    args.stats_json.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    print(
        "Taxonomy prefilter: "
        f"{stats['selected_contigs']}/{stats['total_contigs']} contigs and "
        f"{stats['selected_genes']}/{stats['total_genes']} genes retained "
        f"(minimum {stats['min_genes_per_contig']} predicted genes per contig)"
    )


def cmd_filter(argv):
    """Keep genes on contigs with enough predicted genes"""
    sys.argv = [sys.argv[0] + " filter", *argv]
    return _filter_main()


RANKS = ("domain", "phylum", "class", "order", "family")

UNCLASSIFIED = {
    "", "0", "1", "na", "n/a", "none", "unknown", "unclassified", "no rank",
}

NCBI_DOMAIN_TAXIDS = {
    "2": "Bacteria",
    "2157": "Archaea",
    "2759": "Eukaryota",
}

RANK_ALIASES = {
    "domain": "domain",
    "superkingdom": "domain",
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
}

INSERT_BATCH_SIZE = 50_000

LOOKUP_BATCH_SIZE = 800


def clean_rank(value: str, rank: str) -> str:
    value = value.strip()
    if value.casefold() in UNCLASSIFIED:
        return "Unclassified"

    low = value.casefold()
    prefixes = {
        "domain": ("d__", "d_", "sk__", "sk_", "domain:", "superkingdom:"),
        "phylum": ("p__", "p_", "phylum:"),
        "class": ("c__", "c_", "class:"),
        "order": ("o__", "o_", "order:"),
        "family": ("f__", "f_", "family:"),
    }
    for prefix in prefixes.get(rank, ()):
        if low.startswith(prefix):
            value = value[len(prefix):].strip()
            low = value.casefold()
            break

    if value.casefold() in UNCLASSIFIED or low.startswith("uc_"):
        return "Unclassified"

    if rank == "domain":
        if value in NCBI_DOMAIN_TAXIDS:
            return NCBI_DOMAIN_TAXIDS[value]
        if "bacter" in low:
            return "Bacteria"
        if "archaea" in low or "archae" in low:
            return "Archaea"
        if any(token in low for token in ("eukary", "protist")):
            return "Eukaryota"

    if value.isdigit():
        return f"taxid:{value}"
    return value


def empty_ranks() -> dict[str, str]:
    return {rank: "Unclassified" for rank in RANKS}


def split_rank_field(value: str) -> list[str] | None:
    """Parse the fifth column produced by ``--lca-ranks``."""
    text = value.strip()
    if not text:
        return [""] * len(RANKS)

    for delimiter in (";", ","):
        parts = [part.strip() for part in text.split(delimiter)]
        if len(parts) == len(RANKS):
            return parts

    labelled: dict[str, str] = {}
    rank_pattern = "|".join(RANK_ALIASES)
    for token in re.split(r"[;,]", text):
        token = token.strip()
        match = re.match(
            rf"^({rank_pattern})(?:__|_|=)(.+)$",
            token,
            flags=re.IGNORECASE,
        )
        if match:
            rank = RANK_ALIASES[match.group(1).casefold()]
            labelled[rank] = match.group(2).strip()
    if labelled:
        return [labelled.get(rank, "") for rank in RANKS]
    return None


def parse_taxid_lineage_field(value: str) -> dict[str, str] | None:
    """Parse an MMseqs2 ``--tax-lineage 2`` full-lineage taxid column."""
    text = value.strip()
    if not text or re.fullmatch(r"[0-9]+(?:;[0-9]+)*", text) is None:
        return None

    ranks = empty_ranks()
    for taxid in text.split(";"):
        domain = NCBI_DOMAIN_TAXIDS.get(taxid)
        if domain is not None:
            ranks["domain"] = domain
            break
    return ranks


def parse_name_lineage_field(value: str) -> dict[str, str] | None:
    """Parse a backward-compatible ``--tax-lineage 1`` names column."""
    text = value.strip()
    if not text:
        return None

    ranks = empty_ranks()
    recognized = False
    prefix_to_rank = {
        "d": "domain",
        "sk": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
    }

    for token in text.split(";"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r"^(sk|d|p|c|o|f)(?:__|_)(.+)$", token, re.IGNORECASE)
        if match is None:
            continue
        recognized = True
        rank = prefix_to_rank[match.group(1).casefold()]
        name = clean_rank(match.group(2), rank)
        if name != "Unclassified":
            ranks[rank] = name

    if not recognized:
        return None
    return ranks


def standard_lca_fallback(fields: list[str]) -> dict[str, str]:
    """Recover a requested rank from standard LCA taxid/rank/name fields."""
    ranks = empty_ranks()
    if len(fields) < 4:
        return ranks
    taxid = fields[1].strip()
    rank_text = fields[2].strip().casefold()
    name = fields[3].strip()
    if taxid.casefold() in UNCLASSIFIED or name.casefold() in UNCLASSIFIED:
        return ranks
    target_rank = RANK_ALIASES.get(rank_text)
    if target_rank is not None:
        ranks[target_rank] = clean_rank(name, target_rank)
    return ranks


def merge_rank_sources(
    preferred: dict[str, str], fallback: dict[str, str]
) -> dict[str, str]:
    """Retain preferred values and fill unresolved ranks from fallback."""
    merged = empty_ranks()
    for rank in RANKS:
        primary = preferred.get(rank, "Unclassified")
        secondary = fallback.get(rank, "Unclassified")
        merged[rank] = primary if primary != "Unclassified" else secondary
    return merged


def parse_rank_values(fields: list[str], line_number: int) -> tuple[dict[str, str], bool]:
    """Return domain, phylum, class, order and family from an MMseqs2 TSV."""
    if len(fields) < 4:
        raise SystemExit(
            "ERROR: MMseqs2 taxonomy TSV has fewer than four columns at "
            f"line {line_number}; observed {len(fields)} columns"
        )

    ranks = standard_lca_fallback(fields)
    fallback_was_classified = any(value != "Unclassified" for value in ranks.values())
    found_rank_column = False
    found_lineage_column = False
    nonblank_additional = False

    for value in fields[4:]:
        if not value.strip():
            continue
        nonblank_additional = True

        taxid_lineage = parse_taxid_lineage_field(value)
        if taxid_lineage is not None:
            ranks = merge_rank_sources(ranks, taxid_lineage)
            found_lineage_column = True
            continue

        name_lineage = parse_name_lineage_field(value)
        if name_lineage is not None:
            ranks = merge_rank_sources(ranks, name_lineage)
            found_lineage_column = True
            continue

        resolved = split_rank_field(value)
        if resolved is not None:
            explicit = {
                rank: clean_rank(rank_value, rank)
                for rank, rank_value in zip(RANKS, resolved, strict=True)
            }
            ranks = merge_rank_sources(explicit, ranks)
            found_rank_column = True
            continue

        raise SystemExit(
            "ERROR: could not recognize an additional MMseqs2 taxonomy column "
            f"at line {line_number}: {value!r}. Expected five values from "
            "--lca-ranks superkingdom,phylum,class,order,family or a "
            "--tax-lineage 1/2 field."
        )

    if nonblank_additional and not found_rank_column and not found_lineage_column:
        raise SystemExit(
            "ERROR: no recognized MMseqs2 taxonomy extension column at "
            f"line {line_number}"
        )

    classified = any(value != "Unclassified" for value in ranks.values())
    used_standard_fallback = (
        fallback_was_classified
        and classified
        and not found_rank_column
        and not found_lineage_column
    )
    return ranks, used_standard_fallback


def read_coordinate(line: str, line_number: int) -> tuple[str, str, int, int]:
    fields = line.rstrip("\n\r").split("\t")
    if len(fields) != 4:
        raise SystemExit(f"ERROR: malformed coordinate row {line_number}")
    gene, contig, start_text, end_text = fields
    if not gene or not contig:
        raise SystemExit(f"ERROR: blank gene/contig in coordinate row {line_number}")
    try:
        start, end = int(start_text), int(end_text)
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid coordinate row {line_number}") from exc
    if start < 1 or end < start:
        raise SystemExit(f"ERROR: invalid coordinate range at row {line_number}")
    return gene, contig, start, end


def read_taxonomy(
    line: str, line_number: int
) -> tuple[str, str, str, str, str, str, int]:
    fields = line.rstrip("\n\r").split("\t")
    if not fields or not fields[0]:
        raise SystemExit(f"ERROR: blank query ID in MMseqs2 taxonomy line {line_number}")
    ranks, used_fallback = parse_rank_values(fields, line_number)
    return (
        fields[0],
        ranks["domain"],
        ranks["phylum"],
        ranks["class"],
        ranks["order"],
        ranks["family"],
        int(used_fallback),
    )


def nonempty_lines(handle: TextIO) -> Iterable[tuple[int, str]]:
    for line_number, line in enumerate(handle, start=1):
        if line.strip("\n\r"):
            yield line_number, line


def configure_database(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")
    connection.execute("PRAGMA cache_size=-262144")
    connection.execute("PRAGMA page_size=32768")
    connection.execute(
        """
        CREATE TABLE taxonomy (
            gene TEXT PRIMARY KEY,
            domain TEXT NOT NULL,
            phylum TEXT NOT NULL,
            class_name TEXT NOT NULL,
            order_name TEXT NOT NULL,
            family_name TEXT NOT NULL,
            used_fallback INTEGER NOT NULL
        ) WITHOUT ROWID
        """
    )


def insert_taxonomy_batch(
    connection: sqlite3.Connection,
    rows: list[tuple[str, str, str, str, str, str, int]],
) -> None:
    if not rows:
        return

    batch_seen: set[str] = set()
    for row in rows:
        gene = row[0]
        if gene in batch_seen:
            raise SystemExit(
                f"ERROR: duplicate query ID in MMseqs2 taxonomy TSV: {gene}"
            )
        batch_seen.add(gene)

    try:
        connection.executemany(
            "INSERT INTO taxonomy VALUES (?, ?, ?, ?, ?, ?, ?)", rows
        )
    except sqlite3.IntegrityError as exc:
        connection.rollback()
        for row in rows:
            present = connection.execute(
                "SELECT 1 FROM taxonomy WHERE gene = ?", (row[0],)
            ).fetchone()
            if present is not None:
                raise SystemExit(
                    f"ERROR: duplicate query ID in MMseqs2 taxonomy TSV: {row[0]}"
                ) from exc
        raise SystemExit(
            "ERROR: duplicate query ID in MMseqs2 taxonomy TSV"
        ) from exc


def load_taxonomy(connection: sqlite3.Connection, taxonomy_path: Path, allow_empty: bool = False) -> int:
    n_taxonomy = 0
    batch: list[tuple[str, str, str, str, str, str, int]] = []
    with taxonomy_path.open(encoding="utf-8") as taxonomy:
        for line_number, line in nonempty_lines(taxonomy):
            batch.append(read_taxonomy(line, line_number))
            n_taxonomy += 1
            if len(batch) >= INSERT_BATCH_SIZE:
                insert_taxonomy_batch(connection, batch)
                connection.commit()
                batch.clear()
        insert_taxonomy_batch(connection, batch)
        connection.commit()
    if n_taxonomy == 0 and not allow_empty:
        raise SystemExit("ERROR: MMseqs2 taxonomy TSV is empty")
    return n_taxonomy


def lookup_rows(
    connection: sqlite3.Connection,
    genes: list[str],
) -> dict[str, tuple[str, str, str, str, str, int]]:
    placeholders = ",".join("?" for _ in genes)
    query = (
        "SELECT gene, domain, phylum, class_name, order_name, family_name, "
        "used_fallback FROM taxonomy WHERE gene IN ("
        f"{placeholders})"
    )
    return {
        row[0]: (row[1], row[2], row[3], row[4], row[5], row[6])
        for row in connection.execute(query, genes)
    }


def delete_rows(connection: sqlite3.Connection, genes: list[str]) -> None:
    placeholders = ",".join("?" for _ in genes)
    connection.execute(f"DELETE FROM taxonomy WHERE gene IN ({placeholders})", genes)


def stream_coordinate_batches(
    handle: TextIO,
) -> Iterable[list[tuple[str, str, int, int, int]]]:
    batch: list[tuple[str, str, int, int, int]] = []
    for line_number, line in enumerate(handle, start=1):
        if not line.strip("\n\r"):
            continue
        gene, contig, start, end = read_coordinate(line, line_number)
        batch.append((gene, contig, start, end, line_number))
        if len(batch) >= LOOKUP_BATCH_SIZE:
            yield batch
            batch = []
    if batch:
        yield batch


def indexed_merge(
    coordinates_path: Path, taxonomy_path: Path, output_path: Path,
    allow_missing_as_unclassified: bool = False,
) -> None:
    n_coordinates = 0
    n_any = 0
    n_domain = 0
    n_family = 0
    n_unclassified = 0
    n_lca_fallback = 0
    n_missing_taxonomy = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    work_output = output_path.with_name(f".{output_path.name}.stream.{os.getpid()}")
    database_path = taxonomy_path.parent / f".taxonomy_merge.{os.getpid()}.sqlite3"

    connection: sqlite3.Connection | None = None
    try:
        connection = sqlite3.connect(database_path)
        configure_database(connection)
        n_taxonomy = load_taxonomy(
            connection, taxonomy_path, allow_empty=allow_missing_as_unclassified
        )

        with (
            coordinates_path.open(encoding="utf-8") as coordinates,
            work_output.open("w", encoding="utf-8", newline="") as out,
        ):
            out.write(
                "contig\tstart\tend\tdomain\tphylum\tclass\torder\tfamily\tgene_id\n"
            )

            for batch in stream_coordinate_batches(coordinates):
                genes = [row[0] for row in batch]
                seen: set[str] = set()
                duplicate: str | None = None
                for gene in genes:
                    if gene in seen:
                        duplicate = gene
                        break
                    seen.add(gene)
                if duplicate is not None:
                    raise SystemExit(
                        f"ERROR: duplicate gene coordinate accession: {duplicate}"
                    )

                taxonomy_by_gene = lookup_rows(connection, genes)
                if not allow_missing_as_unclassified:
                    for gene, _contig, _start, _end, line_number in batch:
                        if gene not in taxonomy_by_gene:
                            raise SystemExit(
                                "ERROR: gene ID from coordinate table is absent from the "
                                "remaining MMseqs2 taxonomy rows at coordinate line "
                                f"{line_number}: {gene!r}. This indicates a missing taxonomy "
                                "row, a duplicate coordinate ID, or inputs from different genes.faa files."
                            )

                for gene, contig, start, end, _line_number in batch:
                    tax = taxonomy_by_gene.get(gene)
                    if tax is None:
                        if not allow_missing_as_unclassified:
                            raise SystemExit(
                                f"ERROR: missing taxonomy row for coordinate gene: {gene}"
                            )
                        domain = phylum = class_name = order_name = family_name = "Unclassified"
                        used_fallback = 0
                        n_missing_taxonomy += 1
                    else:
                        (
                            domain,
                            phylum,
                            class_name,
                            order_name,
                            family_name,
                            used_fallback,
                        ) = tax
                    values = (domain, phylum, class_name, order_name, family_name)
                    classified = tuple(value != "Unclassified" for value in values)
                    n_coordinates += 1
                    n_any += int(any(classified))
                    n_domain += int(classified[0])
                    n_family += int(classified[4])
                    n_unclassified += int(not any(classified))
                    n_lca_fallback += used_fallback
                    out.write(
                        f"{contig}\t{start}\t{end}\t"
                        + "\t".join(values)
                        + f"\t{gene}\n"
                    )

                delete_rows(connection, genes)
                connection.commit()

        if n_coordinates == 0:
            raise SystemExit("ERROR: gene coordinate table is empty")

        extra = connection.execute("SELECT gene FROM taxonomy LIMIT 1").fetchone()
        if extra is not None:
            remaining = connection.execute("SELECT COUNT(*) FROM taxonomy").fetchone()[0]
            raise SystemExit(
                "ERROR: MMseqs2 taxonomy TSV contains query IDs absent from the "
                f"coordinate table; first extra query: {extra[0]!r}; remaining rows: {remaining}. "
                "The files may originate from different genes.faa files."
            )

        if n_taxonomy > 0 and (n_any == 0 or n_domain == 0):
            raise SystemExit(
                "ERROR: MMseqs2 produced no usable domain assignments. Verify the "
                "taxonomy database, --lca-ranks output and --tax-lineage 2 column."
            )

        os.replace(work_output, output_path)
    except BaseException:
        work_output.unlink(missing_ok=True)
        raise
    finally:
        if connection is not None:
            connection.close()
        database_path.unlink(missing_ok=True)
        Path(str(database_path) + "-shm").unlink(missing_ok=True)
        Path(str(database_path) + "-wal").unlink(missing_ok=True)
        Path(str(database_path) + "-journal").unlink(missing_ok=True)

    print(
        f"Genes: {n_coordinates}; MMseqs rows: {n_taxonomy}; "
        f"classified at any requested rank: {n_any}; classified at domain: {n_domain}; "
        f"classified at family: {n_family}; fully unclassified: {n_unclassified}; "
        f"classified rows recovered from standard LCA fields: {n_lca_fallback}; "
        f"prefilter-skipped genes written as Unclassified: {n_missing_taxonomy}; "
        "merge key: gene_id (disk-backed)"
    )


def _mmseqs_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("coordinates", type=Path)
    parser.add_argument("taxonomy", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--allow-missing-as-unclassified", action="store_true",
        help="Write coordinate genes absent from the taxonomy TSV as Unclassified. "
             "Used for lossless chimera prefiltering of contigs that cannot meet the minimum gene-count criterion.",
    )
    args = parser.parse_args()
    indexed_merge(
        args.coordinates, args.taxonomy, args.output,
        allow_missing_as_unclassified=args.allow_missing_as_unclassified,
    )


def cmd_mmseqs(argv):
    """Join MMseqs2 taxonomy with gene coordinates"""
    sys.argv = [sys.argv[0] + " mmseqs", *argv]
    return _mmseqs_main()


HEADER = [
    "contig", "start", "end", "domain", "phylum", "class", "order",
    "family", "gene_id",
]


def load(path: Path):
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != HEADER:
            raise SystemExit(
                f"ERROR: unexpected taxonomy header in {path}: {reader.fieldnames!r}; "
                f"expected {HEADER!r}"
            )
        for line_no, row in enumerate(reader, 2):
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError as error:
                raise SystemExit(f"ERROR: invalid coordinates in {path}:{line_no}") from error
            if not row["contig"] or not row["gene_id"] or start < 1 or end < start:
                raise SystemExit(f"ERROR: invalid taxonomy row in {path}:{line_no}")
            rows.append(
                (
                    row["contig"], start, end, row["domain"], row["phylum"],
                    row["class"], row["order"], row["family"], row["gene_id"],
                )
            )
    return rows


def _merge_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()

    rows = []
    genes = set()
    for path in args.inputs:
        for row in load(path):
            gene = row[8]
            if gene in genes:
                raise SystemExit(f"ERROR: duplicate gene_id across taxonomy tables: {gene}")
            genes.add(gene)
            rows.append(row)

    rows.sort(key=lambda row: (row[0], row[1], row[2], row[8]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as out:
        out.write("\t".join(HEADER) + "\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")
    temporary.replace(args.output)
    print(f"Merged {len(rows)} genes from {len(args.inputs)} tables: {args.output}")


def cmd_merge(argv):
    """Merge the per-scope gene taxonomy tables"""
    sys.argv = [sys.argv[0] + " merge", *argv]
    return _merge_main()


COMMANDS = {
    'select': cmd_select,
    'filter': cmd_filter,
    'mmseqs': cmd_mmseqs,
    'merge': cmd_merge,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
