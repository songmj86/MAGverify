#!/usr/bin/env python3
"""Merge MMseqs2 rank-resolved taxonomy with exact gene coordinates by ID.

MMseqs2 ``createtsv`` does not guarantee that query rows are emitted in the
same order as the original protein FASTA. The merge therefore uses the query
accession rather than row position.

The requested MMseqs2 ranks are ``superkingdom,phylum,class,order,family``.
The downstream table renames ``superkingdom`` to ``domain``. To remain robust
to taxonomy dumps in which the three cellular domains are not represented at
the ``superkingdom`` rank, ``--tax-lineage 2`` is also parsed and domain is
recovered from NCBI taxids 2, 2157 or 2759.

For tens of millions of genes, taxonomy rows are indexed in a temporary
on-disk SQLite database. Coordinate rows are processed in bounded batches,
looked up by ``gene_id`` and written in coordinate-table order. Used taxonomy
rows are deleted, allowing strict detection of duplicate coordinate IDs and
extra taxonomy rows without holding all accessions in memory.
"""

from __future__ import annotations

import argparse
import os
import re
import sqlite3
from pathlib import Path
from typing import Iterable, TextIO

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
LOOKUP_BATCH_SIZE = 800  # safely below SQLite's conservative variable limit


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
    """Parse the fifth column produced by ``--lca-ranks``.

    MMseqs2 v18 accepts ``superkingdom`` but not ``domain`` as a requested
    rank. Values are returned in the requested order as a semicolon-separated
    string. Comma-separated input is accepted for compatibility with synthetic
    or older intermediate files.
    """
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
    """Parse an MMseqs2 ``--tax-lineage 2`` full-lineage taxid column.

    Only domain is recovered from the lineage because the other requested
    ranks are supplied directly by ``--lca-ranks``. A one-token lineage such
    as ``0`` is accepted as a valid unclassified lineage field.
    """
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
    connection.execute("PRAGMA cache_size=-262144")  # up to ~256 MiB
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


def main() -> None:
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


if __name__ == "__main__":
    main()
