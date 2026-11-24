#!/usr/bin/env python3
"""Circularity final tables for all assembler groups.

Subcommands:
  merge              Merge the assembler-group junction tables
  refresh            Apply assembler-reported circularity and write circularity.MAG.tsv and its audit
  add-junction       Append the junction read counts to circularity.MAG.tsv
  contig             Write circularity.contig.tsv: one row per assembly contig
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, NoReturn, Optional, Sequence, Set, Tuple


def read_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def cmd_merge(argv):
    """Merge the assembler-group junction tables"""
    plan_path, groups_path, grouped_root, validation_path, output_path, group_outputs_path = argv

    _, plan_rows = read_tsv(plan_path)

    expected = {
        row["MAG_ID"]: {
            "group": row["Group_ID"], "assembler": row["Assembler"],
            "assembly": row["Assembly_FASTA"], "contig": row["Contig_ID"],
        }
        for row in plan_rows
    }

    _, group_rows = read_tsv(groups_path)

    merged_rows = []

    output_fields = None

    seen = set()

    group_output_rows = []

    for group in group_rows:
        group_id = group["Group_ID"]
        final_path = os.path.join(grouped_root, group_id, "circularity.MAG.junction.tsv")
        if not os.path.isfile(final_path) or os.path.getsize(final_path) == 0:
            raise SystemExit(f"ERROR: Missing assembler-group result: {final_path}")
        fields, rows = read_tsv(final_path)
        if output_fields is None:
            output_fields = fields
        elif fields != output_fields:
            raise SystemExit(f"ERROR: Assembler-group final schema mismatch: {final_path}")
        group_output_rows.append((group_id, group["Assembler"], group["Assembly_FASTA"], group["MAG_count"], final_path))
        for row in rows:
            mag_id = row.get("MAG_ID", "")
            contig_id = row.get("Contig_ID", "")
            if mag_id not in expected:
                raise SystemExit(f"ERROR: Unexpected MAG_ID in {final_path}: {mag_id}")
            meta = expected[mag_id]
            if meta["group"] != group_id or meta["contig"] != contig_id:
                raise SystemExit(f"ERROR: Assembler provenance mismatch: {mag_id} / {contig_id}")
            if mag_id in seen:
                raise SystemExit(f"ERROR: MAG_ID occurs in multiple groups: {mag_id}")
            seen.add(mag_id)
            row["Assembler"] = meta["assembler"]
            row["Assembly_FASTA"] = meta["assembly"]
            merged_rows.append(row)

    if output_fields is None:
        raise SystemExit("ERROR: no assembler-group circularity outputs")

    for field in ("Assembler", "Assembly_FASTA"):
        if field not in output_fields:
            output_fields.append(field)

    _, validation_rows = read_tsv(validation_path)

    for item in validation_rows:
        mag_id = item["MAG_ID"]
        n = int(item["Number_of_contigs"])
        if n == 1:
            if mag_id not in seen:
                raise SystemExit(f"ERROR: Single-contig MAG missing from merged evidence: {mag_id}")
            continue
        row = {field: "NA" for field in output_fields}
        row.update({
            "MAG_ID": mag_id,
            "Number_of_contigs": str(n),
            "Contig_ID": "NA",
            "Contig_length_bp": "NA",
            "SINGLE_CONTIG": "N",
            "CIRCULARITY_ANALYZABLE": "N",
            "HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N": "N",
            "HIFI_ALL_MAPPED_ANALYZABLE_Y_N": "N",
            "HIFI_END_START_SUPPORTED": "NA",
            "CIRCULARITY_SUPPORTED": "NA",
            "Assembler": "NA",
            "Assembly_FASTA": "NA",
        })
        merged_rows.append(row)

    merged_rows.sort(key=lambda r: r["MAG_ID"])

    with open(output_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(merged_rows)

    with open(group_outputs_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(("Group_ID", "Assembler", "Assembly_FASTA", "MAG_count", "HiFi_final_result"))
        writer.writerows(group_output_rows)


def fail(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


TRUE_VALUES = {"Y", "YES", "TRUE", "T", "1", "CIRCULAR", "CIRCULAR=Y"}

FALSE_VALUES = {"N", "NO", "FALSE", "F", "0", "LINEAR", "NOT_CIRCULAR"}

ID_COLUMNS = ("id", "mag_id", "genome", "genome_id", "bin", "bin_id")

CALL_COLUMNS = ("assembler_circularity", "circularity", "circular", "is_circular")


def canonical_id(value: str) -> str:
    value = str(value).strip().strip('"').strip("'")
    value = re.split(r"[/\\]", value)[-1]
    lower = value.lower()
    for suffix in (".fasta", ".fna", ".fa"):
        if lower.endswith(suffix):
            return value[: -len(suffix)]
    return value


def normalize_call(value: str) -> Optional[str]:
    upper = str(value).strip().upper()
    if upper in TRUE_VALUES:
        return "Y"
    if upper in FALSE_VALUES:
        return "N"
    return None


def table_header(path: str) -> Optional[Tuple[list[str], int, int]]:
    try:
        with open(path, encoding="utf-8-sig") as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                header = raw.rstrip("\n\r").split("\t")
                lowered = [field.strip().lower() for field in header]
                id_index = next((lowered.index(x) for x in ID_COLUMNS if x in lowered), None)
                call_index = next((lowered.index(x) for x in CALL_COLUMNS if x in lowered), None)
                if id_index is None:
                    return None
                if call_index is not None:
                    return header, id_index, call_index
                if "assembler" not in lowered:
                    return None
                assembler_index = lowered.index("assembler")
                saw_binary = False
                for data_raw in handle:
                    if not data_raw.strip() or data_raw.startswith("#"):
                        continue
                    fields = data_raw.rstrip("\n\r").split("\t")
                    if len(fields) <= assembler_index or not fields[assembler_index].strip():
                        continue
                    if normalize_call(fields[assembler_index]) is None:
                        return None
                    saw_binary = True
                return (header, id_index, assembler_index) if saw_binary else None
    except (OSError, UnicodeDecodeError):
        return None
    return None


def parse_call_table(path: str) -> Dict[str, str]:
    header_info = table_header(path)
    if header_info is None:
        fail(f"no recognized ID/circularity header in {path}")
    _header, id_index, call_index = header_info
    calls: Dict[str, str] = {}
    header_seen = False
    with open(path, encoding="utf-8-sig") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if not header_seen:
                header_seen = True
                continue
            if len(fields) <= max(id_index, call_index):
                continue
            mag_id = canonical_id(fields[id_index])
            call = normalize_call(fields[call_index])
            if not mag_id or call is None:
                continue
            if mag_id in calls and calls[mag_id] != call:
                fail(f"conflicting circularity calls for {mag_id} in {path}")
            calls[mag_id] = call
    return calls


def resolve_call_table(path: str, mag_ids: Set[str]) -> Tuple[str, Dict[str, str]]:
    candidate_path = Path(path)
    if candidate_path.is_file():
        return str(candidate_path), parse_call_table(str(candidate_path))
    if not candidate_path.is_dir():
        fail(f"CIRCULARITY_TSV is neither a file nor directory: {path}")
    scored: list[Tuple[int, int, str, Dict[str, str]]] = []
    for candidate in sorted(candidate_path.rglob("*")):
        if not candidate.is_file():
            continue
        if table_header(str(candidate)) is None:
            continue
        calls = parse_call_table(str(candidate))
        scored.append((len(mag_ids & set(calls)), len(calls), str(candidate), calls))
    if not scored:
        fail(f"no table with recognized circularity calls under {path}")
    scored.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
    best_overlap, best_count = scored[0][0], scored[0][1]
    if best_overlap == 0:
        fail(f"circularity tables under {path} do not overlap target MAG IDs")
    tied = [x for x in scored if x[0] == best_overlap and x[1] == best_count]
    if len(tied) > 1:
        fail("multiple circularity tables tie for best MAG-ID overlap")
    return scored[0][2], scored[0][3]


def _refresh_read_tsv(path: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not fields:
        fail(f"empty TSV header: {path}")
    return fields, rows


def atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        prefix=os.path.basename(path) + ".", suffix=".tmp", dir=directory, text=True
    )
    os.close(fd)
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def load_assembler_calls(circularity_source: str, mag_ids: set[str]) -> Dict[str, str]:
    _selected_call_table, assembler_calls = resolve_call_table(circularity_source, mag_ids)
    missing = sorted(mag_ids - set(assembler_calls))
    if missing:
        preview = ", ".join(missing[:10])
        fail(
            f"assembler circularity call is missing for {len(missing)} MAGs; "
            f"examples: {preview}"
        )
    invalid = {mag: call for mag, call in assembler_calls.items() if mag in mag_ids and call not in {"Y", "N"}}
    if invalid:
        fail(f"invalid assembler circularity calls: {invalid}")
    return assembler_calls


def build_rows(
    canonical_rows: Sequence[Mapping[str, str]], assembler_calls: Mapping[str, str]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    final_rows: List[dict[str, str]] = []
    audit_rows: List[dict[str, str]] = []
    seen: set[str] = set()

    for row in canonical_rows:
        mag = row["MAG_ID"]
        if mag in seen:
            fail(f"duplicate MAG_ID in canonical circularity output: {mag}")
        seen.add(mag)

        try:
            n = int(row["Number_of_contigs"])
        except ValueError:
            fail(f"invalid Number_of_contigs for {mag}: {row['Number_of_contigs']}")
        if n < 1:
            fail(f"Number_of_contigs must be >=1 for {mag}: {n}")

        eligibility = "Y_single_contig" if n == 1 else "N_not_single_contig"
        canonical_analyzable = row.get("CIRCULARITY_ANALYZABLE", "")
        canonical_topology = row.get("HIFI_END_START_SUPPORTED", "")

        if canonical_analyzable not in {"Y", "N"}:
            fail(f"invalid CIRCULARITY_ANALYZABLE for {mag}: {canonical_analyzable}")
        if canonical_topology not in {"Y", "N", "NA"}:
            fail(f"invalid HIFI_END_START_SUPPORTED for {mag}: {canonical_topology}")

        if n != 1:
            if canonical_analyzable != "N" or canonical_topology != "NA":
                fail(
                    f"multi-contig MAG must be non-analyzable with topology NA: {mag} "
                    f"({canonical_analyzable}, {canonical_topology})"
                )
            evaluable = "N"
            topology = "NA"
        else:
            evaluable = canonical_analyzable
            if evaluable == "Y":
                if canonical_topology not in {"Y", "N"}:
                    fail(f"evaluable single-contig MAG lacks Y/N topology call: {mag}")
                topology = canonical_topology
            else:
                if canonical_topology != "NA":
                    fail(f"non-evaluable single-contig MAG must have topology NA: {mag}")
                topology = "NA"

        assembler = assembler_calls[mag]
        if assembler == "NA":
            high_confidence = additional_candidate = "NA"
            circular = topology
        else:
            high_confidence = "Y" if assembler == "Y" and topology == "Y" else "N"
            additional_candidate = "Y" if assembler == "N" and topology == "Y" else "N"
            circular = high_confidence

        final_rows.append(
            {
                "ID": mag,
                "Assembler_reported_circularity": assembler,
                "Circularity_evaluable": evaluable,
                "Circularity_eligibility": eligibility,
                "Read_supported_circular_topology": topology,
                "High_confidence_circularity": high_confidence,
                "Additional_circular_candidate": additional_candidate,
                "Circular": circular,
            }
        )
        audit_rows.append(
            {"ID": mag, "Contig_count": str(n), "Eligibility": eligibility}
        )

    return final_rows, audit_rows


def validate_expected_structure(rows: Sequence[Mapping[str, str]]) -> None:
    high_conf = Counter(row["High_confidence_circularity"] for row in rows)
    candidates = Counter(row["Additional_circular_candidate"] for row in rows)
    both = [row["ID"] for row in rows if row["High_confidence_circularity"] == "Y" and row["Additional_circular_candidate"] == "Y"]
    if both:
        fail(f"MAG cannot be both high-confidence circular and additional candidate: {both[:10]}")
    if high_conf.get("Y", 0) + candidates.get("Y", 0) > sum(row["Read_supported_circular_topology"] == "Y" for row in rows):
        fail("derived circular classes exceed read-supported topology calls")


def _refresh_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical", required=True, help="circularity.MAG.junction.tsv")
    parser.add_argument("--assembler-circularity", required=True,
                        help="Assembler circularity TSV or directory (prep-circularity output)")
    parser.add_argument("--final", required=True, help="circularity.MAG.tsv")
    parser.add_argument("--audit", help="Optional circularity.MAG.audit.tsv")
    args = parser.parse_args()

    fields, canonical_rows = _refresh_read_tsv(args.canonical)
    required = {
        "MAG_ID", "Number_of_contigs", "CIRCULARITY_ANALYZABLE",
        "HIFI_END_START_SUPPORTED",
    }
    missing = required - set(fields)
    if missing:
        fail(f"canonical circularity output lacks {sorted(missing)}: {args.canonical}")
    if not canonical_rows:
        fail(f"canonical circularity output has no data rows: {args.canonical}")

    mag_ids = {row["MAG_ID"] for row in canonical_rows}
    if len(mag_ids) != len(canonical_rows):
        fail("canonical circularity output contains duplicate MAG_ID values")

    assembler_calls = load_assembler_calls(args.assembler_circularity, mag_ids)
    final_rows, audit_rows = build_rows(canonical_rows, assembler_calls)
    validate_expected_structure(final_rows)

    final_fields = [
        "ID", "Assembler_reported_circularity", "Circularity_evaluable",
        "Circularity_eligibility", "Read_supported_circular_topology",
        "High_confidence_circularity", "Additional_circular_candidate", "Circular",
    ]
    atomic_write(args.final, final_fields, final_rows)
    if args.audit:
        atomic_write(args.audit, ["ID", "Contig_count", "Eligibility"], audit_rows)

    counts = {
        "total": len(final_rows),
        "single_contig": sum(row["Circularity_eligibility"] == "Y_single_contig" for row in final_rows),
        "evaluable": sum(row["Circularity_evaluable"] == "Y" for row in final_rows),
        "topology_Y": sum(row["Read_supported_circular_topology"] == "Y" for row in final_rows),
        "topology_N": sum(row["Read_supported_circular_topology"] == "N" for row in final_rows),
        "topology_NA": sum(row["Read_supported_circular_topology"] == "NA" for row in final_rows),
        "assembler_Y": sum(row["Assembler_reported_circularity"] == "Y" for row in final_rows),
        "high_confidence": sum(row["High_confidence_circularity"] == "Y" for row in final_rows),
        "additional_candidates": sum(row["Additional_circular_candidate"] == "Y" for row in final_rows),
    }
    print("Circularity compatibility output refreshed:")
    for key, value in counts.items():
        print(f"  {key}: {value}")


def cmd_refresh(argv):
    """Apply assembler-reported circularity and write circularity.MAG.tsv and its audit"""
    sys.argv = [sys.argv[0] + " refresh", *argv]
    return _refresh_main()


def read(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader.fieldnames or []), list(reader)


def _add_junction_main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    final_path, junction_path = sys.argv[1:]
    final_fields, final_rows = read(final_path)
    junction_fields, junction_rows = read(junction_path)
    if 'ID' not in final_fields or 'MAG_ID' not in junction_fields:
        sys.exit('ERROR: circularity tables lack ID / MAG_ID columns')
    added = [field for field in junction_fields if field != 'MAG_ID']
    base = [field for field in final_fields if field not in added]
    junction = {row['MAG_ID']: row for row in junction_rows}
    missing = [row['ID'] for row in final_rows if row['ID'] not in junction]
    if missing:
        sys.exit(f'ERROR: MAGs absent from the junction table: {missing[:5]}')
    temp = f'{final_path}.tmp.{os.getpid()}'
    with open(temp, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(base + added)
        for row in final_rows:
            extra = junction[row['ID']]
            writer.writerow([row[field] for field in base] + [extra.get(field, 'NA') for field in added])
    os.replace(temp, final_path)


def cmd_add_junction(argv):
    """Append the junction read counts to circularity.MAG.tsv"""
    sys.argv = [sys.argv[0] + " add-junction", *argv]
    return _add_junction_main()


FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fas", ".fa.gz", ".fna.gz", ".fasta.gz", ".fas.gz")
MAG_LEVEL_FIELDS = ("ID", "Assembler_reported_circularity", "Number_of_contigs", "Contig_ID",
                    "Contig_length_bp", "MAG_ID")


def strip_fasta_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in sorted(FASTA_SUFFIXES, key=len, reverse=True):
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return name


def contig_lengths(path: Path) -> List[Tuple[str, int]]:
    """Contig IDs and lengths, from an up-to-date .fai or by reading the FASTA."""
    fai = Path(str(path) + ".fai")
    if fai.is_file() and fai.stat().st_mtime_ns >= path.stat().st_mtime_ns:
        with fai.open(encoding="utf-8") as handle:
            return [(line.split("\t")[0], int(line.split("\t")[1])) for line in handle if line.strip()]
    import gzip
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    out: List[Tuple[str, int]] = []
    name, length = None, 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, length))
                name, length = line[1:].split()[0], 0
            else:
                length += len(line.strip())
    if name is not None:
        out.append((name, length))
    return out


def mag_membership(mag_dir: Path) -> Dict[str, str]:
    import gzip
    membership: Dict[str, str] = {}
    for path in sorted(mag_dir.iterdir()):
        if not (path.is_file() and path.name.lower().endswith(FASTA_SUFFIXES)):
            continue
        mag = strip_fasta_suffix(path.name)
        opener = gzip.open if path.name.lower().endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    contig = line[1:].split()[0]
                    if membership.get(contig, mag) != mag:
                        fail(f"contig {contig} is in two MAGs: {membership[contig]} and {mag}")
                    membership[contig] = mag
    return membership


def cmd_contig(argv):
    """Write circularity.contig.tsv: one row per assembly contig"""
    parser = argparse.ArgumentParser(prog="circularity_finalize.py contig")
    parser.add_argument("--mag-table", required=True, help="circularity.MAG.tsv (with junction columns)")
    parser.add_argument("--mag-dir", required=True)
    parser.add_argument("--assembly-dir", required=True)
    parser.add_argument("--assembly-fastas", required=True, help="Space-separated FASTA basenames")
    parser.add_argument("--assembler-circularity", default="",
                        help="assembler_circularity.contig.tsv (optional)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    fields, mag_rows = _refresh_read_tsv(args.mag_table)
    by_mag = {row["ID"]: row for row in mag_rows}
    evidence = [field for field in fields if field not in MAG_LEVEL_FIELDS]
    membership = mag_membership(Path(args.mag_dir))
    reported: Dict[str, str] = {}
    if args.assembler_circularity and os.path.isfile(args.assembler_circularity):
        _fields, rows = _refresh_read_tsv(args.assembler_circularity)
        for row in rows:
            if "Contig" in row and "Assembler_circularity" in row:
                reported[row["Contig"]] = row["Assembler_circularity"]

    out_fields = ["Contig", "MAG_ID", "Assembler", "Contig_length_bp", "MAG_n_contigs",
                  "Assembler_reported_circularity"] + evidence
    rows_out: List[dict[str, str]] = []
    seen: Set[str] = set()
    for name in args.assembly_fastas.split():
        assembler = strip_fasta_suffix(Path(name).name)
        for contig, length in contig_lengths(Path(args.assembly_dir) / name):
            if contig in seen:
                fail(f"contig ID appears in more than one assembler FASTA: {contig}")
            seen.add(contig)
            mag = membership.get(contig, "NA")
            mag_row = by_mag.get(mag)
            row = {field: "NA" for field in out_fields}
            row.update(Contig=contig, MAG_ID=mag, Assembler=assembler, Contig_length_bp=str(length))
            row["Assembler_reported_circularity"] = reported.get(contig, "NA")
            if mag_row is None:
                row["Circularity_evaluable"] = "N"
                row["Circularity_eligibility"] = "N_unbinned" if mag == "NA" else "NA"
            else:
                row["MAG_n_contigs"] = mag_row.get("Number_of_contigs", "NA")
                if mag_row.get("Contig_ID") == contig:
                    for field in evidence:
                        row[field] = mag_row.get(field, "NA")
                    if contig not in reported:
                        row["Assembler_reported_circularity"] = mag_row.get(
                            "Assembler_reported_circularity", "NA")
                else:
                    row["Circularity_evaluable"] = "N"
                    row["Circularity_eligibility"] = mag_row.get("Circularity_eligibility", "NA")
            rows_out.append(row)
    missing = sorted(set(membership) - seen)
    if missing:
        fail(f"{len(missing)} MAG contigs are absent from the assembler FASTAs; examples: {missing[:5]}")
    atomic_write(args.output, out_fields, rows_out)
    print(f"Circularity contig table: {args.output} ({len(rows_out)} contigs)")

COMMANDS = {
    'merge': cmd_merge,
    'refresh': cmd_refresh,
    'add-junction': cmd_add_junction,
    'contig': cmd_contig,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
