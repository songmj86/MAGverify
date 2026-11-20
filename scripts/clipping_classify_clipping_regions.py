#!/usr/bin/env python3
"""Merge assembler-specific clipping tables and classify sites by contig region."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-manifest", required=True)
    parser.add_argument("--assemblies-fai", required=True)
    parser.add_argument("--terminal-window", type=int, default=100)
    parser.add_argument("--input", action="append", required=True,
                        help="assembler=TSV; may be supplied repeatedly")
    parser.add_argument("--all-output", required=True)
    parser.add_argument("--internal-output", required=True)
    parser.add_argument("--terminal-output", required=True)
    parser.add_argument("--positive-output")
    parser.add_argument("--internal-positive-output")
    parser.add_argument("--terminal-positive-output")
    return parser.parse_args()


def load_fai(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                fail(f"malformed FAI row {line_number}: {path}")
            name = fields[0]
            try:
                length = int(fields[1])
            except ValueError as error:
                fail(f"invalid FAI length at row {line_number}: {error}")
            if not name or length < 1 or name in lengths:
                fail(f"blank/duplicate FAI sequence ID: {name}")
            lengths[name] = length
    if not lengths:
        fail(f"FAI contains no sequences: {path}")
    return lengths


def load_sources(manifest: str) -> Tuple[Dict[str, str], Dict[str, Dict[str, int]]]:
    assembler_to_fasta: Dict[str, str] = {}
    source_lengths: Dict[str, Dict[str, int]] = {}
    with open(manifest, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Assembler", "Assembly_FASTA"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"mapping manifest lacks columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            assembler = (row.get("Assembler") or "").strip()
            fasta = (row.get("Assembly_FASTA") or "").strip()
            if not assembler or not fasta or assembler in assembler_to_fasta:
                fail(f"invalid/duplicate mapping manifest row {line_number}")
            fai = fasta + ".fai"
            if not os.path.isfile(fai) or os.path.getsize(fai) == 0:
                fail(f"source FAI missing for {assembler}: {fai}")
            assembler_to_fasta[assembler] = fasta
            source_lengths[assembler] = load_fai(fai)
    if not assembler_to_fasta:
        fail("mapping manifest contains no assembler rows")
    return assembler_to_fasta, source_lengths


def parse_inputs(values: Sequence[str]) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            fail(f"--input must be assembler=TSV: {value}")
        assembler, path = value.split("=", 1)
        assembler, path = assembler.strip(), path.strip()
        if not assembler or not path or assembler in result:
            fail(f"invalid/duplicate --input: {value}")
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"input table missing or empty: {path}")
        result[assembler] = path
    return result


def atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False,
        dir=str(destination.parent), prefix=destination.name + ".", suffix=".tmp",
    )
    try:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.close()
        os.replace(handle.name, destination)
    finally:
        try:
            handle.close()
        except Exception:
            pass
        if os.path.exists(handle.name):
            os.unlink(handle.name)


def main() -> None:
    args = parse_args()
    if args.terminal_window < 1:
        fail("terminal-window must be >=1")
    if not os.path.isfile(args.mapping_manifest) or os.path.getsize(args.mapping_manifest) == 0:
        fail(f"mapping manifest missing: {args.mapping_manifest}")
    if not os.path.isfile(args.assemblies_fai) or os.path.getsize(args.assemblies_fai) == 0:
        fail(f"all-assemblies FAI missing: {args.assemblies_fai}")

    assembler_to_fasta, source_lengths = load_sources(args.mapping_manifest)
    inputs = parse_inputs(args.input)
    if set(inputs) != set(assembler_to_fasta):
        missing = sorted(set(assembler_to_fasta) - set(inputs))
        extra = sorted(set(inputs) - set(assembler_to_fasta))
        fail(f"input assemblers do not match mapping manifest; missing={missing}, extra={extra}")

    combined = load_fai(args.assemblies_fai)
    union: Dict[str, Tuple[str, int]] = {}
    for assembler, lengths in source_lengths.items():
        for contig, length in lengths.items():
            if contig in union:
                fail(f"contig ID occurs in multiple source assemblers: {contig}")
            union[contig] = (assembler, length)
    if set(union) != set(combined):
        missing = sorted(set(combined) - set(union))[:10]
        extra = sorted(set(union) - set(combined))[:10]
        fail(
            "source-assembly contig union differs from concatenated assembly FAI; "
            f"missing_from_sources={missing}, extra_in_sources={extra}"
        )
    for contig, length in combined.items():
        if union[contig][1] != length:
            fail(f"source/combined contig length mismatch: {contig}")

    base_fields: List[str] | None = None
    rows: List[Dict[str, str]] = []
    seen = set()
    for assembler in sorted(inputs):
        path = inputs[assembler]
        with open(path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            required = {"contig", "length", "anvio_pos_0based", "pos_1based"}
            missing = required - set(fields)
            if missing:
                fail(f"{path} lacks columns: {sorted(missing)}")
            if base_fields is None:
                base_fields = fields
            elif fields != base_fields:
                fail(f"input clipping schemas differ: {path}")
            for line_number, raw in enumerate(reader, start=2):
                contig = (raw.get("contig") or "").strip()
                try:
                    length = int(raw["length"])
                    pos0 = int(raw["anvio_pos_0based"])
                    pos1 = int(raw["pos_1based"])
                except (TypeError, ValueError) as error:
                    fail(f"invalid clipping coordinate at {path}:{line_number}: {error}")
                if contig not in source_lengths[assembler]:
                    fail(f"clipping contig is not in source assembler {assembler}: {contig}")
                if source_lengths[assembler][contig] != length or not 0 <= pos0 < length or pos1 != pos0 + 1:
                    fail(f"stale clipping coordinate: {assembler}:{contig}:{pos1}")
                key = (contig, pos0)
                if key in seen:
                    fail(f"duplicate clipping position across source tables: {contig}:{pos1}")
                seen.add(key)
                distance = min(pos0, length - 1 - pos0)
                region = "TERMINAL" if distance < args.terminal_window else "INTERNAL"
                row = {"source_assembler": assembler}
                row.update(raw)
                row["distance_to_nearest_end_bp"] = str(distance)
                row["region_type"] = region
                row["reference_boundary"] = "Y" if distance == 0 else "N"
                rows.append(row)

    if base_fields is None:
        fail("no clipping table schema was loaded")
    preferred = [
        "source_assembler", "contig", "length", "pos_1based",
        "region_type", "distance_to_nearest_end_bp", "reference_boundary",
        "Clipping_analyzable", "Clipping_analysis_reason",
        "Clipping_ge1", "Clipping", "Stringent_clipping",
        "default_candidate", "stringent_candidate",
        "validated_clipping", "clipping_confidence",
    ]
    all_fields = ["source_assembler", *base_fields,
                  "distance_to_nearest_end_bp", "region_type", "reference_boundary"]
    output_fields = []
    for field in preferred + all_fields:
        if field in all_fields and field not in output_fields:
            output_fields.append(field)
    rows.sort(key=lambda row: (row["source_assembler"], row["contig"], int(row["pos_1based"])))
    internal = [row for row in rows if row["region_type"] == "INTERNAL"]
    terminal = [row for row in rows if row["region_type"] == "TERMINAL"]

    atomic_write(args.all_output, output_fields, rows)
    atomic_write(args.internal_output, output_fields, internal)
    atomic_write(args.terminal_output, output_fields, terminal)

    if args.positive_output or args.internal_positive_output or args.terminal_positive_output:
        if "Clipping" not in output_fields:
            fail("positive outputs requested but input lacks Clipping column")
        positive = [row for row in rows if (row.get("Clipping") or "").strip().upper() == "Y"]
        internal_positive = [row for row in positive if row["region_type"] == "INTERNAL"]
        terminal_positive = [row for row in positive if row["region_type"] == "TERMINAL"]
        if args.positive_output:
            atomic_write(args.positive_output, output_fields, positive)
        if args.internal_positive_output:
            atomic_write(args.internal_positive_output, output_fields, internal_positive)
        if args.terminal_positive_output:
            atomic_write(args.terminal_positive_output, output_fields, terminal_positive)

    print(f"Merged clipping rows: {len(rows)}")
    print(f"Internal rows: {len(internal)}")
    print(f"Terminal rows: {len(terminal)}")


if __name__ == "__main__":
    main()
