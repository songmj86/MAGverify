#!/usr/bin/env python3
"""Build the assembler-reported circularity tables from native assembler outputs.

Accepted files (format detected from the content, gzip allowed):
  hifiasm-meta  p_ctg GFA, or FASTA/header list; contig names end in c (circular) or l (linear)
  metaMDBG      contigs FASTA or header list with circular=yes/no
  metaFlye      assembly_info.txt (circ. column)

Outputs in OUT/00_prep/:
  assembler_circularity.contig.tsv  Contig, Assembler, Assembler_circularity, File
  assembler_circularity.MAG.tsv     MAG_ID, Assembler, Assembler_circularity, Number_of_contigs, Contig
The source files are read in parallel, one process per file.
A MAG is Y when it consists of one contig that its assembler reported as circular.
"""

from __future__ import annotations

import shutil
import argparse
import gzip
import os
import re
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

HIFIASM_NAME = re.compile(r"ctg\d+([cl])$")
MAG_SUFFIXES = (".fa", ".fna", ".fasta")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def open_text(path: Path):
    with path.open("rb") as handle:
        gzipped = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gzipped else path.open()


def open_binary(path: Path):
    with path.open("rb") as handle:
        gzipped = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if gzipped else path.open("rb", buffering=1 << 20)


def extract_lines(path: Path, prefix: bytes) -> Iterator[str]:
    """Yield lines starting with prefix, scanning the file as bytes."""
    with open_binary(path) as handle:
        for raw in handle:
            if raw.startswith(prefix):
                yield raw.decode("utf-8", "replace").rstrip("\r\n")


def first_line(path: Path) -> str:
    with open_text(path) as handle:
        for line in handle:
            if line.strip():
                return line
    fail(f"empty assembler circularity file: {path}")
    return ""


def hifiasm_call(name: str, path: Path) -> str:
    match = HIFIASM_NAME.search(name)
    if not match:
        fail(f"{path}: contig {name} has no hifiasm c/l suffix")
    return "Y" if match.group(1) == "c" else "N"


def read_flye(path: Path) -> Iterator[Tuple[str, str]]:
    with open_text(path) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        try:
            name_i, circ_i = header.index("seq_name"), header.index("circ.")
        except ValueError:
            fail(f"{path}: metaFlye assembly_info.txt needs seq_name and circ. columns")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > max(name_i, circ_i) and fields[name_i]:
                call = fields[circ_i].strip().upper()
                if call not in {"Y", "N"}:
                    fail(f"{path}: invalid circ. value for {fields[name_i]}: {call}")
                yield fields[name_i], call


def read_gfa(path: Path) -> Iterator[Tuple[str, str]]:
    for line in extract_lines(path, b"S\t"):
        name = line.split("\t", 2)[1].strip()
        yield name, hifiasm_call(name, path)


def read_headers(path: Path) -> Iterator[Tuple[str, str]]:
    for line in extract_lines(path, b">"):
        fields = line[1:].split()
        if fields:
            name = fields[0]
            tag = next((f for f in fields[1:] if f.lower().startswith("circular=")), None)
            if tag is not None:
                value = tag.split("=", 1)[1].lower()
                if value not in {"yes", "no", "y", "n", "true", "false"}:
                    fail(f"{path}: invalid {tag} for {name}")
                yield name, "Y" if value in {"yes", "y", "true"} else "N"
            else:
                yield name, hifiasm_call(name, path)


def read_source(path: Path) -> Tuple[str, str, Iterator[Tuple[str, str]]]:
    """Return (assembler, file kind, records)."""
    line = first_line(path)
    if line.startswith("#seq_name") or "\tcirc." in line:
        return "metaFlye", "assembly_info.txt", read_flye(path)
    if line.startswith(("H\t", "S\t", "L\t")):
        return "hifiasm-meta", "GFA", read_gfa(path)
    if line.startswith(">"):
        if "circular=" in line:
            return "metaMDBG", "FASTA", read_headers(path)
        return "hifiasm-meta", "FASTA", read_headers(path)
    fail(f"{path}: not a hifiasm GFA/FASTA, metaMDBG FASTA or metaFlye assembly_info.txt")
    return "", "", iter(())


def load_source(path: Path) -> Tuple[str, str, List[Tuple[str, str]]]:
    assembler, kind, records = read_source(path)
    return assembler, kind, list(records)


def mag_files(mag_dir: Path) -> List[Tuple[str, Path]]:
    found = []
    for path in sorted(mag_dir.iterdir()):
        lower = path.name.lower()
        suffix = next((s for s in MAG_SUFFIXES if lower.endswith(s)), None)
        if path.is_file() and suffix:
            found.append((path.name[: -len(suffix)], path))
    if not found:
        fail(f"no MAG FASTA files in {mag_dir}")
    return found


def mag_contigs(path: Path) -> List[str]:
    with open_text(path) as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">") and line[1:].strip()]


def write_atomic(path: Path, header: str, rows: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    with os.fdopen(fd, "w") as out:
        out.write(header + "\n")
        out.writelines(row + "\n" for row in rows)
    os.replace(temp, path)



class OptionHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Show options without value placeholders, each option on one line."""

    def __init__(self, prog):
        width = max(100, shutil.get_terminal_size((120, 24)).columns)
        super().__init__(prog, max_help_position=40, width=width)

    def _format_args(self, action, default_metavar):
        return '' if action.option_strings else super()._format_args(action, default_metavar)

    def _format_action_invocation(self, action):
        return ', '.join(action.option_strings) or super()._format_action_invocation(action)

    def _format_usage(self, usage, actions, groups, prefix):
        text = super()._format_usage(usage, actions, groups, prefix)
        text = re.sub(r'(?<=\S)  +(?=\S)', ' ', text.replace(' ]', ']'))
        return '\n'.join(line.rstrip() for line in text.split('\n'))

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, add_help=False,
                                     formatter_class=OptionHelpFormatter)
    parser.add_argument("-h", action="help", help=argparse.SUPPRESS)
    parser.add_argument("-c", dest="sources", action="append", required=True, type=Path, metavar="FILE",
                        help="Assembler circularity file (repeat per assembler)")
    parser.add_argument("-g", dest="mag_dir", required=True, type=Path, metavar="DIR", help="Directory of MAGs (*.fa)")
    parser.add_argument("-o", dest="output", required=True, type=Path, metavar="DIR", help="Directory of base output")
    args = parser.parse_args()

    for path in args.sources:
        if not path.is_file():
            fail(f"assembler circularity file not found: {path}")

    with ProcessPoolExecutor(max_workers=len(args.sources)) as pool:
        loaded = list(pool.map(load_source, args.sources))

    calls: Dict[str, Tuple[str, str, str]] = {}
    for path, (assembler, kind, records) in zip(args.sources, loaded):
        circular = 0
        for name, call in records:
            previous = calls.get(name)
            if previous and previous[0] != call:
                fail(f"conflicting circularity for {name}: {previous[2]} vs {path.name}")
            calls[name] = (call, assembler, path.name)
            circular += call == "Y"
        print(f"{assembler} ({path.name}, {kind}): {len(records)} contigs, {circular} circular", flush=True)

    mags = mag_files(args.mag_dir)
    mag_rows, missing = [], []
    for mag_id, path in mags:
        contigs = mag_contigs(path)
        if not contigs:
            fail(f"MAG FASTA without sequences: {path}")
        absent = [c for c in contigs if c not in calls]
        if absent:
            missing.extend(f"{mag_id}:{c}" for c in absent)
            continue
        call = "Y" if len(contigs) == 1 and calls[contigs[0]][0] == "Y" else "N"
        assemblers = ",".join(sorted({calls[c][1] for c in contigs}))
        contig = contigs[0] if len(contigs) == 1 else "NA"
        mag_rows.append(f"{mag_id}\t{assemblers}\t{call}\t{len(contigs)}\t{contig}")
    if missing:
        fail(f"{len(missing)} MAG contigs are absent from the assembler circularity files, e.g. "
             f"{', '.join(missing[:5])}. Give the complete assembler outputs with -c "
             f"(the contig names must match the MAG FASTAs)")

    out_dir = args.output / "00_prep"
    contig_rows = [f"{name}\t{assembler}\t{call}\t{source}"
                   for name, (call, assembler, source) in sorted(calls.items())]
    write_atomic(out_dir / "assembler_circularity.contig.tsv",
                 "Contig\tAssembler\tAssembler_circularity\tFile", contig_rows)
    write_atomic(out_dir / "assembler_circularity.MAG.tsv",
                 "MAG_ID\tAssembler\tAssembler_circularity\tNumber_of_contigs\tContig", mag_rows)
    n_circular = sum(row.split("\t")[2] == "Y" for row in mag_rows)
    print(f"MAGs: {len(mag_rows)}, assembler-reported circular: {n_circular}")
    print(f"MAG table: {out_dir / 'assembler_circularity.MAG.tsv'}")
    print(f"Contig table: {out_dir / 'assembler_circularity.contig.tsv'}")


if __name__ == "__main__":
    main()
