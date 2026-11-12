#!/usr/bin/env python3
"""Minimal helpers for clipped-read contig-connection analysis.

Alignment records are handled by QNAME so multiple records from the same read
are never counted as separate support. Structural cross-contig evidence uses
primary and supplementary alignments only and requires a split alignment
involving at least one supplementary segment. Secondary alignments are excluded.
"""
from __future__ import annotations

import csv
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, NoReturn, Optional, Sequence, Set, Tuple

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class Alignment:
    qname: str
    flag: int
    rname: str
    mapq: int
    cigar: str
    rstart: int
    rend: int
    left_clip: int
    right_clip: int
    matched_ref_bases: int

    @property
    def secondary(self) -> bool:
        return bool(self.flag & 256)

    @property
    def supplementary(self) -> bool:
        return bool(self.flag & 2048)


def fail(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


def atomic_tsv(path: str, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for line_no, row in enumerate(rows, start=2):
                if len(row) != len(header):
                    fail(
                        f"TSV schema mismatch at {path}:{line_no}: "
                        f"{len(row)} != {len(header)}"
                    )
                writer.writerow(row)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def stream_sam(command: Sequence[str]) -> Iterator[str]:
    process = subprocess.Popen(
        list(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="utf-8", errors="replace",
    )
    assert process.stdout is not None and process.stderr is not None
    for line in process.stdout:
        yield line
    process.stdout.close()
    stderr = process.stderr.read()
    rc = process.wait()
    if rc != 0:
        fail(f"command failed ({rc}): {' '.join(command)}\n{stderr}")


def cigar_operations(cigar: str) -> List[Tuple[int, str]]:
    ops = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if not ops or "".join(f"{n}{op}" for n, op in ops) != cigar:
        fail(f"unsupported or malformed CIGAR: {cigar}")
    return ops


def parse_alignment(line: str) -> Optional[Alignment]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11 or fields[2] == "*" or fields[5] == "*":
        return None
    flag = int(fields[1])
    ops = cigar_operations(fields[5])
    ref_span = sum(n for n, op in ops if op in "MDN=X")
    matched = sum(n for n, op in ops if op in "M=X")
    if ref_span <= 0 or matched <= 0:
        return None
    left_clip = 0
    for n, op in ops:
        if op not in "SH":
            break
        left_clip += n
    right_clip = 0
    for n, op in reversed(ops):
        if op not in "SH":
            break
        right_clip += n
    rstart = int(fields[3]) - 1
    rend = rstart + ref_span
    return Alignment(
        fields[0], flag, fields[2], int(fields[4]), fields[5],
        rstart, rend, left_clip, right_clip, matched,
    )


def representative(alignments: Sequence[Alignment]) -> Optional[Alignment]:
    if not alignments:
        return None
    return sorted(
        alignments,
        key=lambda a: (
            a.secondary,
            -a.matched_ref_bases,
            a.supplementary,
            a.rstart,
            a.rend,
            a.flag,
            a.cigar,
        ),
    )[0]


def load_alignments_for_qnames(
    bam: str, qnames: Set[str]
) -> Dict[str, List[Alignment]]:
    result: Dict[str, List[Alignment]] = defaultdict(list)
    if not qnames:
        return result
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".qnames") as qf:
        qpath = qf.name
        for qname in sorted(qnames):
            qf.write(qname + "\n")
    try:
        for line in stream_sam(
            ["samtools", "view", "-F", "4", "-q", "0", "-N", qpath, bam]
        ):
            aln = parse_alignment(line)
            if aln is not None:
                result[aln.qname].append(aln)
    finally:
        os.unlink(qpath)
    return result


def summarize_external_connection(
    source: Dict[str, Alignment],
    all_by_read: Dict[str, List[Alignment]],
    source_contig: str,
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Return primary/supplementary QNAMEs supporting cross-contig split mapping.

    Secondary alignments are excluded from both the source and external segments.
    At least one segment must be supplementary so that an alternative placement
    alone cannot be counted as a structural cross-contig connection.
    """
    supported: Set[str] = set()
    targets: Dict[str, Set[str]] = defaultdict(set)
    for qname, source_aln in source.items():
        if source_aln.secondary:
            continue
        for other in all_by_read.get(qname, []):
            if other.rname == source_contig or other.secondary:
                continue
            if not (source_aln.supplementary or other.supplementary):
                continue
            supported.add(qname)
            targets[other.rname].add(qname)
    return supported, targets


def target_count_text(targets: Dict[str, Set[str]]) -> str:
    if not targets:
        return "NA"
    return ";".join(
        f"{target}:{len(reads)}"
        for target, reads in sorted(
            targets.items(), key=lambda item: (-len(item[1]), item[0])
        )
    )


def dominant_target(
    targets: Dict[str, Set[str]], denominator: int
) -> Tuple[str, int, float]:
    """Audit-only summary; dominance is never a binary evidence gate."""
    if not targets:
        return "NA", 0, 0.0
    ranked = sorted(
        ((target, len(reads)) for target, reads in targets.items()),
        key=lambda item: (-item[1], item[0]),
    )
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        n = ranked[0][1]
        return "AMBIGUOUS", n, n / denominator if denominator else 0.0
    target, n = ranked[0]
    return target, n, n / denominator if denominator else 0.0


def load_contig_to_mag(path: str) -> Dict[str, str]:
    with open(path, encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        lower = {field.lower(): field for field in fields}
        mag_col = lower.get("mag_id") or lower.get("id")
        contig_col = lower.get("contig") or lower.get("original_contig")
        if not mag_col or not contig_col:
            fail(f"MAG contig map lacks MAG_ID/Contig columns: {path}")
        return {
            row[contig_col].strip(): row[mag_col].strip()
            for row in reader
            if row.get(contig_col, "").strip()
        }
