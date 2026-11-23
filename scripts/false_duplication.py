#!/usr/bin/env python3
"""False duplication: candidate pairs from self-alignment, HiFi validation and clipping co-occurrence.

Subcommands:
  detect             Classify duplicated copies and validate them with HiFi reads
  annotate           Report validated clipping at the duplication boundaries
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

MATCH_OPS = {"M", "=", "X"}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-paf", required=True)
    p.add_argument("--min-copy-bp", type=int, default=1000)
    p.add_argument("--min-pident", type=float, default=95.0)
    p.add_argument("--redundant-reciprocal-overlap", type=float, default=0.90)
    p.add_argument("--mapping-manifest", required=True)
    p.add_argument("--mag-map", required=True)
    p.add_argument("--hifi-anchor-bp", type=int, default=500)
    p.add_argument("--hifi-min-spanning-reads", type=int, default=2)
    p.add_argument("--hifi-min-reaching-reads", type=int, default=2)
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--out-all-pairs", required=True)
    p.add_argument("--out-all-boundaries", required=True)
    p.add_argument("--out-all-contigs", required=True)
    p.add_argument("--out-assemblers", required=True)
    p.add_argument("--out-pairs", required=True)
    p.add_argument("--out-boundaries", required=True)
    p.add_argument("--out-contigs", required=True)
    p.add_argument("--out-mags", required=True)
    p.add_argument("--out-idn-events", required=True)
    p.add_argument("--out-idn-summary", required=True)
    p.add_argument("--out-candidate-preprocessing", required=True)
    return p.parse_args()


def read_fai(path: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(path) as h:
        for raw in h:
            f = raw.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            name, length = f[0], int(f[1])
            if not name or length < 1 or name in out:
                fail(f"invalid or duplicate FAI record: {raw.rstrip()}")
            out[name] = length
    if not out:
        fail(f"empty FASTA index: {path}")
    return out


def read_mag_map(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path, newline="") as h:
        r = csv.DictReader(h, delimiter="\t")
        lower = {x.lower(): x for x in (r.fieldnames or [])}
        mag_col = lower.get("mag_id") or lower.get("id")
        contig_col = lower.get("contig") or lower.get("original_contig")
        if not mag_col or not contig_col:
            fail("MAG map requires MAG_ID and Contig columns")
        for row in r:
            mag, contig = row[mag_col].strip(), row[contig_col].strip()
            if not mag or not contig:
                continue
            if contig in out and out[contig] != mag:
                fail(f"contig assigned to multiple MAGs: {contig}")
            out[contig] = mag
    if not out:
        fail(f"MAG map contains no assignments: {path}")
    return out


@dataclass(frozen=True)
class MappingSource:
    assembler: str
    assembly_fasta: str
    lengths: Dict[str, int]
    hifi_bam: str


def read_mapping_manifest(path: str) -> List[MappingSource]:
    required = {"Assembler", "Assembly_FASTA", "HiFi_all_BAM"}
    out: List[MappingSource] = []
    with open(path, newline="") as h:
        r = csv.DictReader(h, delimiter="\t")
        fields = set(r.fieldnames or [])
        missing = required - fields
        if missing:
            fail(f"mapping manifest lacks columns: {sorted(missing)}")
        for row in r:
            assembler = row["Assembler"].strip()
            fasta = row["Assembly_FASTA"].strip()
            hifi = row["HiFi_all_BAM"].strip()
            if not assembler or not fasta or not hifi:
                fail("mapping manifest contains an incomplete assembler row")
            for req, label in ((fasta + ".fai", "FAI"), (hifi, "HiFi BAM")):
                if not os.path.isfile(req) or os.path.getsize(req) == 0:
                    fail(f"{assembler} {label} missing: {req}")
            out.append(MappingSource(assembler, fasta, read_fai(fasta + ".fai"), hifi))
    if not out:
        fail("mapping manifest contains no assembler rows")
    return out


def route_all_contigs(mag_map: Dict[str, str], sources: Sequence[MappingSource]) -> Dict[str, MappingSource]:
    routed: Dict[str, MappingSource] = {}
    for source in sources:
        for contig in source.lengths:
            if contig in routed:
                fail(f"contig occurs in multiple assembler-specific FASTAs: {contig}")
            routed[contig] = source
    for contig in mag_map:
        if contig not in routed:
            fail(f"MAG contig absent from assemblies: {contig}")
    return routed


def interval_overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def reciprocal_overlap_fraction(a1: int, a2: int, b1: int, b2: int) -> float:
    """Minimum reciprocal overlap fraction between two inclusive intervals."""
    ov = interval_overlap(a1, a2, b1, b2)
    if ov <= 0:
        return 0.0
    len_a = a2 - a1 + 1
    len_b = b2 - b1 + 1
    return min(ov / len_a, ov / len_b)


def interval_contains(outer1: int, outer2: int, inner1: int, inner2: int) -> bool:
    """Return True when the inclusive inner interval is fully contained in outer."""
    return outer1 <= inner1 and inner2 <= outer2


@dataclass(frozen=True)
class DuplicationPair:
    pair_id: str
    mag_id: str
    assembler: str
    pair_scope: str
    duplication_type: str
    orientation: str
    copy1_contig: str
    copy1_contig_length: int
    copy1_start: int
    copy1_end: int
    copy2_contig: str
    copy2_contig_length: int
    copy2_start: int
    copy2_end: int
    pident: float
    alignment_length: int
    candidate_cluster_size: int = 1

    @property
    def tandem(self) -> bool:
        return self.duplication_type == "TANDEM"

    def boundary_refs(self) -> Dict[str, Optional[Tuple[str, int]]]:
        def left(contig: str, start: int) -> Optional[Tuple[str, int]]:
            return None if start <= 1 else (contig, start)

        def right(contig: str, end: int, length: int) -> Optional[Tuple[str, int]]:
            return None if end >= length else (contig, end + 1)

        return {
            "copy1_left": left(self.copy1_contig, self.copy1_start),
            "copy1_right": right(self.copy1_contig, self.copy1_end, self.copy1_contig_length),
            "copy2_left": left(self.copy2_contig, self.copy2_start),
            "copy2_right": right(self.copy2_contig, self.copy2_end, self.copy2_contig_length),
        }


def open_text_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def canonical_candidate(
    q: str, qlen: int, qs: int, qe: int,
    t: str, tlen: int, ts: int, te: int,
    strand: str,
) -> Tuple[str, int, int, int, str, int, int, int, str]:
    """Canonicalize reciprocal self-alignment records."""
    orientation = "DIRECT" if strand == "+" else "INVERTED"
    if q == t:
        a = (qs, qe)
        b = (ts, te)
        if a <= b:
            return q, qlen, qs, qe, t, tlen, ts, te, orientation
        return t, tlen, ts, te, q, qlen, qs, qe, orientation
    if (q, qs, qe) <= (t, ts, te):
        return q, qlen, qs, qe, t, tlen, ts, te, orientation
    return t, tlen, ts, te, q, qlen, qs, qe, orientation


def read_duplication_pairs(
    path: str,
    mag_map: Dict[str, str],
    routed: Dict[str, MappingSource],
    min_copy_bp: int,
    min_pident: float,
    redundant_reciprocal_overlap: float,
) -> Tuple[List[DuplicationPair], List[Dict[str, object]]]:
    """Parse Minimap2 self-alignment PAF into discrete duplication candidates."""
    if not 0 < redundant_reciprocal_overlap <= 1:
        fail("redundant_reciprocal_overlap must be >0 and <=1")

    raw: Dict[Tuple[object, ...], Tuple[float, int]] = {}
    stats = defaultdict(int)
    with open_text_maybe_gzip(path) as h:
        for line in h:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                fail(f"malformed PAF row with {len(f)} fields")
            stats["rows"] += 1
            q, t = f[0], f[5]
            if q not in routed or t not in routed:
                fail(f"self-alignment contig absent from assemblies: {q},{t}")
            if routed[q].assembler != routed[t].assembler:
                fail(f"different-assembler self-alignment detected: {q},{t}")
            try:
                qlen = int(f[1]); qs0 = int(f[2]); qe0 = int(f[3])
                strand = f[4]
                tlen = int(f[6]); ts0 = int(f[7]); te0 = int(f[8])
                nmatch = int(f[9]); block_len = int(f[10])
            except ValueError as exc:
                fail(f"malformed numeric PAF field: {exc}")
            if strand not in {"+", "-"} or block_len <= 0:
                continue
            if qlen != routed[q].lengths[q] or tlen != routed[t].lengths[t]:
                fail(f"stale PAF contig length: {q},{t}")
            pident = 100.0 * nmatch / block_len
            if (qe0 - qs0) < min_copy_bp or (te0 - ts0) < min_copy_bp or pident < min_pident:
                continue

            qs, qe = qs0 + 1, qe0
            ts, te = ts0 + 1, te0
            if q == t and qs == ts and qe == te:
                continue

            c1, l1, s1, e1, c2, l2, s2, e2, orientation = canonical_candidate(
                q, qlen, qs, qe, t, tlen, ts, te, strand
            )
            if c1 != c2:
                mag1, mag2 = mag_map.get(c1), mag_map.get(c2)
                if not mag1 or not mag2 or mag1 != mag2:
                    stats["inter_contig_outside_same_mag_excluded"] += 1
                    continue

            key = (c1, s1, e1, c2, s2, e2, orientation)
            prev = raw.get(key)
            score = (pident, block_len)
            if prev is None or score > prev:
                raw[key] = score

    @dataclass(frozen=True)
    class Provisional:
        c1: str; l1: int; s1: int; e1: int
        c2: str; l2: int; s2: int; e2: int
        orientation: str; pident: float; alen: int
        mag: str; scope: str; dtype: str; assembler: str

        @property
        def estimated_matching_bp(self) -> float:
            return self.alen * self.pident / 100.0

        @property
        def total_genomic_span(self) -> int:
            return (self.e1 - self.s1 + 1) + (self.e2 - self.s2 + 1)

    provisional: List[Provisional] = []
    preprocessing_rows: List[Dict[str, object]] = []

    for key, (pident, alen) in raw.items():
        c1, s1, e1, c2, s2, e2, orientation = key
        source = routed[c1]
        mag1 = mag_map.get(c1, "NA")

        if c1 == c2:
            mag = mag1
            scope = "INTRA_CONTIG"
            if interval_overlap(s1, e1, s2, e2) > 0:
                stats["intrinsic_copy_overlap_excluded"] += 1
                preprocessing_rows.append({
                    "preprocessing_action": "EXCLUDED",
                    "reason": "INTRINSIC_COPY_INTERVAL_OVERLAP",
                    "assembler": source.assembler,
                    "MAG_ID": mag,
                    "PAIR_SCOPE": scope,
                    "DUPLICATION_TYPE": "SELF_OVERLAP",
                    "ORIENTATION": orientation,
                    "copy1_contig": c1,
                    "copy1_start": s1,
                    "copy1_end": e1,
                    "copy2_contig": c2,
                    "copy2_start": s2,
                    "copy2_end": e2,
                    "pident": f"{pident:.4f}",
                    "alignment_length": alen,
                    "representative_copy1_contig": "NA",
                    "representative_copy1_start": "NA",
                    "representative_copy1_end": "NA",
                    "representative_copy2_contig": "NA",
                    "representative_copy2_start": "NA",
                    "representative_copy2_end": "NA",
                    "copy1_reciprocal_overlap": "NA",
                    "copy2_reciprocal_overlap": "NA",
                    "redundant_reciprocal_overlap_threshold": f"{redundant_reciprocal_overlap:.4f}",
                })
                continue
            dtype = "TANDEM" if e1 + 1 == s2 else "INTERSPERSED"
        else:
            mag = mag1
            scope = "INTER_CONTIG"
            dtype = "INTER_CONTIG"

        provisional.append(Provisional(
            c1, source.lengths[c1], s1, e1,
            c2, routed[c2].lengths[c2], s2, e2,
            orientation, pident, alen, mag, scope, dtype, source.assembler,
        ))

    grouped: Dict[Tuple[str, str, str, str, str], List[Provisional]] = defaultdict(list)
    for p in provisional:
        grouped[(p.assembler, p.c1, p.c2, p.orientation, p.scope)].append(p)

    kept: List[Provisional] = []
    cluster_size: Dict[Provisional, int] = {}

    for group_key in sorted(grouped):
        members = grouped[group_key]
        ranked = sorted(
            members,
            key=lambda x: (
                -x.total_genomic_span,
                -x.estimated_matching_bp,
                -x.alen,
                -x.pident,
                x.s1, x.e1, x.s2, x.e2,
            ),
        )
        representatives: List[Provisional] = []

        for cand in ranked:
            matched_rep: Optional[Provisional] = None
            matched_ro1 = 0.0
            matched_ro2 = 0.0
            collapse_reason: Optional[str] = None
            for rep in representatives:
                ro1 = reciprocal_overlap_fraction(cand.s1, cand.e1, rep.s1, rep.e1)
                ro2 = reciprocal_overlap_fraction(cand.s2, cand.e2, rep.s2, rep.e2)

                nested_in_rep = (
                    interval_contains(rep.s1, rep.e1, cand.s1, cand.e1)
                    and interval_contains(rep.s2, rep.e2, cand.s2, cand.e2)
                )
                if nested_in_rep:
                    matched_rep = rep
                    matched_ro1 = ro1
                    matched_ro2 = ro2
                    collapse_reason = "REDUNDANT_NESTED_ALIGNMENT"
                    break

                if ro1 >= redundant_reciprocal_overlap and ro2 >= redundant_reciprocal_overlap:
                    matched_rep = rep
                    matched_ro1 = ro1
                    matched_ro2 = ro2
                    collapse_reason = "REDUNDANT_NEAR_DUPLICATE_ALIGNMENT"
                    break

            if matched_rep is None:
                representatives.append(cand)
                cluster_size[cand] = 1
                continue

            cluster_size[matched_rep] += 1
            stats["redundant_candidate_excluded"] += 1
            if collapse_reason == "REDUNDANT_NESTED_ALIGNMENT":
                stats["redundant_nested_excluded"] += 1
            else:
                stats["redundant_near_duplicate_excluded"] += 1
            preprocessing_rows.append({
                "preprocessing_action": "COLLAPSED_TO_REPRESENTATIVE",
                "reason": collapse_reason,
                "assembler": cand.assembler,
                "MAG_ID": cand.mag,
                "PAIR_SCOPE": cand.scope,
                "DUPLICATION_TYPE": cand.dtype,
                "ORIENTATION": cand.orientation,
                "copy1_contig": cand.c1,
                "copy1_start": cand.s1,
                "copy1_end": cand.e1,
                "copy2_contig": cand.c2,
                "copy2_start": cand.s2,
                "copy2_end": cand.e2,
                "pident": f"{cand.pident:.4f}",
                "alignment_length": cand.alen,
                "representative_copy1_contig": matched_rep.c1,
                "representative_copy1_start": matched_rep.s1,
                "representative_copy1_end": matched_rep.e1,
                "representative_copy2_contig": matched_rep.c2,
                "representative_copy2_start": matched_rep.s2,
                "representative_copy2_end": matched_rep.e2,
                "copy1_reciprocal_overlap": f"{matched_ro1:.4f}",
                "copy2_reciprocal_overlap": f"{matched_ro2:.4f}",
                "redundant_reciprocal_overlap_threshold": f"{redundant_reciprocal_overlap:.4f}",
            })

        kept.extend(representatives)

    kept.sort(key=lambda x: (x.assembler, x.c1, x.s1, x.c2, x.s2, x.orientation, x.dtype))
    pairs: List[DuplicationPair] = []
    for i, p in enumerate(kept, 1):
        pairs.append(DuplicationPair(
            pair_id=f"FD{i:08d}", mag_id=p.mag, assembler=p.assembler,
            pair_scope=p.scope, duplication_type=p.dtype, orientation=p.orientation,
            copy1_contig=p.c1, copy1_contig_length=p.l1, copy1_start=p.s1, copy1_end=p.e1,
            copy2_contig=p.c2, copy2_contig_length=p.l2, copy2_start=p.s2, copy2_end=p.e2,
            pident=p.pident, alignment_length=p.alen,
            candidate_cluster_size=cluster_size.get(p, 1),
        ))

    representative_ids = {
        (
            pair.assembler,
            pair.copy1_contig, pair.copy1_start, pair.copy1_end,
            pair.copy2_contig, pair.copy2_start, pair.copy2_end,
            pair.orientation,
        ): pair.pair_id
        for pair in pairs
    }
    for row in preprocessing_rows:
        if row["reason"] in {"REDUNDANT_NESTED_ALIGNMENT", "REDUNDANT_NEAR_DUPLICATE_ALIGNMENT"}:
            rep_key = (
                str(row["assembler"]),
                str(row["representative_copy1_contig"]),
                int(row["representative_copy1_start"]),
                int(row["representative_copy1_end"]),
                str(row["representative_copy2_contig"]),
                int(row["representative_copy2_start"]),
                int(row["representative_copy2_end"]),
                str(row["ORIENTATION"]),
            )
            row["representative_pair_id"] = representative_ids.get(rep_key, "NA")
        else:
            row["representative_pair_id"] = "NA"

    preprocessing_rows.sort(
        key=lambda r: (
            str(r["assembler"]), str(r["copy1_contig"]), int(r["copy1_start"]),
            str(r["copy2_contig"]), int(r["copy2_start"]), str(r["reason"]),
        )
    )

    print(
        "Duplication candidate preprocessing: "
        f"PAF_rows={stats['rows']}; exact_coordinate_candidates={len(raw)}; "
        f"intrinsic_copy_overlap_excluded={stats['intrinsic_copy_overlap_excluded']}; "
        f"redundant_nested_excluded={stats['redundant_nested_excluded']}; "
        f"redundant_near_duplicate_excluded={stats['redundant_near_duplicate_excluded']}; "
        f"redundant_total_excluded={stats['redundant_candidate_excluded']}; "
        f"final_discrete_candidate_pairs={len(pairs)}; "
        f"inter_contig_not_same_MAG_excluded={stats['inter_contig_outside_same_mag_excluded']}; "
        f"redundancy_reciprocal_overlap={redundant_reciprocal_overlap:.3f}"
    )
    return pairs, preprocessing_rows


@dataclass(frozen=True)
class CigarEvent:
    op: str
    length: int
    ref_start: int
    ref_end: int


def cigar_blocks_and_events(pos1: int, cigar: str) -> Tuple[int, List[Tuple[int, int]], List[CigarEvent]]:
    ref = pos1
    blocks: List[Tuple[int, int]] = []
    events: List[CigarEvent] = []
    current_start: Optional[int] = None
    current_end: Optional[int] = None
    for m in CIGAR_RE.finditer(cigar):
        n, op = int(m.group(1)), m.group(2)
        if op in MATCH_OPS:
            start, end = ref, ref + n - 1
            if current_start is not None and current_end is not None and start == current_end + 1:
                current_end = end
            else:
                if current_start is not None and current_end is not None:
                    blocks.append((current_start, current_end))
                current_start, current_end = start, end
            ref += n
        elif op == "I":
            events.append(CigarEvent("I", n, ref, ref))
        elif op in {"D", "N"}:
            if current_start is not None and current_end is not None:
                blocks.append((current_start, current_end))
            current_start = current_end = None
            events.append(CigarEvent(op, n, ref, ref + n - 1))
            ref += n
        elif op in {"S", "H", "P"}:
            if current_start is not None and current_end is not None:
                blocks.append((current_start, current_end))
            current_start = current_end = None
    if current_start is not None and current_end is not None:
        blocks.append((current_start, current_end))
    return ref - 1, blocks, events


def cigar_blocks(pos1: int, cigar: str) -> Tuple[int, List[Tuple[int, int]]]:
    ref_end, blocks, _ = cigar_blocks_and_events(pos1, cigar)
    return ref_end, blocks


def event_hits_boundary_window(event: CigarEvent, boundary: int, anchor: int) -> bool:
    ws, we = boundary - anchor, boundary + anchor - 1
    if event.op == "I":
        return ws < event.ref_start <= we
    return max(event.ref_start, ws) <= min(event.ref_end, we)


def alignment_type(flag: int) -> str:
    if flag & 0x100:
        return "secondary"
    if flag & 0x800:
        return "supplementary"
    return "primary"


def any_block_covers(blocks: Sequence[Tuple[int, int]], start: int, end: int) -> bool:
    return start <= end and any(b1 <= start and b2 >= end for b1, b2 in blocks)


def matched_reference_bases(blocks: Sequence[Tuple[int, int]], start: Optional[int] = None, end: Optional[int] = None) -> int:
    total = 0
    for b1, b2 in blocks:
        lo = b1 if start is None else max(b1, start)
        hi = b2 if end is None else min(b2, end)
        if lo <= hi:
            total += hi - lo + 1
    return total


def boundary_is_spanned(blocks: Sequence[Tuple[int, int]], boundary: int, anchor: int) -> bool:
    """Direct boundary spanning with a continuous M/=X reference block."""
    if boundary <= 1 or anchor < 1 or boundary - anchor < 1:
        return False
    return any_block_covers(blocks, boundary - anchor, boundary + anchor - 1)


def side_has_anchor(blocks: Sequence[Tuple[int, int]], boundary: int, anchor: int, side: str) -> bool:
    """Require a continuous M/=X anchor immediately adjacent to the boundary."""
    if anchor < 1:
        return False
    if side == "left":
        if boundary - anchor < 1:
            return False
        return any_block_covers(blocks, boundary - anchor, boundary - 1)
    if side == "right":
        return any_block_covers(blocks, boundary, boundary + anchor - 1)
    raise ValueError(side)


@dataclass
class BoundaryEvidence:
    spanning_reads: Set[str] = field(default_factory=set)
    left_reaching: Set[str] = field(default_factory=set)
    right_reaching: Set[str] = field(default_factory=set)
    i_qnames: Set[str] = field(default_factory=set)
    d_qnames: Set[str] = field(default_factory=set)
    n_qnames: Set[str] = field(default_factory=set)
    i_lengths: List[int] = field(default_factory=list)
    d_lengths: List[int] = field(default_factory=list)
    n_lengths: List[int] = field(default_factory=list)


def all_boundary_positions(pairs: Sequence[DuplicationPair]) -> Tuple[Dict[str, List[int]], Dict[Tuple[str, int], BoundaryEvidence]]:
    positions: Dict[str, Set[int]] = defaultdict(set)
    evidence: Dict[Tuple[str, int], BoundaryEvidence] = {}
    for pair in pairs:
        for ref in pair.boundary_refs().values():
            if ref is None:
                continue
            contig, pos = ref
            positions[contig].add(pos)
            evidence.setdefault((contig, pos), BoundaryEvidence())
    return {c: sorted(v) for c, v in positions.items()}, evidence


def process_bam(
    bam: str, samtools: str, anchor: int,
    positions: Dict[str, List[int]],
    evidence: Dict[Tuple[str, int], BoundaryEvidence],
    spanning_by_class: Dict[Tuple[str, int], Dict[str, Set[str]]],
    idn_event_rows: List[Dict[str, object]],
) -> None:
    """Scan HiFi alignments once, using primary + supplementary records only."""
    if not positions:
        return
    command = [samtools, "view", "-F", "260", bam] + sorted(positions)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for raw in process.stdout:
        f = raw.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        qname, rname, cigar = f[0], f[2], f[5]
        if rname not in positions or cigar == "*":
            continue
        try:
            flag, pos1 = int(f[1]), int(f[3])
        except ValueError:
            continue
        if flag & 0x4:
            continue
        aln_type = alignment_type(flag)
        if aln_type not in {"primary", "supplementary"}:
            continue
        ref_end, blocks, events = cigar_blocks_and_events(pos1, cigar)
        if ref_end < pos1 or not blocks:
            continue
        plist = positions[rname]
        lo = bisect.bisect_left(plist, max(2, pos1 - 1))
        hi = bisect.bisect_right(plist, ref_end + 1)
        for boundary in plist[lo:hi]:
            ev = evidence[(rname, boundary)]
            for event in events:
                if not event_hits_boundary_window(event, boundary, anchor):
                    continue
                if event.op == "I":
                    ev.i_qnames.add(qname); ev.i_lengths.append(event.length)
                elif event.op == "D":
                    ev.d_qnames.add(qname); ev.d_lengths.append(event.length)
                elif event.op == "N":
                    ev.n_qnames.add(qname); ev.n_lengths.append(event.length)
                idn_event_rows.append({
                    "assembler": "", "contig": rname,
                    "boundary_position_1based": boundary, "QNAME": qname,
                    "alignment_type": aln_type, "CIGAR": cigar,
                    "event_type": event.op, "event_length_bp": event.length,
                    "event_ref_start_1based": event.ref_start,
                    "event_ref_end_1based": event.ref_end,
                })
            left_ok = side_has_anchor(blocks, boundary, anchor, "left")
            right_ok = side_has_anchor(blocks, boundary, anchor, "right")
            span_ok = boundary_is_spanned(blocks, boundary, anchor)
            if left_ok:
                ev.left_reaching.add(qname)
            if right_ok:
                ev.right_reaching.add(qname)
            if span_ok:
                ev.spanning_reads.add(qname)
                spanning_by_class[(rname, boundary)][aln_type].add(qname)
    process.stdout.close()
    rc = process.wait()
    if rc != 0:
        fail(f"samtools view failed with exit code {rc}")


def boundary_status(ev: BoundaryEvidence, min_spanning: int, min_reaching: int) -> str:
    """Y: spanned by enough reads; N: reads reach both sides but none span; NA: not analyzable."""
    if len(ev.spanning_reads) >= min_spanning:
        return "Y"
    if ev.spanning_reads:
        return "NA"
    if len(ev.left_reaching) >= min_reaching and len(ev.right_reaching) >= min_reaching:
        return "N"
    return "NA"


def pair_boundary_statuses(
    pair: DuplicationPair,
    evidence: Dict[Tuple[str, int], BoundaryEvidence],
    min_spanning: int,
    min_reaching: int,
) -> Dict[str, str]:
    statuses: Dict[str, str] = {}
    for name, ref in pair.boundary_refs().items():
        if ref is None:
            statuses[name] = "NA"
        else:
            statuses[name] = boundary_status(evidence.get(ref, BoundaryEvidence()), min_spanning, min_reaching)
    return statuses


def classify_pair(pair: DuplicationPair, statuses: Dict[str, str]) -> Tuple[str, str, str]:
    """Return detailed state, assessment-complete Y/N, and Y/N/NA LR call."""
    if pair.duplication_type == "OVERLAPPING":
        return "NA_overlapping_copies", "N", "NA"
    if pair.duplication_type == "TANDEM":
        internal = statuses["copy1_right"]
        if internal == "Y":
            return "N_tandem_junction_spanned", "Y", "N"
        if pair.orientation == "INVERTED":
            return "NA_inverted_tandem", "N", "NA"
        outer_left = statuses["copy1_left"]
        outer_right = statuses["copy2_right"]
        required = (outer_left, internal, outer_right)
        if not all(x in {"Y", "N"} for x in required):
            return "NA_insufficient_reads", "N", "NA"
        if outer_left == "Y" and internal == "N" and outer_right == "Y":
            return "Y_false_duplication_pattern", "Y", "Y"
        return "N_no_false_duplication_pattern", "Y", "N"

    required = tuple(statuses.values())
    if not all(x in {"Y", "N"} for x in required):
        return "NA_insufficient_reads", "N", "NA"
    c1s = statuses["copy1_left"] == statuses["copy1_right"] == "Y"
    c2s = statuses["copy2_left"] == statuses["copy2_right"] == "Y"
    c1u = statuses["copy1_left"] == statuses["copy1_right"] == "N"
    c2u = statuses["copy2_left"] == statuses["copy2_right"] == "N"
    if (c1s and c2u) or (c2s and c1u):
        return "Y_false_duplication_pattern", "Y", "Y"
    return "N_no_false_duplication_pattern", "Y", "N"


def aggregate_lr_calls(calls: Sequence[str]) -> Tuple[str, str, str]:
    """Return assessment completeness, final-call analyzability, and Y/N/NA call."""
    if not calls:
        return "Y", "Y", "N"
    complete = "Y" if all(x in {"Y", "N"} for x in calls) else "N"
    if "Y" in calls:
        return complete, "Y", "Y"
    if "NA" in calls:
        return "N", "N", "NA"
    return "Y", "Y", "N"


def write_tsv(path: str, fields: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    os.replace(tmp, path)


def _detect_main() -> None:
    args = parse_args()
    for name in (
        "min_copy_bp", "hifi_anchor_bp",
        "hifi_min_spanning_reads", "hifi_min_reaching_reads",
    ):
        if getattr(args, name) < 1:
            fail(f"{name} must be positive")
    if not 0 < args.min_pident <= 100:
        fail("min_pident must be >0 and <=100")
    if not 0 < args.redundant_reciprocal_overlap <= 1:
        fail("redundant_reciprocal_overlap must be >0 and <=1")
    for path, label in (
        (args.candidate_paf, "candidate PAF"),
        (args.mapping_manifest, "mapping manifest"),
        (args.mag_map, "MAG map"),
    ):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"{label} missing or empty: {path}")

    mag_map = read_mag_map(args.mag_map)
    sources = read_mapping_manifest(args.mapping_manifest)
    routed = route_all_contigs(mag_map, sources)
    pairs, preprocessing_rows = read_duplication_pairs(
        args.candidate_paf, mag_map, routed,
        args.min_copy_bp, args.min_pident, args.redundant_reciprocal_overlap,
    )
    preprocessing_fields = [
        "preprocessing_action", "reason", "assembler", "MAG_ID", "PAIR_SCOPE",
        "DUPLICATION_TYPE", "ORIENTATION",
        "copy1_contig", "copy1_start", "copy1_end",
        "copy2_contig", "copy2_start", "copy2_end",
        "pident", "alignment_length", "representative_pair_id",
        "representative_copy1_contig", "representative_copy1_start", "representative_copy1_end",
        "representative_copy2_contig", "representative_copy2_start", "representative_copy2_end",
        "copy1_reciprocal_overlap", "copy2_reciprocal_overlap",
        "redundant_reciprocal_overlap_threshold",
    ]
    write_tsv(args.out_candidate_preprocessing, preprocessing_fields, preprocessing_rows)

    positions, hifi_evidence = all_boundary_positions(pairs)
    hifi_spanning_by_class: Dict[Tuple[str, int], Dict[str, Set[str]]] = defaultdict(
        lambda: {"primary": set(), "supplementary": set()}
    )
    idn_events: List[Dict[str, object]] = []

    for source in sources:
        source_positions = {c: p for c, p in positions.items() if routed[c].assembler == source.assembler}
        if not source_positions:
            continue
        print(f"HiFi boundary validation: {source.assembler}: {len(source_positions)} contigs")
        start_event = len(idn_events)
        process_bam(
            source.hifi_bam, args.samtools, args.hifi_anchor_bp,
            source_positions, hifi_evidence, hifi_spanning_by_class, idn_events,
        )
        for row in idn_events[start_event:]:
            row["assembler"] = source.assembler

    pair_rows: List[Dict[str, object]] = []
    boundary_rows: List[Dict[str, object]] = []
    for pair in pairs:
        statuses = pair_boundary_statuses(
            pair, hifi_evidence, args.hifi_min_spanning_reads, args.hifi_min_reaching_reads
        )
        state, complete, call = classify_pair(pair, statuses)

        for name, ref in pair.boundary_refs().items():
            if ref is None:
                contig = pair.copy1_contig if name.startswith("copy1") else pair.copy2_contig
                pos: object = "NA"
                ev = BoundaryEvidence()
                bstatus = "NA"
                class_spans = {"primary": set(), "supplementary": set()}
            else:
                contig, ipos = ref
                pos = ipos
                ev = hifi_evidence.get(ref, BoundaryEvidence())
                bstatus = statuses[name]
                class_spans = hifi_spanning_by_class.get(ref, {"primary": set(), "supplementary": set()})
            boundary_rows.append({
                "pair_id": pair.pair_id, "MAG_ID": pair.mag_id, "assembler": pair.assembler,
                "PAIR_SCOPE": pair.pair_scope, "DUPLICATION_TYPE": pair.duplication_type,
                "ORIENTATION": pair.orientation, "boundary": name, "contig": contig,
                "boundary_position_1based": pos,
                "HIFI_PRIMARY_SUPPLEMENTARY_SPANNING_UNIQUE_READS": len(ev.spanning_reads),
                "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_REACHING_UNIQUE_READS": len(ev.left_reaching),
                "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_REACHING_UNIQUE_READS": len(ev.right_reaching),
                "HIFI_PRIMARY_SUPPLEMENTARY_BOUNDARY_STATUS": bstatus,
                "HIFI_PRIMARY_SPANNING_UNIQUE_READS": len(class_spans["primary"]),
                "HIFI_SUPPLEMENTARY_SPANNING_UNIQUE_READS": len(class_spans["supplementary"]),
                "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_I_UNIQUE_READS": len(ev.i_qnames),
                "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_D_UNIQUE_READS": len(ev.d_qnames),
                "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_N_UNIQUE_READS": len(ev.n_qnames),
                "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_I_BP": max(ev.i_lengths, default=0),
                "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_D_BP": max(ev.d_lengths, default=0),
                "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_N_BP": max(ev.n_lengths, default=0),
                "LR_n_spanning_qnames": len(ev.spanning_reads),
                "LR_n_left_reaching_qnames": len(ev.left_reaching),
                "LR_n_right_reaching_qnames": len(ev.right_reaching),
                "LR_BOUNDARY_STATUS": bstatus,
                "LR_n_I_qnames": len(ev.i_qnames), "LR_n_D_qnames": len(ev.d_qnames), "LR_n_N_qnames": len(ev.n_qnames),
                "LR_max_I_bp": max(ev.i_lengths, default=0), "LR_max_D_bp": max(ev.d_lengths, default=0), "LR_max_N_bp": max(ev.n_lengths, default=0),
            })

        pair_rows.append({
            "pair_id": pair.pair_id, "MAG_ID": pair.mag_id, "assembler": pair.assembler,
            "PAIR_SCOPE": pair.pair_scope, "DUPLICATION_TYPE": pair.duplication_type, "ORIENTATION": pair.orientation,
            "TANDEM_ORIENTATION": pair.orientation if pair.duplication_type == "TANDEM" else "NA",
            "copy1_contig": pair.copy1_contig, "copy1_contig_length": pair.copy1_contig_length,
            "copy1_start": pair.copy1_start, "copy1_end": pair.copy1_end,
            "copy2_contig": pair.copy2_contig, "copy2_contig_length": pair.copy2_contig_length,
            "copy2_start": pair.copy2_start, "copy2_end": pair.copy2_end,
            "pident": f"{pair.pident:.4f}", "alignment_length": pair.alignment_length,
            "candidate_cluster_size": pair.candidate_cluster_size,
            "LR_COPY1_LEFT_SUPPORTED": statuses["copy1_left"], "LR_COPY1_RIGHT_SUPPORTED": statuses["copy1_right"],
            "LR_COPY2_LEFT_SUPPORTED": statuses["copy2_left"], "LR_COPY2_RIGHT_SUPPORTED": statuses["copy2_right"],
            "LR_ASSESSMENT": state, "LR_ASSESSMENT_COMPLETE_Y_N": complete,
            "LR_ANALYZABLE_Y_N": complete, "LR_FALSE_DUPLICATION": call,
            "FALSE_DUPLICATION": call,
        })

    pair_fields = [
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION", "TANDEM_ORIENTATION",
        "copy1_contig", "copy1_contig_length", "copy1_start", "copy1_end",
        "copy2_contig", "copy2_contig_length", "copy2_start", "copy2_end",
        "pident", "alignment_length", "candidate_cluster_size",
        "LR_COPY1_LEFT_SUPPORTED", "LR_COPY1_RIGHT_SUPPORTED", "LR_COPY2_LEFT_SUPPORTED", "LR_COPY2_RIGHT_SUPPORTED",
        "LR_ASSESSMENT", "LR_ASSESSMENT_COMPLETE_Y_N", "LR_ANALYZABLE_Y_N", "LR_FALSE_DUPLICATION",
        "FALSE_DUPLICATION",
    ]
    boundary_fields = [
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION",
        "boundary", "contig", "boundary_position_1based",
        "HIFI_PRIMARY_SUPPLEMENTARY_SPANNING_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_REACHING_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_REACHING_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_BOUNDARY_STATUS",
        "HIFI_PRIMARY_SPANNING_UNIQUE_READS", "HIFI_SUPPLEMENTARY_SPANNING_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_I_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_D_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_N_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_I_BP",
        "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_D_BP",
        "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_N_BP",
        "LR_n_spanning_qnames", "LR_n_left_reaching_qnames", "LR_n_right_reaching_qnames", "LR_BOUNDARY_STATUS",
        "LR_n_I_qnames", "LR_n_D_qnames", "LR_n_N_qnames", "LR_max_I_bp", "LR_max_D_bp", "LR_max_N_bp",
    ]
    write_tsv(args.out_all_pairs, pair_fields, pair_rows)
    write_tsv(args.out_all_boundaries, boundary_fields, boundary_rows)
    mag_pair_rows = [r for r in pair_rows if str(r["MAG_ID"]) not in {"", "NA"}]
    mag_pair_ids = {str(r["pair_id"]) for r in mag_pair_rows}
    mag_boundary_rows = [r for r in boundary_rows if str(r["pair_id"]) in mag_pair_ids]
    write_tsv(args.out_pairs, pair_fields, mag_pair_rows)
    write_tsv(args.out_boundaries, boundary_fields, mag_boundary_rows)

    idn_fields = [
        "assembler", "contig", "boundary_position_1based", "QNAME", "alignment_type", "CIGAR",
        "event_type", "event_length_bp", "event_ref_start_1based", "event_ref_end_1based",
    ]
    uniq_idn: List[Dict[str, object]] = []
    seen_idn: Set[Tuple[object, ...]] = set()
    for row in idn_events:
        key = tuple(row[x] for x in idn_fields)
        if key not in seen_idn:
            seen_idn.add(key)
            uniq_idn.append(row)
    write_tsv(args.out_idn_events, idn_fields, uniq_idn)
    idn_summary_rows: List[Dict[str, object]] = []
    for (contig, pos), ev in sorted(hifi_evidence.items()):
        idn_summary_rows.append({
            "assembler": routed[contig].assembler, "contig": contig, "boundary_position_1based": pos,
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_I_UNIQUE_READS": len(ev.i_qnames),
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_D_UNIQUE_READS": len(ev.d_qnames),
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_N_UNIQUE_READS": len(ev.n_qnames),
            "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_I_BP": max(ev.i_lengths, default=0),
            "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_D_BP": max(ev.d_lengths, default=0),
            "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_N_BP": max(ev.n_lengths, default=0),
            "n_I_qnames": len(ev.i_qnames), "n_D_qnames": len(ev.d_qnames), "n_N_qnames": len(ev.n_qnames),
            "max_I_bp": max(ev.i_lengths, default=0), "max_D_bp": max(ev.d_lengths, default=0), "max_N_bp": max(ev.n_lengths, default=0),
        })
    idn_summary_fields = [
        "assembler", "contig", "boundary_position_1based",
        "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_I_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_D_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_N_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_I_BP",
        "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_D_BP",
        "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_N_BP",
        "n_I_qnames", "n_D_qnames", "n_N_qnames", "max_I_bp", "max_D_bp", "max_N_bp",
    ]
    write_tsv(args.out_idn_summary, idn_summary_fields, idn_summary_rows)

    pairs_by_contig: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        pairs_by_contig[str(row["copy1_contig"])].append(row)
        if row["copy2_contig"] != row["copy1_contig"]:
            pairs_by_contig[str(row["copy2_contig"])].append(row)

    contig_rows: List[Dict[str, object]] = []
    for contig in sorted(routed):
        rows = pairs_by_contig.get(contig, [])
        calls = [str(r["LR_FALSE_DUPLICATION"]) for r in rows]
        complete, lr_analyzable, lr_call = aggregate_lr_calls(calls)
        contig_rows.append({
            "MAG_ID": mag_map.get(contig, "NA"), "assembler": routed[contig].assembler,
            "contig": contig, "length": routed[contig].lengths[contig],
            "n_duplicated_pairs": len(rows),
            "n_intra_contig_pairs": sum(r["PAIR_SCOPE"] == "INTRA_CONTIG" for r in rows),
            "n_inter_contig_pairs": sum(r["PAIR_SCOPE"] == "INTER_CONTIG" for r in rows),
            "n_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" for r in rows),
            "n_direct_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" and r["ORIENTATION"] == "DIRECT" for r in rows),
            "n_inverted_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" and r["ORIENTATION"] == "INVERTED" for r in rows),
            "n_interspersed_pairs": sum(r["DUPLICATION_TYPE"] == "INTERSPERSED" for r in rows),
            "n_direct_pairs": sum(r["ORIENTATION"] == "DIRECT" for r in rows),
            "n_inverted_pairs": sum(r["ORIENTATION"] == "INVERTED" for r in rows),
            "n_false_duplication_pairs": sum(r["LR_FALSE_DUPLICATION"] == "Y" for r in rows),
            "n_no_false_duplication_pairs": sum(r["LR_FALSE_DUPLICATION"] == "N" for r in rows),
            "n_not_analyzable_pairs": sum(r["LR_FALSE_DUPLICATION"] == "NA" for r in rows),
            "HIFI_PRIMARY_SUPPLEMENTARY_ASSESSMENT_COMPLETE_Y_N": complete,
            "HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N": lr_analyzable,
            "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_DUPLICATION_Y_N": lr_call,
            "LR_ASSESSMENT_COMPLETE_Y_N": complete, "LR_ANALYZABLE_Y_N": lr_analyzable,
            "LR_FALSE_DUPLICATION": lr_call, "FALSE_DUPLICATION": lr_call,
            "n_fully_evaluable_pairs": sum(r["LR_FALSE_DUPLICATION"] in {"Y", "N"} for r in rows),
            "n_false_duplication_candidates": sum(r["LR_FALSE_DUPLICATION"] == "Y" for r in rows),
            "n_no_false_duplication_signal_pairs": sum(r["LR_FALSE_DUPLICATION"] == "N" for r in rows),
            "n_insufficient_evidence_pairs": sum(r["LR_FALSE_DUPLICATION"] == "NA" for r in rows),
            "structural_assessment_status": (
                "Not_applicable_no_candidate_pairs" if not rows else
                "Complete" if all(r["LR_FALSE_DUPLICATION"] in {"Y", "N"} for r in rows) else
                "Partial" if any(r["LR_FALSE_DUPLICATION"] in {"Y", "N"} for r in rows) else
                "Insufficient_for_all_candidate_pairs"
            ),
            "structural_assessment_complete": "Y" if (not rows or all(r["LR_FALSE_DUPLICATION"] in {"Y", "N"} for r in rows)) else "N",
            "false_duplication_candidate_detected": "Y" if any(r["LR_FALSE_DUPLICATION"] == "Y" for r in rows) else "N",
            "false_duplication_assessment_call": (
                "Y" if any(r["LR_FALSE_DUPLICATION"] == "Y" for r in rows) else
                "U" if any(r["LR_FALSE_DUPLICATION"] == "NA" for r in rows) else "N"
            ),
            "false_duplication_call": (
                "False_duplication_candidate" if any(r["LR_FALSE_DUPLICATION"] == "Y" for r in rows) else
                "No_duplicated_structure_detected" if not rows else
                "No_false_duplication_candidate_detected_assessment_incomplete" if any(r["LR_FALSE_DUPLICATION"] == "NA" for r in rows) else
                "No_false_duplication_signal"
            ),
        })
    contig_fields = [
        "MAG_ID", "assembler", "contig", "length", "n_duplicated_pairs", "n_intra_contig_pairs", "n_inter_contig_pairs",
        "n_tandem_pairs", "n_direct_tandem_pairs", "n_inverted_tandem_pairs", "n_interspersed_pairs",
        "n_direct_pairs", "n_inverted_pairs",
        "LR_ASSESSMENT_COMPLETE_Y_N", "LR_ANALYZABLE_Y_N", "LR_FALSE_DUPLICATION", "FALSE_DUPLICATION",
        "n_fully_evaluable_pairs", "n_false_duplication_candidates", "n_no_false_duplication_signal_pairs", "n_insufficient_evidence_pairs",
        "structural_assessment_status", "structural_assessment_complete", "false_duplication_candidate_detected",
        "false_duplication_assessment_call", "false_duplication_call",
    ]
    write_tsv(args.out_all_contigs, contig_fields, contig_rows)
    mag_contig_rows = [r for r in contig_rows if str(r["MAG_ID"]) not in {"", "NA"}]
    write_tsv(args.out_contigs, contig_fields, mag_contig_rows)

    source_by_assembler = {x.assembler: x for x in sources}
    assembler_rows: List[Dict[str, object]] = []
    for assembler in sorted(source_by_assembler):
        pr = [r for r in pair_rows if r["assembler"] == assembler]
        cr = [r for r in contig_rows if r["assembler"] == assembler]
        assembler_rows.append({
            "assembler": assembler, "n_contigs": len(cr), "total_length_bp": sum(int(r["length"]) for r in cr),
            "n_contigs_with_duplicated_pairs": sum(int(r["n_duplicated_pairs"]) > 0 for r in cr),
            "n_duplicated_pairs": len(pr),
            "n_intra_contig_pairs": sum(r["PAIR_SCOPE"] == "INTRA_CONTIG" for r in pr),
            "n_inter_contig_pairs": sum(r["PAIR_SCOPE"] == "INTER_CONTIG" for r in pr),
            "n_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" for r in pr),
            "n_direct_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" and r["ORIENTATION"] == "DIRECT" for r in pr),
            "n_inverted_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" and r["ORIENTATION"] == "INVERTED" for r in pr),
            "n_interspersed_pairs": sum(r["DUPLICATION_TYPE"] == "INTERSPERSED" for r in pr),
            "n_false_duplication_pairs": sum(r["LR_FALSE_DUPLICATION"] == "Y" for r in pr),
            "n_no_false_duplication_pairs": sum(r["LR_FALSE_DUPLICATION"] == "N" for r in pr),
            "n_not_analyzable_pairs": sum(r["LR_FALSE_DUPLICATION"] == "NA" for r in pr),
            "n_false_duplication_contigs": sum(r["LR_FALSE_DUPLICATION"] == "Y" for r in cr),
        })
    assembler_fields = [
        "assembler", "n_contigs", "total_length_bp", "n_contigs_with_duplicated_pairs", "n_duplicated_pairs",
        "n_intra_contig_pairs", "n_inter_contig_pairs", "n_tandem_pairs", "n_direct_tandem_pairs",
        "n_inverted_tandem_pairs", "n_interspersed_pairs", "n_false_duplication_pairs",
        "n_no_false_duplication_pairs", "n_not_analyzable_pairs", "n_false_duplication_contigs",
    ]
    write_tsv(args.out_assemblers, assembler_fields, assembler_rows)

    mag_rows: List[Dict[str, object]] = []
    mags = sorted(set(mag_map.values()))
    for mag in mags:
        pr = [r for r in mag_pair_rows if r["MAG_ID"] == mag]
        calls = [str(r["LR_FALSE_DUPLICATION"]) for r in pr]
        complete, lr_analyzable, lr_call = aggregate_lr_calls(calls)
        mag_rows.append({
            "MAG_ID": mag, "n_contigs": sum(1 for c, m in mag_map.items() if m == mag),
            "n_duplicated_pairs": len(pr),
            "n_intra_contig_pairs": sum(r["PAIR_SCOPE"] == "INTRA_CONTIG" for r in pr),
            "n_inter_contig_pairs": sum(r["PAIR_SCOPE"] == "INTER_CONTIG" for r in pr),
            "n_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" for r in pr),
            "n_direct_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" and r["ORIENTATION"] == "DIRECT" for r in pr),
            "n_inverted_tandem_pairs": sum(r["DUPLICATION_TYPE"] == "TANDEM" and r["ORIENTATION"] == "INVERTED" for r in pr),
            "n_interspersed_pairs": sum(r["DUPLICATION_TYPE"] == "INTERSPERSED" for r in pr),
            "n_false_duplication_pairs": sum(r["LR_FALSE_DUPLICATION"] == "Y" for r in pr),
            "n_no_false_duplication_pairs": sum(r["LR_FALSE_DUPLICATION"] == "N" for r in pr),
            "n_not_analyzable_pairs": sum(r["LR_FALSE_DUPLICATION"] == "NA" for r in pr),
            "HIFI_PRIMARY_SUPPLEMENTARY_ASSESSMENT_COMPLETE_Y_N": complete,
            "HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N": lr_analyzable,
            "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_DUPLICATION_Y_N": lr_call,
            "LR_ASSESSMENT_COMPLETE_Y_N": complete, "LR_ANALYZABLE_Y_N": lr_analyzable,
            "LR_FALSE_DUPLICATION": lr_call, "FALSE_DUPLICATION": lr_call,
        })
    mag_fields = [
        "MAG_ID", "n_contigs", "n_duplicated_pairs", "n_intra_contig_pairs", "n_inter_contig_pairs",
        "n_tandem_pairs", "n_direct_tandem_pairs", "n_inverted_tandem_pairs", "n_interspersed_pairs",
        "n_false_duplication_pairs", "n_no_false_duplication_pairs", "n_not_analyzable_pairs",
        "LR_ASSESSMENT_COMPLETE_Y_N", "LR_ANALYZABLE_Y_N", "LR_FALSE_DUPLICATION", "FALSE_DUPLICATION",
    ]
    write_tsv(args.out_mags, mag_fields, mag_rows)

    print(f"All assembly contigs summarized: {len(contig_rows)}")
    print(f"Final-MAG contigs summarized: {len(mag_contig_rows)}")
    print(f"Duplication candidate pairs: {len(pair_rows)}")
    print(f"Intra-contig pairs: {sum(r['PAIR_SCOPE'] == 'INTRA_CONTIG' for r in pair_rows)}")
    print(f"Same-MAG inter-contig pairs: {sum(r['PAIR_SCOPE'] == 'INTER_CONTIG' for r in pair_rows)}")
    print(f"Tandem pairs: {sum(r['DUPLICATION_TYPE'] == 'TANDEM' for r in pair_rows)}")
    print(f"Direct tandem pairs: {sum(r['DUPLICATION_TYPE'] == 'TANDEM' and r['ORIENTATION'] == 'DIRECT' for r in pair_rows)}")
    print(f"Inverted tandem pairs: {sum(r['DUPLICATION_TYPE'] == 'TANDEM' and r['ORIENTATION'] == 'INVERTED' for r in pair_rows)}")
    print(f"HiFi false-duplication pairs: {sum(r['LR_FALSE_DUPLICATION'] == 'Y' for r in pair_rows)}")
    print(f"HiFi unresolved pairs: {sum(r['LR_FALSE_DUPLICATION'] == 'NA' for r in pair_rows)}")


def cmd_detect(argv):
    """Classify duplicated copies and validate them with HiFi reads"""
    sys.argv = [sys.argv[0] + " detect", *argv]
    return _detect_main()


def _annotate_parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--boundaries", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument(
        "--clipping",
        required=True,
        help="Validated clipping TSV, or NA when clipping was not run.",
    )
    parser.add_argument("--out-boundaries", required=True)
    parser.add_argument("--out-pairs", required=True)
    parser.add_argument("--out-contigs", required=True)
    parser.add_argument("--out-mags", required=True)
    return parser.parse_args()


def read_rows(path: str, required: Set[str], label: str) -> Tuple[List[str], List[Dict[str, str]]]:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        fail(f"{label} missing or empty: {path}")
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = required - set(fields)
        if missing:
            fail(f"{label} lacks columns: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    return fields, rows


def parse_positive_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        fail(f"invalid integer for {label}: {value!r}")
    if parsed < 1:
        fail(f"{label} must be positive: {parsed}")
    return parsed


@dataclass(frozen=True)
class ClippingSite:
    position: int
    coverage: int
    clipping_count: int
    clipping_ratio: float
    validated: bool


def read_contigs(path: str) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    _fields, rows = read_rows(
        path,
        {"MAG_ID", "assembler", "contig", "length"},
        "false-duplication contig table",
    )
    lengths: Dict[str, int] = {}
    for line_number, row in enumerate(rows, start=2):
        contig = row["contig"].strip()
        if not contig or contig in lengths:
            fail(f"blank or duplicate contig in contig table at line {line_number}: {contig!r}")
        lengths[contig] = parse_positive_int(row["length"], f"contig length at line {line_number}")
    if not rows:
        fail("false-duplication contig table contains no rows")
    return rows, lengths


def read_clipping(
    path: str,
    lengths: Mapping[str, int],
) -> Tuple[bool, Dict[str, Dict[int, ClippingSite]]]:
    if not path or path.upper() == "NA":
        return False, {}
    _fields, rows = read_rows(
        path,
        {
            "contig", "length", "pos_1based", "coverage",
            "clipping_count", "clipping_ratio", "Clipping",
        },
        "validated clipping table",
    )
    result: Dict[str, Dict[int, ClippingSite]] = defaultdict(dict)
    for line_number, row in enumerate(rows, start=2):
        contig = row["contig"].strip()
        if contig not in lengths:
            fail(f"clipping contig absent from false-duplication scope at line {line_number}: {contig}")
        length = parse_positive_int(row["length"], f"clipping contig length at line {line_number}")
        if length != lengths[contig]:
            fail(f"stale clipping contig length at line {line_number}: {contig}")
        position = parse_positive_int(row["pos_1based"], f"clipping position at line {line_number}")
        if position > length or position in result[contig]:
            fail(f"invalid or duplicate clipping position at line {line_number}: {contig}:{position}")
        coverage = parse_positive_int(row["coverage"], f"clipping coverage at line {line_number}")
        clipping_count = parse_positive_int(
            row["clipping_count"], f"clipping count at line {line_number}"
        )
        if clipping_count > coverage:
            fail(f"clipping count exceeds coverage at line {line_number}")
        try:
            ratio = float(row["clipping_ratio"])
        except (TypeError, ValueError) as exc:
            fail(f"invalid clipping ratio at line {line_number}: {row['clipping_ratio']!r}")
        expected = clipping_count / coverage
        if not 0 < ratio <= 1 or abs(ratio - expected) > 5.0e-7:
            fail(f"inconsistent clipping ratio at line {line_number}: {contig}:{position}")
        call = row["Clipping"].strip().upper()
        if call not in {"Y", "N", "NA"}:
            fail(f"invalid Clipping call at line {line_number}: {call!r}")
        result[contig][position] = ClippingSite(
            position=position,
            coverage=coverage,
            clipping_count=clipping_count,
            clipping_ratio=ratio,
            validated=(call == "Y"),
        )
    return True, result


def defining_boundaries(pair: Mapping[str, str]) -> Set[str]:
    """Return the LR-unsupported boundary names that establish a positive call."""
    if pair["LR_FALSE_DUPLICATION"] != "Y":
        return set()
    duplication_type = pair["DUPLICATION_TYPE"]
    if duplication_type == "TANDEM":
        if pair["LR_COPY1_RIGHT_SUPPORTED"] != "N":
            fail(f"positive tandem pair lacks an unsupported internal boundary: {pair['pair_id']}")
        return {"copy1_right"}

    copy1_unsupported = (
        pair["LR_COPY1_LEFT_SUPPORTED"] == "N"
        and pair["LR_COPY1_RIGHT_SUPPORTED"] == "N"
    )
    copy2_unsupported = (
        pair["LR_COPY2_LEFT_SUPPORTED"] == "N"
        and pair["LR_COPY2_RIGHT_SUPPORTED"] == "N"
    )
    copy1_supported = (
        pair["LR_COPY1_LEFT_SUPPORTED"] == "Y"
        and pair["LR_COPY1_RIGHT_SUPPORTED"] == "Y"
    )
    copy2_supported = (
        pair["LR_COPY2_LEFT_SUPPORTED"] == "Y"
        and pair["LR_COPY2_RIGHT_SUPPORTED"] == "Y"
    )
    if copy1_unsupported and copy2_supported:
        return {"copy1_left", "copy1_right"}
    if copy2_unsupported and copy1_supported:
        return {"copy2_left", "copy2_right"}
    fail(f"positive non-tandem pair lacks the expected supported/unsupported pattern: {pair['pair_id']}")


def matched_sites(
    contig: str,
    boundary_position: int,
    clipping: Mapping[str, Mapping[int, ClippingSite]],
) -> List[ClippingSite]:
    positions = {boundary_position - 1, boundary_position}
    return [
        clipping[contig][position]
        for position in sorted(positions)
        if position >= 1 and position in clipping.get(contig, {})
    ]


def format_positions(sites: Sequence[ClippingSite]) -> str:
    return ";".join(str(site.position) for site in sites) or "NA"


def format_offsets(sites: Sequence[ClippingSite], boundary_position: int) -> str:
    return ";".join(str(site.position - boundary_position) for site in sites) or "NA"


def max_or_na(values: Iterable[float | int], decimal: bool = False) -> str:
    materialized = list(values)
    if not materialized:
        return "NA"
    maximum = max(materialized)
    return f"{float(maximum):.6f}" if decimal else str(maximum)


def _annotate_write_tsv(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp, path)


def _annotate_main() -> None:
    args = _annotate_parse_args()
    pair_required = {
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION",
        "copy1_contig", "copy2_contig", "LR_COPY1_LEFT_SUPPORTED", "LR_COPY1_RIGHT_SUPPORTED",
        "LR_COPY2_LEFT_SUPPORTED", "LR_COPY2_RIGHT_SUPPORTED", "LR_FALSE_DUPLICATION",
    }
    _pair_fields, pairs = read_rows(args.pairs, pair_required, "false-duplication pair table")
    _boundary_fields, boundaries = read_rows(
        args.boundaries,
        {
            "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE",
            "ORIENTATION", "boundary", "contig", "boundary_position_1based",
        },
        "false-duplication boundary table",
    )
    contigs, lengths = read_contigs(args.contigs)
    clipping_available, clipping = read_clipping(args.clipping, lengths)

    pair_by_id: Dict[str, Dict[str, str]] = {}
    for pair in pairs:
        pair_id = pair["pair_id"]
        if not pair_id or pair_id in pair_by_id:
            fail(f"blank or duplicate pair_id: {pair_id!r}")
        if pair["LR_FALSE_DUPLICATION"] not in {"Y", "N", "NA"}:
            fail(f"invalid LR false-duplication call for {pair_id}")
        pair_by_id[pair_id] = pair

    boundary_by_pair: Dict[str, Dict[str, Dict[str, str]]] = defaultdict(dict)
    for row in boundaries:
        pair_id = row["pair_id"]
        boundary_name = row["boundary"]
        if pair_id not in pair_by_id:
            fail(f"boundary row references an unknown pair: {pair_id}")
        if boundary_name not in {"copy1_left", "copy1_right", "copy2_left", "copy2_right"}:
            fail(f"invalid boundary name for {pair_id}: {boundary_name}")
        if boundary_name in boundary_by_pair[pair_id]:
            fail(f"duplicate boundary row for {pair_id}/{boundary_name}")
        if row["contig"] not in lengths:
            fail(f"boundary contig absent from contig table: {row['contig']}")
        boundary_by_pair[pair_id][boundary_name] = row
    for pair_id in pair_by_id:
        if set(boundary_by_pair[pair_id]) != {
            "copy1_left", "copy1_right", "copy2_left", "copy2_right"
        }:
            fail(f"pair does not have exactly four boundary rows: {pair_id}")

    boundary_output: List[Dict[str, object]] = []
    boundary_annotations: Dict[Tuple[str, str], Dict[str, object]] = {}
    for pair in pairs:
        pair_id = pair["pair_id"]
        defining = defining_boundaries(pair)
        for boundary_name in ("copy1_left", "copy1_right", "copy2_left", "copy2_right"):
            row = boundary_by_pair[pair_id][boundary_name]
            position_text = row["boundary_position_1based"].strip()
            is_defining = boundary_name in defining
            if position_text.upper() == "NA":
                sites: List[ClippingSite] = []
                position: int | None = None
            else:
                position = parse_positive_int(
                    position_text, f"boundary position for {pair_id}/{boundary_name}"
                )
                if position > lengths[row["contig"]]:
                    fail(f"boundary position exceeds contig length: {pair_id}/{boundary_name}")
                sites = matched_sites(row["contig"], position, clipping) if clipping_available else []
            validated = [site for site in sites if site.validated]

            if not is_defining:
                interpretation = "NOT_A_FALSE_DUPLICATION_DEFINING_BOUNDARY"
            elif not clipping_available:
                interpretation = "CLIPPING_DATA_NOT_AVAILABLE"
            elif validated:
                interpretation = "FALSE_DUPLICATION_BOUNDARY_WITH_VALIDATED_CLIPPING"
            elif sites:
                interpretation = "CLIPPING_CANDIDATE_NOT_VALIDATED_AT_FALSE_DUPLICATION_BOUNDARY"
            else:
                interpretation = "NO_CLIPPING_CANDIDATE_AT_FALSE_DUPLICATION_BOUNDARY"

            annotated: Dict[str, object] = {
                "pair_id": pair_id,
                "MAG_ID": row["MAG_ID"],
                "assembler": row["assembler"],
                "PAIR_SCOPE": row["PAIR_SCOPE"],
                "DUPLICATION_TYPE": row["DUPLICATION_TYPE"],
                "ORIENTATION": row["ORIENTATION"],
                "boundary": boundary_name,
                "contig": row["contig"],
                "boundary_position_1based": position_text,
                "FALSE_DUPLICATION_PAIR_Y_N_NA": pair["LR_FALSE_DUPLICATION"],
                "FALSE_DUPLICATION_DEFINING_BOUNDARY_Y_N": "Y" if is_defining else "N",
                "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
                "CLIPPING_CANDIDATE_AT_BOUNDARY_Y_N_NA": (
                    "Y" if sites else "N" if clipping_available and position is not None else "NA"
                ),
                "VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA": (
                    "Y" if validated else "N" if clipping_available and position is not None else "NA"
                ),
                "matched_clipping_positions_1based": format_positions(sites),
                "matched_clipping_offsets_from_boundary_bp": (
                    format_offsets(sites, position) if position is not None else "NA"
                ),
                "matched_validated_clipping_positions_1based": format_positions(validated),
                "max_matched_clipping_coverage": max_or_na(site.coverage for site in sites),
                "max_matched_clipping_count": max_or_na(site.clipping_count for site in sites),
                "max_matched_clipping_ratio": max_or_na(
                    (site.clipping_ratio for site in sites), decimal=True
                ),
                "max_validated_clipping_ratio": max_or_na(
                    (site.clipping_ratio for site in validated), decimal=True
                ),
                "FALSE_DUPLICATION_CLIPPING_INTERPRETATION": interpretation,
            }
            boundary_output.append(annotated)
            boundary_annotations[(pair_id, boundary_name)] = annotated

    boundary_output_fields = [
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION",
        "boundary", "contig", "boundary_position_1based", "FALSE_DUPLICATION_PAIR_Y_N_NA",
        "FALSE_DUPLICATION_DEFINING_BOUNDARY_Y_N", "CLIPPING_DATA_AVAILABLE_Y_N",
        "CLIPPING_CANDIDATE_AT_BOUNDARY_Y_N_NA", "VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA",
        "matched_clipping_positions_1based", "matched_clipping_offsets_from_boundary_bp",
        "matched_validated_clipping_positions_1based", "max_matched_clipping_coverage",
        "max_matched_clipping_count", "max_matched_clipping_ratio",
        "max_validated_clipping_ratio", "FALSE_DUPLICATION_CLIPPING_INTERPRETATION",
    ]
    _annotate_write_tsv(args.out_boundaries, boundary_output_fields, boundary_output)

    pair_output: List[Dict[str, object]] = []
    for pair in pairs:
        pair_id = pair["pair_id"]
        defining = defining_boundaries(pair)
        defining_rows = [boundary_annotations[(pair_id, name)] for name in sorted(defining)]
        candidate_count = sum(
            row["CLIPPING_CANDIDATE_AT_BOUNDARY_Y_N_NA"] == "Y" for row in defining_rows
        )
        validated_count = sum(
            row["VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA"] == "Y" for row in defining_rows
        )
        positive = pair["LR_FALSE_DUPLICATION"] == "Y"
        if not positive:
            any_call = all_call = "NA"
            interpretation = "NOT_APPLICABLE_FALSE_DUPLICATION_NOT_DETECTED"
        elif not clipping_available:
            any_call = all_call = "NA"
            interpretation = "CLIPPING_DATA_NOT_AVAILABLE"
        elif validated_count == len(defining):
            any_call = all_call = "Y"
            interpretation = "ALL_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING"
        elif validated_count > 0:
            any_call, all_call = "Y", "N"
            interpretation = "SOME_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING"
        else:
            any_call = all_call = "N"
            interpretation = "NO_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARIES"
        pair_output.append({
            "pair_id": pair_id,
            "MAG_ID": pair["MAG_ID"],
            "assembler": pair["assembler"],
            "PAIR_SCOPE": pair["PAIR_SCOPE"],
            "DUPLICATION_TYPE": pair["DUPLICATION_TYPE"],
            "ORIENTATION": pair["ORIENTATION"],
            "FALSE_DUPLICATION_PAIR_Y_N_NA": pair["LR_FALSE_DUPLICATION"],
            "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
            "n_false_duplication_defining_boundaries": len(defining),
            "n_defining_boundaries_with_clipping_candidate": candidate_count,
            "n_defining_boundaries_with_validated_clipping": validated_count,
            "ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA": any_call,
            "ALL_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING_Y_N_NA": all_call,
            "FALSE_DUPLICATION_CLIPPING_INTERPRETATION": interpretation,
        })

    pair_output_fields = [
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION",
        "FALSE_DUPLICATION_PAIR_Y_N_NA", "CLIPPING_DATA_AVAILABLE_Y_N",
        "n_false_duplication_defining_boundaries",
        "n_defining_boundaries_with_clipping_candidate",
        "n_defining_boundaries_with_validated_clipping",
        "ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA",
        "ALL_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING_Y_N_NA",
        "FALSE_DUPLICATION_CLIPPING_INTERPRETATION",
    ]
    _annotate_write_tsv(args.out_pairs, pair_output_fields, pair_output)

    defining_rows_by_contig: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in boundary_output:
        if (
            row["FALSE_DUPLICATION_PAIR_Y_N_NA"] == "Y"
            and row["FALSE_DUPLICATION_DEFINING_BOUNDARY_Y_N"] == "Y"
        ):
            defining_rows_by_contig[str(row["contig"])].append(row)

    summary_fields = [
        "MAG_ID", "assembler", "contig", "length", "CLIPPING_DATA_AVAILABLE_Y_N",
        "n_false_duplication_pairs", "n_false_duplication_pairs_with_any_validated_clipping",
        "n_false_duplication_defining_boundaries",
        "n_false_duplication_defining_boundaries_with_validated_clipping",
        "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA",
        "FALSE_DUPLICATION_CLIPPING_SUMMARY",
    ]
    contig_output: List[Dict[str, object]] = []
    for contig_row in contigs:
        defining_rows = defining_rows_by_contig.get(contig_row["contig"], [])
        positive_pair_ids = {str(row["pair_id"]) for row in defining_rows}
        pairs_with_clipping = {
            str(row["pair_id"])
            for row in defining_rows
            if row["VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA"] == "Y"
        }
        pair_with_clipping = len(pairs_with_clipping)
        if not clipping_available:
            any_call = "NA"
            summary = "CLIPPING_DATA_NOT_AVAILABLE"
        elif pair_with_clipping:
            any_call = "Y"
            summary = "FALSE_DUPLICATION_WITH_VALIDATED_CLIPPING_CO_OCCURRENCE"
        elif positive_pair_ids:
            any_call = "N"
            summary = "FALSE_DUPLICATION_WITHOUT_VALIDATED_CLIPPING_CO_OCCURRENCE"
        else:
            any_call = "N"
            summary = "NO_FALSE_DUPLICATION_PAIR"
        contig_output.append({
            "MAG_ID": contig_row["MAG_ID"],
            "assembler": contig_row["assembler"],
            "contig": contig_row["contig"],
            "length": contig_row["length"],
            "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
            "n_false_duplication_pairs": len(positive_pair_ids),
            "n_false_duplication_pairs_with_any_validated_clipping": pair_with_clipping,
            "n_false_duplication_defining_boundaries": len(defining_rows),
            "n_false_duplication_defining_boundaries_with_validated_clipping": sum(
                row["VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA"] == "Y" for row in defining_rows
            ),
            "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA": any_call,
            "FALSE_DUPLICATION_CLIPPING_SUMMARY": summary,
        })
    _annotate_write_tsv(args.out_contigs, summary_fields, contig_output)

    mag_pair_rows: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in pair_output:
        mag = str(row["MAG_ID"])
        if mag not in {"", "NA"} and row["FALSE_DUPLICATION_PAIR_Y_N_NA"] == "Y":
            mag_pair_rows[mag].append(row)
    mag_contig_count: Dict[str, int] = defaultdict(int)
    for row in contigs:
        mag = row["MAG_ID"]
        if mag not in {"", "NA"}:
            mag_contig_count[mag] += 1

    mag_fields = [
        "MAG_ID", "n_contigs", "CLIPPING_DATA_AVAILABLE_Y_N", "n_false_duplication_pairs",
        "n_false_duplication_pairs_with_any_validated_clipping",
        "n_false_duplication_defining_boundaries",
        "n_false_duplication_defining_boundaries_with_validated_clipping",
        "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA",
        "FALSE_DUPLICATION_CLIPPING_SUMMARY",
    ]
    mag_output: List[Dict[str, object]] = []
    for mag in sorted(mag_contig_count):
        positive_rows = mag_pair_rows.get(mag, [])
        pair_with_clipping = sum(
            row["ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA"] == "Y"
            for row in positive_rows
        )
        if not clipping_available:
            any_call = "NA"
            summary = "CLIPPING_DATA_NOT_AVAILABLE"
        elif pair_with_clipping:
            any_call = "Y"
            summary = "FALSE_DUPLICATION_WITH_VALIDATED_CLIPPING_CO_OCCURRENCE"
        elif positive_rows:
            any_call = "N"
            summary = "FALSE_DUPLICATION_WITHOUT_VALIDATED_CLIPPING_CO_OCCURRENCE"
        else:
            any_call = "N"
            summary = "NO_FALSE_DUPLICATION_PAIR"
        mag_output.append({
            "MAG_ID": mag,
            "n_contigs": mag_contig_count[mag],
            "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
            "n_false_duplication_pairs": len(positive_rows),
            "n_false_duplication_pairs_with_any_validated_clipping": pair_with_clipping,
            "n_false_duplication_defining_boundaries": sum(
                int(row["n_false_duplication_defining_boundaries"]) for row in positive_rows
            ),
            "n_false_duplication_defining_boundaries_with_validated_clipping": sum(
                int(row["n_defining_boundaries_with_validated_clipping"]) for row in positive_rows
            ),
            "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA": any_call,
            "FALSE_DUPLICATION_CLIPPING_SUMMARY": summary,
        })
    _annotate_write_tsv(args.out_mags, mag_fields, mag_output)

    print(f"Clipping data available: {'Y' if clipping_available else 'N'}")
    print(f"False-duplication pairs annotated: {len(pair_output)}")
    print(
        "Positive false-duplication pairs with any validated clipping: "
        f"{sum(row['ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA'] == 'Y' for row in pair_output)}"
    )


def cmd_annotate(argv):
    """Report validated clipping at the duplication boundaries"""
    sys.argv = [sys.argv[0] + " annotate", *argv]
    return _annotate_main()


COMMANDS = {
    'detect': cmd_detect,
    'annotate': cmd_annotate,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
