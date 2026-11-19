"""Ribosomal-protein family counts per contig from the anvi'o SCG table.

The SCG table is written by modules/lib_scg.sh with the columns
contig, hmm_name, gene_callers_id and source. hmm_name may carry the HMM
collection as a prefix (``Bacteria_71|Ribosomal_L2``).
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from typing import Dict, Set, Tuple

RIBOSOMAL_PATTERN = re.compile(r"^ribosom", re.IGNORECASE)
SCG_SOURCES = ("Bacteria_71", "Archaea_76")


def _source_and_family(row: Dict[str, str]) -> Tuple[str, str]:
    marker = (row.get("hmm_name") or row.get("marker") or "").strip()
    source = (row.get("source") or row.get("hmm_source") or "").strip()
    if "|" in marker:
        prefix, rest = marker.split("|", 1)
        if prefix in SCG_SOURCES:
            source, marker = prefix, rest
    return source, marker


def load_ribosomal_families(path: str) -> Dict[str, Tuple[str, ...]]:
    """Return contig -> sorted distinct ribosomal-protein HMM families.

    Bacteria_71 and Archaea_76 both contain ribosomal-protein models, so the
    same protein can be hit twice under different names. Families are counted
    within each collection and the collection with the most families is kept.
    """
    per_source: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if "contig" not in fields or not ({"hmm_name", "marker"} & fields):
            raise SystemExit(f"ERROR: SCG table needs contig and hmm_name columns: {path}")
        for row in reader:
            contig = (row.get("contig") or "").strip()
            source, family = _source_and_family(row)
            if contig and family and RIBOSOMAL_PATTERN.match(family):
                per_source[contig][source or "NA"].add(family)

    families: Dict[str, Tuple[str, ...]] = {}
    for contig, by_source in per_source.items():
        best = max(by_source.items(), key=lambda item: (len(item[1]), item[0]))[1]
        families[contig] = tuple(sorted(best))
    return families
