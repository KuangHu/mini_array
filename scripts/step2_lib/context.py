"""Classify each pair's locus by its spatial relationship to known arrays.

Buckets (user requested in Step-2 instructions):
  - `known_array_interior`      — pair's genomic footprint overlaps a
                                    bona-fide + coord-pass array by ≥50%
  - `suspected_terminal_extension` — abuts a known array (overlap <50% but
                                      touches within 50 bp)
  - `possible_soft_context`     — overlaps a non-bona-fide-coord-pass but
                                    coord-valid candidate (Alternative,
                                    Possible, Possible Discarded, Low score)
  - `contig_edge`               — pair footprint lies within 100 bp of a
                                    contig end
  - `standalone`                — none of the above
"""
from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class KnownRegion:
    contig: str
    start0: int
    end0: int
    source: str      # "known_array_bed" | "soft_context"


def load_known_array_bed(bed: Path) -> dict[str, list[tuple[int, int]]]:
    """Return {contig: [(start0,end0), …]} from a BED6 file."""
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with bed.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            out[p[0]].append((int(p[1]), int(p[2])))
    for v in out.values():
        v.sort()
    return out


def load_soft_context(arrays_all_tsv: Path) -> dict[str, list[tuple[int, int]]]:
    """Load candidates that are coord-valid but not bona-fide — these form the
    soft-context mask (Alternative/Possible/Possible Discarded/Low score)."""
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with arrays_all_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["coordinate_check_pass"] != "true":
                continue
            if r["category"] == "bona_fide":
                continue
            out[r["contig_uid"]].append((int(r["start0"]), int(r["end0"])))
    for v in out.values():
        v.sort()
    return out


def contig_lengths_from_staged(staged_dir: Path) -> dict[str, int]:
    """Walk every staged FASTA, return {contig_uid: length}."""
    lengths: dict[str, int] = {}
    for f in sorted(staged_dir.glob("*.fna")):
        cur = None
        total = 0
        with f.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur is not None:
                        lengths[cur] = total
                    cur = line[1:].split()[0]
                    total = 0
                else:
                    total += len(line.strip())
            if cur is not None:
                lengths[cur] = total
    return lengths


def classify_locus(contig: str, start0: int, end0: int,
                   hard: dict[str, list[tuple[int, int]]],
                   soft: dict[str, list[tuple[int, int]]],
                   contig_lengths: dict[str, int],
                   edge_bp: int = 100) -> str:
    pair_len = max(1, end0 - start0)
    best_overlap_hard = 0
    for h0, h1 in hard.get(contig, []):
        ov = max(0, min(end0, h1) - max(start0, h0))
        if ov > best_overlap_hard:
            best_overlap_hard = ov
    if best_overlap_hard / pair_len >= 0.50:
        return "known_array_interior"
    # terminal extension: touches a known array within 50 bp but not interior
    for h0, h1 in hard.get(contig, []):
        if (h0 - 50) <= start0 <= h1 + 50 or (h0 - 50) <= end0 <= h1 + 50:
            return "suspected_terminal_extension"
    # soft context (non-bona-fide valid candidates)
    for h0, h1 in soft.get(contig, []):
        if max(0, min(end0, h1) - max(start0, h0)) > 0:
            return "possible_soft_context"
    # contig edge
    L = contig_lengths.get(contig)
    if L is not None and (start0 < edge_bp or end0 > L - edge_bp):
        return "contig_edge"
    return "standalone"
