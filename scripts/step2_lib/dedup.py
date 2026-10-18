"""Dedup pair results to unique candidate loci.

If multiple queries hit the same genomic region (within `merge_bp`), collapse
them to ONE locus and preserve all query-to-locus attributions.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .pair import RepeatPair


@dataclass
class Locus:
    contig: str
    start0: int
    end0: int
    strand: str
    n_pairs: int = 0
    queries: set[str] = field(default_factory=set)
    attributions: list[dict] = field(default_factory=list)
    context: str = ""   # filled by context module


def dedup(pairs: list[RepeatPair], merge_bp: int) -> list[Locus]:
    by_contig: dict[tuple[str, str], list[RepeatPair]] = defaultdict(list)
    for p in pairs:
        # pair footprint = [anchor_start0, second_end0)
        start = min(p.anchor_start0, p.second_start0)
        end = max(p.anchor_end0, p.second_end0)
        by_contig[(p.contig, p.strand)].append(p)

    loci: list[Locus] = []
    for (contig, strand), lst in by_contig.items():
        # sort by footprint start
        lst.sort(key=lambda p: min(p.anchor_start0, p.second_start0))
        cur: Locus | None = None
        for p in lst:
            start = min(p.anchor_start0, p.second_start0)
            end = max(p.anchor_end0, p.second_end0)
            if cur is None or start - cur.end0 > merge_bp:
                cur = Locus(contig=contig, start0=start, end0=end, strand=strand)
                loci.append(cur)
            cur.start0 = min(cur.start0, start)
            cur.end0 = max(cur.end0, end)
            cur.n_pairs += 1
            cur.queries.add(p.anchor_qseqid)
            cur.queries.add(p.second_qseqid)
            cur.attributions.append({
                "anchor_qseqid": p.anchor_qseqid,
                "second_qseqid": p.second_qseqid,
                "anchor_pident": p.anchor_pident,
                "second_pident": p.second_pident,
                "spacer_len": p.spacer_len,
                "pair_kind": p.pair_kind,
                "anchor_start0": p.anchor_start0,
                "anchor_end0": p.anchor_end0,
                "second_start0": p.second_start0,
                "second_end0": p.second_end0,
            })
    return loci


LOCUS_COLS = ["locus_id", "contig", "strand", "start0", "end0",
               "length", "n_pairs", "n_unique_queries", "query_set"]


def write_loci_tsv(loci: list[Locus], out: Path, strict_queries: set[str]) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(LOCUS_COLS + ["strict_support", "context", "pair_kinds", "n_mini_array_pairs", "n_half_repeat_pairs"])
        for idx, L in enumerate(sorted(loci, key=lambda x: (x.contig, x.start0, x.strand))):
            lid = f"loc_{idx+1:06d}"
            strict = any(q in strict_queries for q in L.queries)
            w.writerow([
                lid, L.contig, L.strand, L.start0, L.end0,
                L.end0 - L.start0, L.n_pairs, len(L.queries),
                ",".join(sorted(L.queries)),
                "true" if strict else "false",
                L.context,
                ",".join(sorted({a["pair_kind"] for a in L.attributions})),
                sum(a["pair_kind"] == "mini_array" for a in L.attributions),
                sum(a["pair_kind"] == "half_repeat" for a in L.attributions),
            ])


ATTR_COLS = ["locus_id", "contig", "strand", "anchor_qseqid",
              "second_qseqid", "anchor_pident", "second_pident",
              "spacer_len", "pair_kind", "anchor_start0", "anchor_end0", "second_start0", "second_end0"]


def write_locus_attributions(loci: list[Locus], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ATTR_COLS)
        for idx, L in enumerate(sorted(loci, key=lambda x: (x.contig, x.start0, x.strand))):
            lid = f"loc_{idx+1:06d}"
            for a in L.attributions:
                w.writerow([
                    lid, L.contig, L.strand,
                    a["anchor_qseqid"], a["second_qseqid"],
                    a["anchor_pident"], a["second_pident"],
                    a["spacer_len"],
                    a["pair_kind"], a["anchor_start0"], a["anchor_end0"],
                    a["second_start0"], a["second_end0"],
                ])
