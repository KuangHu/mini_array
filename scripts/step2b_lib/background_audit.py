"""Re-bucket background hits and audit the rep1 anomaly.

Three audits:
  1. Count real vs each bg replicate using the SAME pipeline and the SAME
     bucket semantics (so pair-count and locus-count are not mixed).
  2. Trace rep1's 113 pairs → query → contig → locus.
  3. Check shuffle validity (dinuc composition preserved) and similarity to
     the original query / its reverse complement.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import edlib


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    rows = []
    cur_h, cur_s = None, []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur_h is not None:
                    rows.append((cur_h, "".join(cur_s)))
                cur_h = ln[1:].strip().split()[0]
                cur_s = []
            else:
                cur_s.append(ln.strip().upper())
        if cur_h is not None:
            rows.append((cur_h, "".join(cur_s)))
    return rows


def dinuc_counts(seq: str) -> Counter:
    return Counter(seq[i:i+2] for i in range(len(seq) - 1))


def mono_counts(seq: str) -> Counter:
    return Counter(seq)


def audit_shuffle(real_fa: Path, bg_fa: Path) -> list[dict]:
    """For each (real, bg) query pair with matching qseqid, report whether
    the shuffle preserves mono and di composition, and how similar bg is to
    the original (or its RC)."""
    real = {h.split("|")[0]: s for h, s in read_fasta(real_fa)}
    bg = {h.split("|")[0]: s for h, s in read_fasta(bg_fa)}
    rows = []
    for qid in sorted(real):
        rs, bs = real[qid], bg.get(qid)
        if bs is None:
            continue
        mono_ok = mono_counts(rs) == mono_counts(bs)
        dinuc_ok = dinuc_counts(rs) == dinuc_counts(bs)
        d_fwd = edlib.align(rs, bs, task="distance", mode="NW")["editDistance"]
        d_rc = edlib.align(rs, revcomp(bs), task="distance", mode="NW")["editDistance"]
        min_edit = min(d_fwd, d_rc)
        sim = 1 - min_edit / max(len(rs), len(bs))
        rows.append({
            "qseqid": qid,
            "mono_preserved": mono_ok,
            "dinuc_preserved": dinuc_ok,
            "edit_distance_vs_original": d_fwd,
            "edit_distance_vs_rc": d_rc,
            "min_edit_distance": min_edit,
            "similarity_to_closest": round(sim, 3),
            "orig_len": len(rs),
            "bg_len": len(bs),
        })
    return rows


def rebucket_pairs(pairs_tsv: Path,
                   contig_to_genome: dict[str, str],
                   context_fn) -> dict:
    """Walk a pairs TSV, classify each pair by context_fn, return counters."""
    if not pairs_tsv.is_file():
        return {"n_pairs": 0, "bucket_counts": {}}
    buckets = Counter()
    uniq_loci = {}
    pair_rows = []
    with pairs_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                r["anchor_start0"] = int(r["anchor_start0"])
                r["second_end0"] = int(r["second_end0"])
            except ValueError:
                continue
            pair_rows.append(r)
            bucket = context_fn(r["contig"],
                                 min(int(r["anchor_start0"]),
                                     int(r["second_start0"])),
                                 max(int(r["anchor_end0"]),
                                     int(r["second_end0"])))
            buckets[bucket] += 1
            key = (r["contig"], r["strand"],
                   int(r["anchor_start0"]) // 60,     # cheap dedup proxy
                   int(r["second_end0"]) // 60)
            uniq_loci[key] = bucket
    uniq_counts = Counter(uniq_loci.values())
    return {
        "n_pairs": len(pair_rows),
        "pair_bucket_counts": dict(buckets),
        "n_unique_loci_proxy": len(uniq_loci),
        "locus_bucket_counts": dict(uniq_counts),
    }
