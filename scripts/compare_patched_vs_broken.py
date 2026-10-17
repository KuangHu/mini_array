#!/usr/bin/env python3
"""Pre/post comparison report for the patched identify re-run.

Reads the preserved broken-feature-extraction outputs (results/
step1_broken_feature_extraction_min_score_0.90/arrays.all.tsv) and the newly
regenerated patched outputs (results/step1_patched_min_score_0.90/arrays.all.tsv)
and emits a diff report:

  - category count deltas
  - QC tier count deltas
  - score distribution shifts
  - per-array promotions (locus-matched)

Writes `results/step1_patched_min_score_0.90/pre_post_comparison.md`.

Note: the broken tree's arrays.all.tsv has NEW qc_tier columns because it was
overwritten once during Phase E (see its BROKEN.md). We ignore qc_tier when
comparing — only categories and scores are reliable evidence from that tree.
"""
from __future__ import annotations

import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
BROKEN = PROJECT / "results" / "step1_broken_feature_extraction_min_score_0.90"
PATCHED = PROJECT / "results" / "step1_patched_min_score_0.90"


def _read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _bucket_by_locus(rows: list[dict]) -> dict[tuple, list[dict]]:
    """Index rows by (genome_id, contig_uid, start0, end0, category) for pairing."""
    out: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        key = (r["genome_id"], r["contig_uid"],
               int(r["start0"]), int(r["end0"]))
        out[key].append(r)
    return out


def main() -> int:
    import sys as _s
    _s.path.insert(0, "/home/kuangh/tools/CRISPR_array_related/scripts")
    from _slurm_gate import require_slurm
    require_slurm()
    broken = _read(BROKEN / "arrays.all.tsv")
    patched = _read(PATCHED / "arrays.all.tsv")

    out = []
    out.append("# Patched vs broken-feature-extraction comparison")
    out.append("")
    out.append(f"- broken: `{BROKEN.name}` — {len(broken)} candidates")
    out.append(f"- patched: `{PATCHED.name}` — {len(patched)} candidates")
    out.append("")

    out.append("## Category distribution")
    out.append("")
    out.append("| category | broken | patched | Δ |")
    out.append("|---|---:|---:|---:|")
    bc = Counter(r["category"] for r in broken)
    pc = Counter(r["category"] for r in patched)
    for cat in sorted(set(bc) | set(pc)):
        b, p = bc.get(cat, 0), pc.get(cat, 0)
        out.append(f"| `{cat}` | {b} | {p} | {p-b:+d} |")
    out.append("")

    out.append("## QC tier distribution (patched tree only — broken tree's QC was overwritten)")
    out.append("")
    out.append("| qc_tier | patched |")
    out.append("|---|---:|")
    for tier, n in sorted(Counter(r["qc_tier"] for r in patched).items(),
                          key=lambda kv: -kv[1]):
        out.append(f"| `{tier}` | {n} |")
    out.append("")

    # Score deltas on locus-matched pairs. We match on (genome, contig, start0, end0)
    # — the patched run may re-call boundaries slightly, so this is a conservative
    # inner join; mismatches are reported separately.
    bi = _bucket_by_locus(broken)
    pi = _bucket_by_locus(patched)
    matched_pairs: list[tuple[dict, dict]] = []
    for key, prows in pi.items():
        for pr in prows:
            for br in bi.get(key, []):
                if br["category_raw"] == pr["category_raw"]:
                    matched_pairs.append((br, pr))
                    break
            else:
                # try: same locus, any category
                if bi.get(key):
                    matched_pairs.append((bi[key][0], pr))
    out.append(f"## Locus-matched score shifts")
    out.append("")
    out.append(f"- exact locus+category matches: **{len(matched_pairs)}**")
    deltas = [float(p["score"]) - float(b["score"]) for b, p in matched_pairs]
    if deltas:
        out.append(f"- score Δ: mean=**{statistics.mean(deltas):+.3f}** "
                   f"median={statistics.median(deltas):+.3f} "
                   f"min={min(deltas):+.3f} max={max(deltas):+.3f}")
    out.append("")

    # Promotions: candidates that were non-bona-fide in broken and are now bona-fide
    out.append("## Promotions to Bona-fide")
    out.append("")
    promoted = [(b, p) for b, p in matched_pairs
                if b["category"] != "bona_fide" and p["category"] == "bona_fide"]
    out.append(f"**{len(promoted)}** candidates promoted from non-bona-fide → bona-fide.")
    out.append("")
    if promoted:
        out.append("| genome | contig | start | pre_cat | pre_score | post_score | n_rep |")
        out.append("|---|---|---:|---|---:|---:|---:|")
        for b, p in sorted(promoted, key=lambda bp: -float(bp[1]["score"]))[:30]:
            out.append(f"| {b['genome_id']} | {b['contig_uid'][:14]} | "
                       f"{b['start0']} | {b['category']} | {float(b['score']):.3f} | "
                       f"{float(p['score']):.3f} | {p['n_repeats']} |")
        if len(promoted) > 30:
            out.append(f"... ({len(promoted)-30} more)")
    out.append("")

    # Regressions: candidates that went DOWN in category or score
    regressions = [(b, p) for b, p in matched_pairs
                   if b["category"] == "bona_fide" and p["category"] != "bona_fide"]
    out.append("## Regressions")
    out.append("")
    out.append(f"**{len(regressions)}** candidates demoted from bona-fide → something else.")
    if regressions:
        for b, p in regressions:
            out.append(f"  - {b['genome_id']} {b['contig_uid'][:14]} [{b['start0']}:{b['end0']}) "
                       f"was bona_fide({float(b['score']):.3f}), now {p['category']}({float(p['score']):.3f})")
    out.append("")

    # HQ / Reserve seed library impact
    hq = sum(1 for r in patched if r["qc_tier"] == "hq")
    reserve = sum(1 for r in patched if r["qc_tier"] == "reserve")
    out.append("## Seed library impact (patched)")
    out.append("")
    out.append(f"- HQ seeds: **{hq}**")
    out.append(f"- Reserve seeds: **{reserve}**")
    out.append("")

    report_path = PATCHED / "pre_post_comparison.md"
    report_path.write_text("\n".join(out) + "\n")
    print(f"wrote {report_path}  (promoted={len(promoted)}, regressions={len(regressions)}, "
          f"matched_pairs={len(matched_pairs)}, HQ={hq}, reserve={reserve})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
