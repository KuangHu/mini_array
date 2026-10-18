#!/usr/bin/env python3
"""Render candidate_cards.md — one card per structurally-surviving candidate
with genomic context, overlay evidence, ±10 kb neighborhood summary, and
module family assignment."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
STEP2B = PROJECT / "results" / "step2b_v2"


def main() -> int:
    global PROJECT, STEP2B
    import argparse
    from step2 import load_cfg, require_v2_run
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=Path(__file__).resolve().parents[1] / "config/step2.yaml")
    cfg = load_cfg(ap.parse_args().config)
    PROJECT, STEP2B = cfg["project_root"], cfg["project_root"] / cfg["review_b_dir"]
    require_v2_run(PROJECT / cfg["output_dir"])
    import sys as _s
    _s.path.insert(0, "/home/kuangh/tools/CRISPR_array_related/scripts")
    from _slurm_gate import require_slurm
    require_slurm()
    evidence = list(csv.DictReader((STEP2B/"candidate_evidence.tsv").open(),
                                     delimiter="\t"))
    modules = list(csv.DictReader((STEP2B/"candidate_module_families.tsv").open(),
                                    delimiter="\t"))

    # Map locus → module via the hmm_signature field
    sig_to_module = {m["hmm_signature"]: m["module_id"] for m in modules}

    survived = [r for r in evidence if r["structural_class"] not in ("ordinary_array_fragment", "edge_excluded")]
    survived.sort(key=lambda r: (
        0 if r["strict_support"] == "true" else 1,         # strict first
        -int(r["n_orfs_with_hmm"] or 0),
        r["genome_id"], r["locus_id"],
    ))

    md = []
    md.append("# Step 2B candidate cards\n")
    md.append(f"Rendered: {len(survived)} candidates that survived structural filter\n")
    md.append("Ordering: strict_support first, then by HMM-rich neighborhood, then by genome.\n")
    md.append("---\n")

    for r in survived:
        md.append(f"## {r['locus_id']} — {r['structural_class']}\n")
        md.append(f"- **genome**: `{r['genome_id']}` ({r['phylum']})")
        md.append(f"- **contig**: `{r['contig_uid']}` ({r['contig_length']} bp)")
        md.append(f"- **candidate**: [{r['candidate_start0']}, {r['candidate_end0']}), "
                   f"strand `{r['strand']}`, length {r['candidate_length']} bp")
        md.append(f"- **flank available**: left {r['left_flank_bp']} bp, "
                   f"right {r['right_flank_bp']} bp")
        md.append(f"- **step2 context (before explain)**: `{r['step2_context']}`")
        md.append(f"- **structural reasons**: {r['structural_reasons']}")
        md.append(f"- **pairs forming this locus**: {r['n_pairs']} "
                   f"({r['n_unique_queries']} unique queries); "
                   f"strict support `{r['strict_support']}`")
        md.append(f"- **ORFs in ±10 kb**: {r['n_orfs_in_neighborhood']} "
                   f"({r['n_orfs_with_hmm']} with HMM hits)")
        md.append(f"- **HMM signature**: `{r['hmm_signature'] or 'n/a'}`")
        if r["hmm_families_in_window"]:
            md.append(f"- **HMM families nearby**: `{r['hmm_families_in_window']}`")
        md.append(f"- **genome-wide Cas (high conf.)**: `{r['genome_cas_types_high']}`")
        md.append(f"- **genome-wide Cas (putative)**: `{r['genome_cas_types_putative']}`")
        md.append(f"- **module family**: `{sig_to_module.get(r['hmm_signature'], 'unassigned')}`")
        md.append(f"- **query repeats**: `{r['query_set']}`")
        md.append("")

    (STEP2B / "candidate_cards.md").write_text("\n".join(md))
    print(f"wrote candidate_cards.md with {len(survived)} entries")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
