#!/usr/bin/env python3
"""Step 2B — structural + neighborhood explainer for Step 2 candidate loci.

Pipeline:
  1. Load candidate loci (standalone + soft-context + contig-edge from
     Step 2's loci.tsv) and the Step 1 context layers.
  2. For each candidate, pull ±2 kb flank from staged FASTA and overlay
     identify arrays/units + cctyper arrays/ORFs/HMM hits.
  3. Assign a structural_class.
  4. For candidates that survive (not ordinary_array_fragment), pull ±10 kb
     ORF neighborhood with HMM-signature.
  5. Cluster by HMM signature → module families.
  6. Audit background bg_rep{0,1,2} using the same bucket pipeline + check
     the dinuc shuffle validity + similarity to original.
  7. Emit all deliverables.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step2b_lib import background_audit, classify, flanks, modules, neighborhood, overlay  # noqa: E402

try:
    import yaml
except ImportError:
    sys.stderr.write("PyYAML required — run inside envs/driver.\n")
    raise


FLANK_BP = 2000
NEIGHBORHOOD_BP = 10000


def load_step2_config(path: Path) -> dict:
    from step2 import load_cfg
    return load_cfg(path)


def main(argv: list[str]) -> int:
    from _slurm_gate import require_slurm
    require_slurm()
    ap = argparse.ArgumentParser(prog="step2b.py")
    ap.add_argument("--config", required=True, type=Path,
                    help="config/step2.yaml — reused for paths")
    args = ap.parse_args(argv)
    global FLANK_BP
    cfg = load_step2_config(args.config)
    FLANK_BP = cfg["edge"]["review_flank_bp"]
    p = cfg["project_root"]

    step2_dir = p / cfg["output_dir"]
    from step2 import require_v2_run, sha256
    run_info = require_v2_run(step2_dir, cfg)
    out_dir = p / cfg["review_b_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[load] contig map, identify + cctyper layers…")
    contigs = flanks.load_contig_lookup(p / "results" / "staging" / "contigs.1000.tsv")
    id_arrays = overlay.load_identify_arrays(p / cfg["arrays_all_tsv"])
    id_units = overlay.load_identify_units(
        p / "results" / "step1_patched_min_score_0.75" / "array_units.tsv")
    cct = overlay.load_cctyper_crisprs(p / "raw" / "cctyper")

    # Manifest for phylum
    manifest = {}
    with (p / "input" / "genomes.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            manifest[r["genome_id"]] = r

    # 1. Load candidate loci from Step 2
    loci = []
    with (step2_dir / "loci.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["context"] not in ("standalone", "possible_soft_context"):
                continue
            loci.append(r)
    print(f"[load] candidate loci: {len(loci)}")

    # 2. Explain each candidate
    staged_cache = flanks.StagedGenomeCache(p / "input" / "staged")
    rows = []
    cas_summary_cache = {}
    orf_cache = {}
    hmm_cache = {}

    for L in loci:
        cuid = L["contig"]
        if cuid not in contigs:
            raise ValueError(f"Missing contig metadata: {cuid}")
        meta = contigs[cuid]
        guid = meta["genome_uid"]
        gid = meta["genome_id"]
        start0 = int(L["start0"])
        end0 = int(L["end0"])
        try:
            seq = staged_cache.get_contig_seq(guid, cuid)
        except KeyError:
            raise ValueError(f"Missing staged sequence: {cuid}")
        window = flanks.extract_with_flanks(seq, start0, end0, FLANK_BP)
        if guid not in orf_cache:
            orf_cache[guid] = overlay.load_cctyper_genes(guid, p / "raw" / "cctyper")
            hmm_cache[guid] = overlay.load_cctyper_hmmer(guid, p / "raw" / "cctyper")
        ov = overlay.overlay_candidate(
            cuid, start0, end0,
            identify_arrays=id_arrays, identify_units=id_units,
            cctyper_crisprs_contig=cct.get((guid, cuid), []),
            cctyper_genes_contig=orf_cache[guid].get(cuid, []),
            cctyper_hmm_contig=hmm_cache[guid].get(cuid, []),
            flank_bp=FLANK_BP,
        )
        cl = classify.classify(start0, end0, ov,
                                left_flank_bp=window["left_flank_bp"],
                                right_flank_bp=window["right_flank_bp"],
                                flank_bp=FLANK_BP)

        # Neighborhood (±10kb) for candidates that are NOT ordinary_array_fragment
        if cl["structural_class"] not in ("ordinary_array_fragment", "edge_excluded"):
            all_orfs = orf_cache[guid].get(cuid, [])
            hmm_by_pos = defaultdict(list)
            for h in hmm_cache[guid].get(cuid, []):
                hmm_by_pos[h["Pos"]].append({
                    "hmm": h["Hmm"], "evalue": h["Eval"],
                    "score": h["score"], "cov_hmm": h.get("Cov_hmm"),
                })
            nb = neighborhood.neighborhood_for_candidate(
                cuid, start0, end0,
                neighborhood_bp=NEIGHBORHOOD_BP, orfs=all_orfs,
                hmm_by_pos=hmm_by_pos,
            )
            # Genome-level Cas presence (any operon, any confidence)
            if guid not in cas_summary_cache:
                cas_summary_cache[guid] = neighborhood.load_cas_operons_for_genome(
                    guid, p / "raw" / "cctyper")
            genome_cas = cas_summary_cache[guid]
            cas_types = sorted({co["prediction"] for co in genome_cas
                                if co["confidence"] == "high"
                                and co["prediction"] not in ("", "False", "None", None)})
            cas_subtypes_putative = sorted({co["prediction"] for co in genome_cas
                                             if co["confidence"] == "putative"
                                             and co["prediction"]
                                             not in ("", "False", "None", None)})
        else:
            nb = {"n_orfs_in_window": 0, "n_orfs_with_hmm": 0,
                   "hmm_families_in_window": [], "hmm_signature": "n/a (ordinary)",
                   "orfs": []}
            genome_cas = []
            cas_types = []
            cas_subtypes_putative = []

        phy = ""
        if gid in manifest:
            phy = manifest[gid].get("independent_group_id", "")

        rows.append({
            "locus_id": L["locus_id"],
            "genome_id": gid,
            "genome_uid": guid,
            "contig_uid": cuid,
            "contig_length": meta["length"],
            "phylum": phy,
            "candidate_start0": start0,
            "candidate_end0": end0,
            "candidate_length": end0 - start0,
            "strand": L["strand"],
            "n_pairs": L["n_pairs"],
            "n_unique_queries": L["n_unique_queries"],
            "query_set": L["query_set"],
            "pair_kinds": L["pair_kinds"],
            "strict_support": L["strict_support"],
            "step2_context": L["context"],
            "structural_class": cl["structural_class"],
            "structural_reasons": "|".join(cl["reasons"]),
            "n_units_inside_candidate": cl["n_units_inside"],
            "n_units_in_flank": cl["n_units_in_flank"],
            "n_identify_arrays_in_window": len(ov["nearby_identify_arrays"]),
            "n_cctyper_arrays_in_window": len(ov["nearby_cctyper_crisprs"]),
            "n_orfs_in_neighborhood": nb["n_orfs_in_window"],
            "n_orfs_with_hmm": nb["n_orfs_with_hmm"],
            "hmm_signature": nb["hmm_signature"],
            "hmm_families_in_window": ",".join(nb["hmm_families_in_window"]),
            "genome_cas_types_high": ",".join(cas_types) or "none",
            "genome_cas_types_putative": ",".join(cas_subtypes_putative) or "none",
            "left_flank_bp": window["left_flank_bp"],
            "right_flank_bp": window["right_flank_bp"],
        })

    print(f"[explain] evaluated {len(rows)} candidates")

    # Write candidate_evidence.tsv
    cols = list(rows[0].keys()) if rows else ["locus_id", "genome_id", "contig_uid", "candidate_start0", "candidate_end0", "strand", "structural_class", "pair_kinds"]
    with (out_dir / "candidate_evidence.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                            lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Protein links (one row per locus × HMM family for the kept ones)
    with (out_dir / "candidate_protein_links.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["locus_id", "hmm_family", "n_orfs_with_this_hmm"])
        for r in rows:
            if r["structural_class"] in ("ordinary_array_fragment", "edge_excluded"):
                continue
            counts = Counter(r["hmm_families_in_window"].split(","))
            for fam, n in counts.items():
                if fam:
                    w.writerow([r["locus_id"], fam, n])

    # Module families
    mods = modules.group_by_signature(
        [r for r in rows if r["structural_class"] not in ("ordinary_array_fragment", "edge_excluded")]
    )
    modules.write_modules_tsv(mods, out_dir / "candidate_module_families.tsv")

    # 3. Re-bucket background
    print("[background] re-bucketing bg_rep{0,1,2} …")
    # Build a context_fn mapping (contig, s0, e0) → bucket, reusing Step 2 context
    # Need: known_array_regions.bed + arrays.all.tsv for soft context, + contig lengths
    from step1_lib import stage as _stage  # noqa: E402
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from step2_lib import context as context_mod  # noqa: E402

    hard = context_mod.load_known_array_bed(p / cfg["known_array_bed"])
    soft = context_mod.load_soft_context(p / cfg["arrays_all_tsv"])
    contig_lengths = {k: v["length"] for k, v in contigs.items()}

    def _ctx(contig, s0, e0):
        return context_mod.classify_locus(contig, s0, e0, hard, soft, contig_lengths,
                                          cfg["edge"]["min_flank_bp"])

    bg_bucket = {}
    for name in ("real", "bg_rep0", "bg_rep1", "bg_rep2"):
        pairs_tsv = step2_dir / f"pairs_{name}.tsv"
        if name == "real":
            pairs_tsv = step2_dir / "pairs_real.tsv"
        if not pairs_tsv.is_file():
            bg_bucket[name] = {"n_pairs": 0}
            continue
        bg_bucket[name] = background_audit.rebucket_pairs(
            pairs_tsv, {}, _ctx,
        )
    with (out_dir / "background_by_bucket.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["replicate", "metric", "known_array_interior",
                     "possible_soft_context", "standalone", "contig_edge",
                     "suspected_terminal_extension", "n_pairs",
                     "n_unique_loci_proxy"])
        for name, info in bg_bucket.items():
            pbc = info.get("pair_bucket_counts", {}) or {}
            lbc = info.get("locus_bucket_counts", {}) or {}
            w.writerow([name, "pair_counts",
                         pbc.get("known_array_interior", 0),
                         pbc.get("possible_soft_context", 0),
                         pbc.get("standalone", 0),
                         pbc.get("contig_edge", 0),
                         pbc.get("suspected_terminal_extension", 0),
                         info.get("n_pairs", 0),
                         info.get("n_unique_loci_proxy", 0)])
            w.writerow([name, "locus_counts",
                         lbc.get("known_array_interior", 0),
                         lbc.get("possible_soft_context", 0),
                         lbc.get("standalone", 0),
                         lbc.get("contig_edge", 0),
                         lbc.get("suspected_terminal_extension", 0),
                         "",  # no pair counter here
                         info.get("n_unique_loci_proxy", 0)])

    # 4. Background query audit — shuffle validity + similarity to original
    audit_rows = []
    real_fa = p / cfg["queries_fna"]
    for rep in range(3):
        bg_fa = step2_dir / f"bg_queries_rep{rep}.fna"
        if not bg_fa.is_file():
            continue
        for ar in background_audit.audit_shuffle(real_fa, bg_fa):
            ar["replicate"] = f"bg_rep{rep}"
            audit_rows.append(ar)
    audit_cols = ["replicate", "qseqid", "mono_preserved", "dinuc_preserved",
                   "orig_len", "bg_len",
                   "edit_distance_vs_original", "edit_distance_vs_rc",
                   "min_edit_distance", "similarity_to_closest"]
    with (out_dir / "background_query_audit.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=audit_cols, delimiter="\t",
                            lineterminator="\n")
        w.writeheader()
        for ar in audit_rows:
            w.writerow(ar)

    # Trace rep1 anomalous pairs → query/contig
    rep1_pairs = step2_dir / "pairs_bg_rep1.tsv"
    rep1_trace: list[dict] = []
    if rep1_pairs.is_file():
        with rep1_pairs.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rep1_trace.append({
                    "anchor_qseqid": r["anchor_qseqid"],
                    "second_qseqid": r["second_qseqid"],
                    "contig": r["contig"],
                    "strand": r["strand"],
                    "anchor_start0": r["anchor_start0"],
                    "second_end0": r["second_end0"],
                    "spacer_len": r["spacer_len"],
                })
    (out_dir / "bg_rep1_pair_trace.json").write_text(
        json.dumps({"n_pairs": len(rep1_trace), "trace": rep1_trace[:50]},
                    indent=2) + "\n"
    )

    # Summary
    summary = {
        "method_version": run_info["method_version"],
        "step2_manifest_sha256": sha256(step2_dir / "run_manifest.json"),
        "candidate_evidence_sha256": sha256(out_dir / "candidate_evidence.tsv"),
        "review_flank_bp": FLANK_BP,
        "n_candidates_evaluated": len(rows),
        "structural_class_counts": dict(Counter(r["structural_class"] for r in rows)),
        "step2_context_counts": dict(Counter(r["step2_context"] for r in rows)),
        "n_candidates_with_hmm_hit_in_window": sum(1 for r in rows if r["n_orfs_with_hmm"] > 0),
        "n_candidates_with_genome_cas_high": sum(1 for r in rows if r["genome_cas_types_high"] != "none"),
        "n_module_families": len(mods),
        "largest_modules_top5": [
            {"module_id": m.module_id, "signature": m.hmm_signature,
              "n_members": m.n_members, "n_genomes": len(m.genomes)}
            for m in sorted(mods, key=lambda x: -x.n_members)[:5]
        ],
        "background_bucket": bg_bucket,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
