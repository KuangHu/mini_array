#!/usr/bin/env python3
"""Export v2 pairs and reviewed candidate loci, split by repeat completeness.

The old all-array catalogs remain frozen. This exporter consumes only v2
pair evidence and v2 structural review; it never relabels historical calls.
All coordinates and DNA are genomic-forward; alignment strand is evidence,
not an inferred biological transcription direction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path

from step2 import load_cfg, require_v2_run, sha256
from step2_lib.pair import METHOD_VERSION, normalized_record
from step2b_lib.flanks import StagedGenomeCache, load_contig_lookup


def read_table(path):
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_table(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def pair_key(row):
    return (row["contig"], row["strand"], row["anchor_qseqid"], row["second_qseqid"],
            *(int(row[f"{role}_{bound}0"]) for role in ("anchor", "second") for bound in ("start", "end")))


REPEAT_FIELDS = [f"repeat_{i}_{field}" for i in (1, 2) for field in
                 ("role_in_pair", "query_id", "start0", "end0", "pident", "cov",
                  "qstart", "qend", "qlen", "gapopen", "completeness", "sequence_quality")]
PAIR_FIELDS = ["pair_id", "contig_uid", "genome_id", "genome_uid", "alignment_strand",
               "locus_id", "context", "method_version", *REPEAT_FIELDS,
               "spacer_start0", "spacer_end0", "spacer_len", "pair_kind", "edge_excluded"]
LOCUS_FIELDS = ["physical_locus_id", "genome_id", "genome_uid", "contig_uid", "start0", "end0",
                "strand", "strand_hypotheses", "structural_class", "pair_kinds", "n_mini_array_pairs",
                "n_half_repeat_pairs", "pair_ids", "left_flank_bp", "right_flank_bp"]
INDEX_FIELDS = ["record_id", "record_type", "pair_id", "physical_locus_id", "pair_kind", "contig_uid",
                "start0", "end0", "length", "sha256", "has_ambiguous_bases"]


def export_sequences(out, stem, pairs, loci, cache):
    """Write genomic-forward physical/pair/unit sequences with one hash per record."""
    index = []
    seen = set()
    with (out / f"{stem}.fna").open("w", encoding="utf-8") as fh:
        def emit(r, record_id, record_type, start, end):
            seq = cache.get_contig_seq(r["genome_uid"], r["contig_uid"])
            if record_id in seen or not 0 <= start < end <= len(seq):
                raise ValueError(f"Duplicate ID or invalid sequence coordinates: {record_id}")
            seen.add(record_id)
            dna = seq[start:end]
            fh.write(f">{record_id}\n{dna}\n")
            index.append(dict(record_id=record_id, record_type=record_type,
                              pair_id=r.get("pair_id", ""), physical_locus_id=r.get("physical_locus_id", ""),
                              pair_kind=r.get("pair_kind", r.get("pair_kinds", "")), contig_uid=r["contig_uid"],
                              start0=start, end0=end, length=len(dna), sha256=hashlib.sha256(dna.encode()).hexdigest(),
                              has_ambiguous_bases=str(bool(set(dna.upper()) - set("ACGT"))).lower()))
        for r in loci:
            emit(r, r["physical_locus_id"], "physical_locus", r["start0"], r["end0"])
        for r in pairs:
            base = (r.get("physical_locus_id", "") + "|" if r.get("physical_locus_id") else "") + r["pair_id"]
            emit(r, base, "pair", r["repeat_1_start0"], r["repeat_2_end0"])
            for unit in ("repeat_1", "spacer", "repeat_2"):
                emit(r, base + "|" + unit, unit, r[unit + "_start0"], r[unit + "_end0"])
    write_table(out / f"{stem}.index.tsv", index, INDEX_FIELDS)


def build(cfg):
    root = cfg["project_root"]
    step2 = root / cfg["output_dir"]
    step2c = root / cfg["review_c_dir"]
    require_v2_run(step2, cfg)
    status = json.loads((step2c / "step2c_status.json").read_text())
    if (status.get("method_version") != METHOD_VERSION or "review" not in status.get("phases_completed", [])
            or status.get("step2_manifest_sha256") != sha256(step2 / "run_manifest.json")
            or status.get("review_flank_bp") != cfg["edge"]["review_flank_bp"]):
        raise ValueError("A matching completed v2 structural review is required")
    for name in ("physical_loci.tsv", "physical_candidate_evidence.tsv"):
        if sha256(step2c / name) != status["outputs_sha256"][name]:
            raise ValueError(f"Structural review changed since completion: {name}")
    out = root / cfg["catalog_dir"]
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError(f"Catalog directory must be empty: {out}")
    lookup = load_contig_lookup(root / "results/staging/contigs.1000.tsv")
    cache = StagedGenomeCache(root / cfg["staged_dir"])
    loci2 = {r["locus_id"]: r for r in read_table(step2 / "loci.tsv")}
    attr = {}
    for r in read_table(step2 / "loci_attributions.tsv"):
        key = pair_key(r)
        if key in attr and attr[key] != r["locus_id"]:
            raise ValueError("Pair attributed to multiple Step 2 loci")
        attr[key] = r["locus_id"]
    pairs, by_locus = [], defaultdict(list)
    for i, raw in enumerate(read_table(step2 / "pairs_real.tsv"), 1):
        nr = normalized_record(raw)
        if nr["edge_excluded"]:
            raise ValueError("Edge pair in accepted input")
        meta = lookup[raw["contig"]]
        seq = cache.get_contig_seq(meta["genome_uid"], raw["contig"])
        if len(seq) != meta["length"] or len(seq) != int(raw["contig_length"]):
            raise ValueError("Contig length mismatch")
        lid = attr[pair_key(raw)]
        row = dict(pair_id=f"pair_{i:07d}", contig_uid=raw["contig"], genome_id=meta["genome_id"],
                   genome_uid=meta["genome_uid"], alignment_strand=raw["strand"], locus_id=lid,
                   context=loci2[lid]["context"], **nr)
        pairs.append(row)
        by_locus[lid].append(row)
    physical = {r["physical_locus_id"]: r for r in read_table(step2c / "physical_loci.tsv")}
    reviewed = read_table(step2c / "physical_candidate_evidence.tsv")
    candidates, candidate_pairs, excluded = [], [], []
    for ev in reviewed:
        s, e, length = (int(ev[k]) for k in ("span_start0", "span_end0", "contig_length"))
        flank = cfg["edge"]["review_flank_bp"]
        if not 0 <= s < e <= length:
            raise ValueError("Invalid reviewed locus coordinates")
        if ev["structural_class"] == "edge_excluded" or min(s, length - e) < flank:
            excluded.append(ev)
            continue
        if ev["structural_class"] != "candidate_mini_array":
            continue
        pid = ev["physical_locus_id"]
        phys = physical[pid]
        supporting = {r["pair_id"]: r for lid in phys["original_locus_ids"].split(",") for r in by_locus[lid]}
        if not supporting:
            raise ValueError(f"No pair evidence for {pid}")
        for r in supporting.values():
            if r["contig_uid"] != ev["contig_uid"] or not s <= r["repeat_1_start0"] < r["repeat_2_end0"] <= e:
                raise ValueError("Pair outside its attributed physical locus")
        counts = Counter(r["pair_kind"] for r in supporting.values())
        candidates.append(dict(physical_locus_id=pid, genome_id=ev["genome_id"], genome_uid=ev["genome_uid"],
                               contig_uid=ev["contig_uid"], start0=s, end0=e, strand=".",
                               strand_hypotheses=phys["strand_hypotheses"], structural_class=ev["structural_class"],
                               pair_kinds=",".join(sorted(counts)), n_mini_array_pairs=counts["mini_array"],
                               n_half_repeat_pairs=counts["half_repeat"], pair_ids=",".join(supporting),
                               left_flank_bp=s, right_flank_bp=length-e))
        candidate_pairs.extend(dict(r, physical_locus_id=pid) for r in supporting.values())
    write_table(out / "all_pairs.tsv", pairs, PAIR_FIELDS)
    export_sequences(out, "all_pairs", pairs, [], cache)
    write_table(out / "candidate_loci.tsv", candidates, LOCUS_FIELDS)
    write_table(out / "candidate_pairs.tsv", candidate_pairs, PAIR_FIELDS + ["physical_locus_id"])
    for kind, stem in (("mini_array", "mini_arrays"), ("half_repeat", "half_repeats")):
        selected = [r for r in candidate_pairs if r["pair_kind"] == kind]
        ids = {r["physical_locus_id"] for r in selected}
        physical_view = [r for r in candidates if r["physical_locus_id"] in ids]
        write_table(out / f"{stem}.tsv", physical_view, LOCUS_FIELDS)
        write_table(out / f"{stem}.pairs.tsv", selected, PAIR_FIELDS + ["physical_locus_id"])
        export_sequences(out, stem, selected, physical_view, cache)
        with (out / f"{stem}.bed").open("w") as fh:
            for r in physical_view:
                fh.write(f"{r['contig_uid']}\t{r['start0']}\t{r['end0']}\t{r['physical_locus_id']}\t0\t.\n")
    edge_fields = list(reviewed[0]) if reviewed else ["physical_locus_id", "structural_class"]
    write_table(out / "edge_excluded.review.tsv", excluded, edge_fields)
    shutil.copyfile(step2 / "pairs_real.edge_excluded.tsv", out / "edge_excluded.pairs.tsv")
    counts = Counter(r["pair_kind"] for r in candidate_pairs)
    summary = {"method_version": METHOD_VERSION, "step2_screened_pairs": len(pairs),
               "step2_min_flank_bp": cfg["edge"]["min_flank_bp"],
               "candidate_min_flank_bp": cfg["edge"]["review_flank_bp"],
               "candidate_physical_loci_union": len(candidates), "candidate_pair_counts": dict(counts),
               "mini_array_loci": sum(r["n_mini_array_pairs"] > 0 for r in candidates),
               "half_repeat_loci": sum(r["n_half_repeat_pairs"] > 0 for r in candidates),
               "loci_with_both_kinds": sum(r["n_mini_array_pairs"] > 0 and r["n_half_repeat_pairs"] > 0 for r in candidates),
               "review_edge_excluded": len(excluded),
               "step2_edge_excluded": len(read_table(step2 / "pairs_real.edge_excluded.tsv")),
               "inputs_sha256": {str(f): sha256(f) for f in [step2 / "run_manifest.json", step2c / "step2c_status.json"]},
               "outputs_sha256": {f.name: sha256(f) for f in sorted(out.iterdir()) if f.is_file()}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    from _slurm_gate import require_slurm
    require_slurm()
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True)
    args = ap.parse_args()
    print(json.dumps(build(load_cfg(args.config)), indent=2))


if __name__ == "__main__":
    main()
