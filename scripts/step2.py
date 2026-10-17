#!/usr/bin/env python3
"""Step 2 pilot driver.

Pipeline:
  1. Read config/step2.yaml
  2. Build a combined BLAST DB from all staged *.fna (once)
  3. Run blastn with the real queries → hits_real.tsv
  4. (If background enabled) shuffle queries and run blastn → hits_bg_*.tsv
  5. Parse hits, apply anchor + secondary filters, find pairs → pairs_real.tsv
     (and pairs_bg_*.tsv)
  6. Dedup pair footprints → loci.tsv (+ attributions)
  7. Classify each locus by context mask (known_array_bed + arrays.all)
  8. Write loci_classified.tsv + SUMMARY.md
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step2_lib import context as context_mod  # noqa: E402
from step2_lib import dedup as dedup_mod  # noqa: E402
from step2_lib import dinuc_shuffle  # noqa: E402
from step2_lib import pair as pair_mod  # noqa: E402
from step2_lib import search  # noqa: E402

try:
    import yaml
except ImportError:
    sys.stderr.write("PyYAML required — run inside envs/driver.\n")
    raise


def load_cfg(path: Path) -> dict:
    cfg = yaml.safe_load(path.read_text())
    cfg["project_root"] = Path(cfg["project_root"]).resolve()
    cfg["code_root"] = Path(cfg["code_root"]).resolve()
    return cfg


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    cur_h = None
    cur_s = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_h is not None:
                    records.append((cur_h, "".join(cur_s)))
                cur_h = line[1:].strip()
                cur_s = []
            else:
                cur_s.append(line.strip())
        if cur_h is not None:
            records.append((cur_h, "".join(cur_s)))
    return records


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as fh:
        for h, s in records:
            fh.write(f">{h}\n{s}\n")


def _blastdb_env_bin(cfg: dict) -> Path:
    return cfg["project_root"] / "envs" / "identify" / "bin"


def do_run(cfg: dict, *, run_background: bool = True,
            keep_blastdb: bool = True) -> None:
    p = cfg["project_root"]
    out_dir = p / cfg["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = out_dir / "subject_db" / "contigs"
    env_bin = _blastdb_env_bin(cfg)

    log = out_dir / "step2.log"
    def _log(msg: str) -> None:
        print(msg)
        with log.open("a") as fh:
            fh.write(f"{time.strftime('%H:%M:%S')} {msg}\n")

    _log(f"env_bin={env_bin}, output={out_dir}")

    # 1. Build combined DB if not present
    db_meta = db_prefix.parent / ".built"
    if not db_meta.is_file():
        _log("[build_subject_db] concatenating staged FASTAs and running makeblastdb…")
        search.build_subject_db(p / cfg["staged_dir"], db_prefix, env_bin)
        db_meta.write_text("ok\n")
        _log("[build_subject_db] done")
    else:
        _log("[build_subject_db] reusing existing DB")

    # 2. Run blastn with real queries
    real_queries = p / cfg["queries_fna"]
    hits_real = out_dir / "hits_real.tsv"
    _log(f"[blastn real] {real_queries} → {hits_real}")
    search.run_blastn(real_queries, db_prefix, hits_real, env_bin,
                       threads=cfg["execution"]["cpus_per_job"],
                       word_size=cfg["blast"]["word_size"],
                       evalue=cfg["blast"]["evalue"],
                       dust=cfg["blast"]["dust"],
                       outfmt=cfg["blast"]["outfmt"])
    _log(f"[blastn real] done ({hits_real.stat().st_size} bytes)")

    # 3. (Optional) background
    bg_hits_files: list[Path] = []
    if run_background:
        bg_cfg = cfg.get("background", {}) or {}
        n_rep = int(bg_cfg.get("n_replicates", 0))
        if n_rep:
            seed = int(bg_cfg.get("seed", 42))
            records = _read_fasta(real_queries)
            for rep in range(n_rep):
                shuffled = dinuc_shuffle.shuffle_fasta_records(records, seed, rep)
                bg_fa = out_dir / f"bg_queries_rep{rep}.fna"
                _write_fasta(shuffled, bg_fa)
                bg_hits = out_dir / f"hits_bg_rep{rep}.tsv"
                _log(f"[blastn bg {rep}] {bg_fa.name} → {bg_hits.name}")
                search.run_blastn(bg_fa, db_prefix, bg_hits, env_bin,
                                   threads=cfg["execution"]["cpus_per_job"],
                                   word_size=cfg["blast"]["word_size"],
                                   evalue=cfg["blast"]["evalue"],
                                   dust=cfg["blast"]["dust"],
                                   outfmt=cfg["blast"]["outfmt"])
                bg_hits_files.append(bg_hits)
                _log(f"[blastn bg {rep}] done ({bg_hits.stat().st_size} bytes)")

    # 4. Parse hits and pair repeats
    hits = search.read_hits(hits_real)
    _log(f"[pairs real] {len(hits)} raw hits")
    anchors, secondary = pair_mod.classify_hits(
        hits,
        anchor_min_cov=cfg["anchor"]["min_cov_query"],
        anchor_min_id=cfg["anchor"]["min_identity"],
        second_min_cov=cfg["second_repeat"]["min_cov_query"],
        second_min_id=cfg["second_repeat"]["min_identity"],
    )
    pairs = pair_mod.find_pairs(
        hits,
        anchor_min_cov=cfg["anchor"]["min_cov_query"],
        anchor_min_id=cfg["anchor"]["min_identity"],
        second_min_cov=cfg["second_repeat"]["min_cov_query"],
        second_min_id=cfg["second_repeat"]["min_identity"],
        spacer_min=cfg["spacer"]["min_len"],
        spacer_max=cfg["spacer"]["max_len"],
        overlap_merge_bp=cfg["overlap_merge_bp"],
    )
    _log(f"[pairs real] {len(anchors)} anchors, {len(secondary)} secondary-passing, "
          f"{len(pairs)} pairs")
    pair_mod.write_pairs_tsv(pairs, out_dir / "pairs_real.tsv")

    # Background pair counts (for comparison only)
    bg_summary = []
    for bg_hits_file in bg_hits_files:
        bg_hits = search.read_hits(bg_hits_file)
        bg_anchors, _bg_sec = pair_mod.classify_hits(
            bg_hits,
            anchor_min_cov=cfg["anchor"]["min_cov_query"],
            anchor_min_id=cfg["anchor"]["min_identity"],
            second_min_cov=cfg["second_repeat"]["min_cov_query"],
            second_min_id=cfg["second_repeat"]["min_identity"],
        )
        bg_pairs = pair_mod.find_pairs(
            bg_hits,
            anchor_min_cov=cfg["anchor"]["min_cov_query"],
            anchor_min_id=cfg["anchor"]["min_identity"],
            second_min_cov=cfg["second_repeat"]["min_cov_query"],
            second_min_id=cfg["second_repeat"]["min_identity"],
            spacer_min=cfg["spacer"]["min_len"],
            spacer_max=cfg["spacer"]["max_len"],
            overlap_merge_bp=cfg["overlap_merge_bp"],
        )
        bg_pairs_tsv = out_dir / f"pairs_{bg_hits_file.stem.replace('hits_','')}.tsv"
        pair_mod.write_pairs_tsv(bg_pairs, bg_pairs_tsv)
        bg_summary.append({
            "replicate": bg_hits_file.stem,
            "raw_hits": len(bg_hits),
            "anchors": len(bg_anchors),
            "pairs": len(bg_pairs),
        })
        _log(f"[pairs {bg_hits_file.stem}] {len(bg_hits)} hits / {len(bg_anchors)} anchors / {len(bg_pairs)} pairs")

    # 5. Dedup to loci and attach context
    loci = dedup_mod.dedup(pairs, cfg["locus_dedup_merge_bp"])
    hard = context_mod.load_known_array_bed(p / cfg["known_array_bed"])
    soft = context_mod.load_soft_context(p / cfg["arrays_all_tsv"])
    contig_lengths = context_mod.contig_lengths_from_staged(p / cfg["staged_dir"])
    for L in loci:
        L.context = context_mod.classify_locus(
            L.contig, L.start0, L.end0, hard, soft, contig_lengths,
        )

    strict_qids = set()
    strict_fa = p / cfg["strict_queries_fna"]
    if strict_fa.is_file():
        for h, _ in _read_fasta(strict_fa):
            strict_qids.add(h.split()[0])

    dedup_mod.write_loci_tsv(loci, out_dir / "loci.tsv", strict_qids)
    dedup_mod.write_locus_attributions(loci, out_dir / "loci_attributions.tsv")

    # 6. Summarise
    bucket_counts = Counter(L.context for L in loci)
    summary = {
        "total_raw_hits": len(hits),
        "anchors": len(anchors),
        "secondary_passing": len(secondary),
        "pairs_real": len(pairs),
        "loci_real": len(loci),
        "genomes_with_locus": len({L.contig.rsplit("_", 1)[0] for L in loci}),
        "bucket_counts": dict(bucket_counts),
        "background_summary": bg_summary,
        "strict_support_loci": sum(1 for L in loci if any(q in strict_qids for q in L.queries)),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _log(f"[summary] loci={len(loci)}  buckets={dict(bucket_counts)}  "
          f"strict_support={summary['strict_support_loci']}  "
          f"bg_runs={len(bg_summary)}")

    print(json.dumps(summary, indent=2))


def main(argv: list[str]) -> int:
    from _slurm_gate import require_slurm
    require_slurm()
    p = argparse.ArgumentParser(prog="step2.py")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--config", required=True, type=Path)
    r.add_argument("--no-background", action="store_true")
    args = p.parse_args(argv)
    cfg = load_cfg(args.config)
    if args.cmd == "run":
        do_run(cfg, run_background=not args.no_background)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
