#!/usr/bin/env python3
"""Step 2 v2: coherent bidirectional pairing; edge records are audit-only."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
from collections import Counter
from pathlib import Path

import yaml
from step2_lib import context, dedup, dinuc_shuffle, pair, search


def load_cfg(path: Path) -> dict:
    cfg = yaml.safe_load(path.read_text())
    if cfg.get("schema_version") != 2:
        raise ValueError("Step 2 v2 requires schema_version: 2")
    cfg["project_root"] = Path(os.environ.get("PROJECT_ROOT", cfg["project_root"])).resolve()
    cfg["code_root"] = Path(cfg["code_root"]).resolve()
    return cfg


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append([line[1:].strip(), ""])
        elif line.strip():
            records[-1][1] += line.strip()
    return [(h, s) for h, s in records]


def _write_fasta(records, path):
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require_v2_run(directory: Path, cfg: dict | None = None) -> dict:
    meta = json.loads((directory / "run_manifest.json").read_text())
    if meta.get("method_version") != pair.METHOD_VERSION or meta.get("status") != "complete":
        raise ValueError("A completed Step 2 v2 run is required; rerun raw hits first")
    for name, digest in meta["outputs_sha256"].items():
        if sha256(directory / name) != digest:
            raise ValueError(f"Step 2 output changed since completion: {name}")
    if cfg is not None:
        for key in ("anchor", "second_repeat", "spacer", "edge", "locus_dedup_merge_bp",
                    "queries_fna", "strict_queries_fna", "known_array_bed", "arrays_all_tsv", "staged_dir"):
            if meta["config"][key] != cfg[key]:
                raise ValueError(f"Config differs from the completed Step 2 run: {key}")
        for name, digest in meta["input_sha256"].items():
            if sha256(cfg["project_root"] / name) != digest:
                raise ValueError(f"Step 2 input changed since completion: {name}")
    return meta


def process_hits(hits, cfg, output: Path, name: str):
    """Identical rules for real and shuffled controls. Return accepted pairs only."""
    params = dict(anchor_min_cov=cfg["anchor"]["min_cov_query"],
                  anchor_min_id=cfg["anchor"]["min_identity"],
                  second_min_cov=cfg["second_repeat"]["min_cov_query"],
                  second_min_id=cfg["second_repeat"]["min_identity"])
    anchors, secondary = pair.classify_hits(hits, **params)
    found = pair.find_pairs(hits, **params,
                           spacer_min=cfg["spacer"]["min_len"],
                           spacer_max=cfg["spacer"]["max_len"],
                           min_flank_bp=cfg["edge"]["min_flank_bp"])
    accepted = [p for p in found if not p.edge_excluded]
    excluded = [p for p in found if p.edge_excluded]
    pair.write_pairs_tsv(accepted, output / f"pairs_{name}.tsv")
    pair.write_pairs_tsv(excluded, output / f"pairs_{name}.edge_excluded.tsv")
    kinds = Counter(p.pair_kind for p in accepted)
    for kind in ("mini_array", "half_repeat"):
        pair.write_pairs_tsv([p for p in accepted if p.pair_kind == kind],
                             output / f"pairs_{name}.{kind}.tsv")
    return accepted, dict(raw_hits=len(hits), anchors=len(anchors),
                         secondary_passing=len(secondary), pairs=len(accepted),
                         mini_array=kinds["mini_array"], half_repeat=kinds["half_repeat"],
                         edge_excluded=len(excluded))


def do_run(cfg: dict, *, run_background=True, hits_dir: Path | None = None) -> None:
    p = cfg["project_root"]
    out = p / cfg["output_dir"]
    if out.name == "step2_pilot":
        raise ValueError("Historical step2_pilot is frozen; select a new output directory")
    out.mkdir(parents=True, exist_ok=True)
    # Never combine new real pairs with stale controls or downstream artifacts.
    if any(out.iterdir()):
        raise ValueError(f"Output must be empty: {out}; choose a fresh run directory")
    manifest = {"method_version": pair.METHOD_VERSION, "status": "running",
                "config": {k: str(v) if isinstance(v, Path) else v for k, v in cfg.items()},
                "source_hits": {}, "replayed": hits_dir is not None,
                "input_sha256": {cfg[k]: sha256(p / cfg[k]) for k in
                                 ("queries_fna", "strict_queries_fna", "known_array_bed", "arrays_all_tsv")},
                "code_sha256": {str(f.relative_to(Path(__file__).parent)): sha256(f) for f in
                                [Path(__file__), *sorted((Path(__file__).parent / "step2_lib").glob("*.py"))]}}
    manifest_path = out / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    names = ["real"]
    if run_background:
        names += [f"bg_rep{i}" for i in range(cfg["background"]["n_replicates"])]
    real_queries = p / cfg["queries_fna"]
    if hits_dir is not None:
        for name in names:
            src = hits_dir / f"hits_{name}.tsv"
            shutil.copyfile(src, out / src.name)  # missing requested controls are an error
            manifest["source_hits"][name] = {"path": str(src.resolve()), "sha256": sha256(src)}
        for name in names[1:]:
            src = hits_dir / f"bg_queries_{name.removeprefix('bg_')}.fna"
            if src.is_file():
                shutil.copyfile(src, out / src.name)
    else:
        env_bin = p / "envs" / "identify" / "bin"
        db = out / "subject_db" / "contigs"
        search.build_subject_db(p / cfg["staged_dir"], db, env_bin)
        records = _read_fasta(real_queries)
        for name in names:
            queries = real_queries
            if name != "real":
                rep = int(name.removeprefix("bg_rep"))
                queries = out / f"bg_queries_rep{rep}.fna"
                _write_fasta(dinuc_shuffle.shuffle_fasta_records(records, cfg["background"]["seed"], rep), queries)
            dest = out / f"hits_{name}.tsv"
            search.run_blastn(queries, db, dest, env_bin,
                             threads=cfg["execution"]["cpus_per_job"],
                             word_size=cfg["blast"]["word_size"], evalue=cfg["blast"]["evalue"],
                             dust=cfg["blast"]["dust"], outfmt=cfg["blast"]["outfmt"])
            manifest["source_hits"][name] = {"path": str(dest), "sha256": sha256(dest)}
    stats = {}
    for name in names:
        accepted, stats[name] = process_hits(search.read_hits(out / f"hits_{name}.tsv"), cfg, out, name)
        if name == "real":
            pairs = accepted
    loci = dedup.dedup(pairs, cfg["locus_dedup_merge_bp"])
    hard = context.load_known_array_bed(p / cfg["known_array_bed"])
    soft = context.load_soft_context(p / cfg["arrays_all_tsv"])
    lengths = context.contig_lengths_from_staged(p / cfg["staged_dir"])
    for pp in pairs:
        if lengths.get(pp.contig) != pp.contig_length:
            raise ValueError(f"BLAST/staged contig length mismatch: {pp.contig}")
    for locus in loci:
        locus.context = context.classify_locus(locus.contig, locus.start0, locus.end0,
                                               hard, soft, lengths, cfg["edge"]["min_flank_bp"])
        if locus.context == "contig_edge":
            raise ValueError("Edge pair leaked into accepted loci")
    strict = {h.split()[0] for h, _ in _read_fasta(p / cfg["strict_queries_fna"])}
    dedup.write_loci_tsv(loci, out / "loci.tsv", strict)
    dedup.write_locus_attributions(loci, out / "loci_attributions.tsv")
    summary = {"method_version": pair.METHOD_VERSION, "pair_counts": stats,
               "loci_real": len(loci), "bucket_counts": dict(Counter(l.context for l in loci)),
               "strict_support_loci": sum(bool(l.queries & strict) for l in loci),
               "background_is_calibrated_fdr": False}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    manifest["status"] = "complete"
    manifest["outputs_sha256"] = {f.name: sha256(f) for f in sorted(out.glob("*.tsv"))}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main(argv=None):
    from _slurm_gate import require_slurm
    require_slurm()
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run")
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--no-background", action="store_true")
    run.add_argument("--hits-dir", type=Path, help="Replay raw BLAST HSPs; never replay old pairs.tsv")
    args = ap.parse_args(argv)
    do_run(load_cfg(args.config), run_background=not args.no_background, hits_dir=args.hits_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
