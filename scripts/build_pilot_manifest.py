#!/usr/bin/env python3
"""Build the pilot manifest for step 1 from the GTDB R220 reps dump.

Policy:
  - Source  : /groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/
  - Scope   : bacteria (bac120_metadata.tsv), gtdb_representative == 't'.
  - Sampling: deterministic round-robin across GTDB phyla to N (default 100),
              seeded by `seed` from config/step1.yaml for reproducibility.
  - FASTA   : must exist as <accession>_genomic.fna under
              gtdb_genomes_reps_r220/allfilestogether/.

Writes <project_root>/input/genomes.tsv matching step1.md §1 columns:
  genome_id, fasta_path, source_type, taxon, sample_id,
  independent_group_id, completeness, contamination

Also writes <project_root>/provenance/pilot_manifest.build.json with the source
versions, file hashes, and the exact sampling parameters, so the pilot can be
regenerated identically.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import yaml

GTDB_ROOT = Path("/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest")
BAC_META = GTDB_ROOT / "bac120_metadata.tsv"
REPS_DIR = GTDB_ROOT / "gtdb_genomes_reps_r220" / "allfilestogether"


def strip_prefix(acc: str) -> str:
    """`GB_GCA_000008085.1` / `RS_GCF_...` → `GCA_000008085.1` / `GCF_...`."""
    if acc.startswith(("GB_", "RS_")):
        return acc[3:]
    return acc


def parse_phylum(tax: str) -> str:
    """`d__Bacteria;p__Pseudomonadota;c__...` → `p__Pseudomonadota`."""
    for token in tax.split(";"):
        token = token.strip()
        if token.startswith("p__"):
            return token
    return "p__NA"


def hash_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def load_config(path: Path) -> dict:
    with path.open() as fh:
        cfg = yaml.safe_load(fh)
    cfg["project_root"] = Path(cfg["project_root"]).resolve()
    return cfg


def classify_source_type(ncbi_genome_category: str) -> str:
    """Map GTDB ncbi_genome_category → step1.md source_type values.

    step1.md §1 enumerates: isolate, MAG, phage, plasmid, unknown.
    - 'derived from metagenome' / 'single cell' / 'environmental sample' → MAG
      (all are non-cultured drafts and want prodigal --meta per §8).
    - 'none' is NCBI's label for cultured isolates → isolate.
    """
    cat = (ncbi_genome_category or "").strip().lower()
    if "metagenome" in cat or "single cell" in cat or "environmental" in cat:
        return "MAG"
    if cat in ("none", ""):
        return "isolate"
    return "unknown"


def collect_candidates() -> dict[str, list[dict]]:
    """Return {phylum: [row, ...]} for bac120 representatives with existing FASTA."""
    buckets: dict[str, list[dict]] = defaultdict(list)
    seen_without_fasta = 0
    total_reps = 0
    with BAC_META.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("gtdb_representative") != "t":
                continue
            total_reps += 1
            acc = strip_prefix(row["accession"])
            fasta = REPS_DIR / f"{acc}_genomic.fna"
            if not fasta.is_file():
                seen_without_fasta += 1
                continue
            phylum = parse_phylum(row.get("gtdb_taxonomy", ""))
            buckets[phylum].append({
                "accession": acc,
                "fasta_path": str(fasta),
                "gtdb_taxonomy": row["gtdb_taxonomy"],
                "phylum": phylum,
                "source_type": classify_source_type(row.get("ncbi_genome_category", "")),
                "ncbi_genome_category": row.get("ncbi_genome_category", "none"),
                "checkm_completeness": row.get("checkm_completeness", "NA"),
                "checkm_contamination": row.get("checkm_contamination", "NA"),
            })
    sys.stderr.write(
        f"[info] bac120 reps scanned: {total_reps}; with FASTA: "
        f"{total_reps - seen_without_fasta}; missing FASTA: {seen_without_fasta}; "
        f"phyla: {len(buckets)}\n"
    )
    return buckets


def round_robin_sample(buckets: dict[str, list[dict]], n: int,
                      seed: int) -> list[dict]:
    """Deterministic round-robin across phyla, within-phylum shuffle by seed.

    `n` can be larger than the number of phyla — in that case the loop makes
    multiple passes, picking one additional member from each pool per pass.
    """
    rng = random.Random(seed)
    phyla = sorted(buckets.keys())
    pools = {p: list(buckets[p]) for p in phyla}
    for p in phyla:
        rng.shuffle(pools[p])
    picks: list[dict] = []
    while len(picks) < n:
        progressed = False
        for p in phyla:
            if not pools[p]:
                continue
            picks.append(pools[p].pop())
            progressed = True
            if len(picks) >= n:
                break
        if not progressed:
            break
    return picks


def write_manifest(picks: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["genome_id", "fasta_path", "source_type", "taxon", "sample_id",
            "independent_group_id", "completeness", "contamination"]
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for row in picks:
            w.writerow([
                row["accession"],            # genome_id = GCA accession
                row["fasta_path"],           # absolute path to .fna
                row["source_type"],          # isolate | MAG derived from ncbi_genome_category
                row["gtdb_taxonomy"],
                "NA",                        # sample_id not provided by GTDB
                row["phylum"],               # phylum is the sampling stratum
                row.get("checkm_completeness", "NA"),
                row.get("checkm_contamination", "NA"),
            ])


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--n", type=int, default=100,
                   help="pilot size (default 100)")
    p.add_argument("--out", type=Path, default=None,
                   help="override manifest output path (default: project_root/input/genomes.tsv)")
    args = p.parse_args(argv)

    cfg = load_config(args.config)
    seed = int(cfg.get("seed", 42))
    out_path = args.out or (cfg["project_root"] / cfg["manifest"])

    for required in (BAC_META, REPS_DIR):
        if not required.exists():
            raise SystemExit(f"GTDB source missing: {required}")

    buckets = collect_candidates()
    picks = round_robin_sample(buckets, args.n, seed)
    if len(picks) < args.n:
        sys.stderr.write(
            f"[warn] only {len(picks)} candidates after round-robin "
            f"(requested {args.n})\n"
        )
    write_manifest(picks, out_path)

    build = {
        "generator": "scripts/build_pilot_manifest.py",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(args.config.resolve()),
        "seed": seed,
        "requested_n": args.n,
        "written_n": len(picks),
        "gtdb_release": "r220",
        "sources": {
            "bac120_metadata.tsv": {
                "path": str(BAC_META),
                "size": BAC_META.stat().st_size,
                "sha256": hash_file(BAC_META),
            },
            "reps_dir": str(REPS_DIR),
        },
        "sampling": "round-robin across p__ strata, deterministic shuffle(seed)",
        "manifest": str(out_path),
        "phyla_sampled": sorted({r["phylum"] for r in picks}),
    }
    prov_dir = cfg["project_root"] / "provenance"
    prov_dir.mkdir(parents=True, exist_ok=True)
    (prov_dir / "pilot_manifest.build.json").write_text(
        json.dumps(build, indent=2) + "\n"
    )

    sys.stderr.write(
        f"[ok] manifest -> {out_path} ({len(picks)} rows, "
        f"{len(build['phyla_sampled'])} phyla)\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
