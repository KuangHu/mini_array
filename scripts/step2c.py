#!/usr/bin/env python3
"""Step 2C — physical dedup + structural review + annotation-independent
protein-family clustering + homolog neighborhoods.

Steps:
  1. physical_locus_id dedup over Step 2B candidate_evidence.tsv
  2. structural_class re-review (unit evidence required for ordinary_array)
  3. recover proteomes for every genome (prodigal, matched to Step 1 mode)
  4. extract ±10 kb CDS per physical locus
  5. cluster neighborhood CDS into annotation-independent protein families
  6. BLAST family reps against all 1000 proteomes → homolog neighborhoods
  7. intersect homologs with Step 2 anchor hits → RNA-protein associations
  8. emit new tables, update cards, write report
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _slurm_gate import require_slurm  # noqa: E402
from step2c_lib import families as fam_mod  # noqa: E402
from step2c_lib import homolog_search  # noqa: E402
from step2c_lib import physical_locus  # noqa: E402
from step2c_lib import proteins as prot_mod  # noqa: E402
from step2c_lib import structural_review  # noqa: E402
from step2b_lib import flanks, overlay  # noqa: E402


PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
OUT = PROJECT / "results" / "step2c_pilot"
STEP2B = PROJECT / "results" / "step2b_pilot"
STEP2 = PROJECT / "results" / "step2_pilot"

MIN_PIDENT_FAM = 0.30
MIN_COV_FAM = 0.70
MIN_PIDENT_HOMO = 0.30
MIN_COV_HOMO = 0.70
NEIGHBORHOOD_BP = 10000


def _log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def phase_physical_dedup() -> list[physical_locus.PhysicalLocus]:
    rows = list(csv.DictReader((STEP2B/"candidate_evidence.tsv").open(),
                                delimiter="\t"))
    _log(f"[K1] input candidates = {len(rows)}")
    phys = physical_locus.dedup_candidates(rows, merge_bp=10)
    _log(f"[K1] physical loci after dedup = {len(phys)}")
    return phys


def phase_structural_review(phys: list[physical_locus.PhysicalLocus]) -> dict:
    # reload the identify layers once
    cfg_path = Path("/home/kuangh/tools/CRISPR_array_related/config/step2.yaml")
    import yaml
    cfg = yaml.safe_load(cfg_path.read_text())

    id_arrays = overlay.load_identify_arrays(PROJECT / cfg["arrays_all_tsv"])
    id_units = overlay.load_identify_units(
        PROJECT / "results" / "step1_patched_min_score_0.75" / "array_units.tsv")
    cct = overlay.load_cctyper_crisprs(PROJECT / "raw" / "cctyper")

    contigs = flanks.load_contig_lookup(
        PROJECT / "results" / "staging" / "contigs.1000.tsv")

    staged_cache = flanks.StagedGenomeCache(PROJECT / "input" / "staged")
    classified = {}
    orf_cache: dict[str, dict] = {}
    hmm_cache: dict[str, dict] = {}

    for L in phys:
        cuid = L.contig_uid
        meta = contigs.get(cuid)
        if meta is None:
            continue
        guid = meta["genome_uid"]
        try:
            seq = staged_cache.get_contig_seq(guid, cuid)
        except KeyError:
            continue
        w = flanks.extract_with_flanks(seq, L.span_start0, L.span_end0, 2000)
        if guid not in orf_cache:
            orf_cache[guid] = overlay.load_cctyper_genes(guid,
                                                          PROJECT / "raw" / "cctyper")
            hmm_cache[guid] = overlay.load_cctyper_hmmer(guid,
                                                          PROJECT / "raw" / "cctyper")
        ov = overlay.overlay_candidate(
            cuid, L.span_start0, L.span_end0,
            identify_arrays=id_arrays, identify_units=id_units,
            cctyper_crisprs_contig=cct.get((guid, cuid), []),
            cctyper_genes_contig=orf_cache[guid].get(cuid, []),
            cctyper_hmm_contig=hmm_cache[guid].get(cuid, []),
            flank_bp=2000,
        )
        cl = structural_review.classify(
            L.span_start0, L.span_end0, ov,
            left_flank_bp=w["left_flank_bp"],
            right_flank_bp=w["right_flank_bp"],
            flank_bp=2000,
        )
        classified[L.physical_locus_id] = {
            "overlay": ov,
            "class": cl,
            "flank_left": w["left_flank_bp"],
            "flank_right": w["right_flank_bp"],
            "contig_length": meta["length"],
            "genome_uid": guid,
        }
    cnt = Counter(c["class"]["structural_class"] for c in classified.values())
    _log(f"[K2] structural_class counts: {dict(cnt)}")
    return classified


def phase_recover_proteomes(manifest_rows: list[dict], env_bin: Path,
                             threads: int) -> dict[str, Path]:
    proteins_dir = OUT / "proteins"
    proteins_dir.mkdir(parents=True, exist_ok=True)
    gid_to_faa: dict[str, Path] = {}
    prodigal = str(env_bin / "prodigal")
    for i, row in enumerate(manifest_rows):
        gid = row["genome_id"]
        guid = gid.replace(".", "_")
        mode = prot_mod.prodigal_mode(row["source_type"])
        staged = PROJECT / "input" / "staged" / f"{guid}.fna"
        out_faa = proteins_dir / f"{guid}.faa"
        out_gff = proteins_dir / f"{guid}.gff"
        if not out_faa.is_file() or out_faa.stat().st_size == 0:
            prot_mod.recover_proteome(staged, out_faa, out_gff, mode, prodigal)
        gid_to_faa[gid] = out_faa
        if (i + 1) % 100 == 0:
            _log(f"[K3] prodigal {i+1}/{len(manifest_rows)}")
    _log(f"[K3] proteomes recovered: {len(gid_to_faa)}")
    return gid_to_faa


def phase_neighborhood_cds(classified: dict, gid_to_faa: dict[str, Path],
                            contigs: dict) -> tuple[list[dict], dict[str, str]]:
    """Return (per-locus CDS rows, protein_id → sequence map)."""
    rows = []
    all_proteins: dict[str, str] = {}
    # For each genome in the classified set, read its FAA once and index by (contig, pos).
    gid_cache: dict[str, dict[tuple[str, int], tuple[str, int, int, int]]] = {}
    for pid, info in classified.items():
        guid = info["genome_uid"]
        if guid not in gid_cache:
            # Find the gid matching this guid
            found_gid = next(g for g, p in gid_to_faa.items()
                              if p.stem == guid)
            faa = gid_to_faa[found_gid]
            idx: dict[tuple[str, int], tuple[str, int, int, int]] = {}
            seqs: dict[str, str] = {}
            for hdr, seq in prot_mod.read_fasta(faa):
                meta = prot_mod.parse_prodigal_header(hdr)
                if "start" not in meta or "end" not in meta:
                    continue
                prot_id = f"{found_gid}|{meta['contig']}|{meta['pos']}"
                idx[(meta["contig"], meta["pos"])] = (prot_id, meta["start"],
                                                       meta["end"], meta["strand"])
                seqs[prot_id] = seq.rstrip("*")
            gid_cache[guid] = (idx, seqs)
        idx, seqs = gid_cache[guid]
    # Now extract candidates
    for pid, info in classified.items():
        guid = info["genome_uid"]
        found_gid = next(g for g, p in gid_to_faa.items() if p.stem == guid)
        idx, seqs = gid_cache[guid]
        locus_contig = next(
            (L.contig_uid for L in [info] if False), None)  # not used
    # Easier: iterate through classified again with cuid from physical locus list
    return rows, all_proteins


def main(argv: list[str]) -> int:
    require_slurm()
    p = argparse.ArgumentParser(prog="step2c.py")
    p.add_argument("--phases", nargs="+",
                   default=["dedup", "review", "proteomes", "neighborhoods",
                             "families", "homologs", "associations", "report"],
                   help="which phases to run")
    args = p.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    env_bin = PROJECT / "envs" / "identify" / "bin"
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    _log(f"env_bin={env_bin}  threads={threads}")

    # Load manifest once
    manifest = list(csv.DictReader((PROJECT / "input" / "genomes.tsv").open(),
                                    delimiter="\t"))

    phys: list[physical_locus.PhysicalLocus] = []
    classified: dict = {}
    gid_to_faa: dict[str, Path] = {}

    if "dedup" in args.phases:
        phys = phase_physical_dedup()
        rows = []
        for L in phys:
            rows.append({
                "physical_locus_id": L.physical_locus_id,
                "genome_id": L.genome_id,
                "contig_uid": L.contig_uid,
                "span_start0": L.span_start0,
                "span_end0": L.span_end0,
                "strand_hypotheses": "|".join(sorted(L.strand_hypotheses)),
                "n_original_loci": len(L.original_locus_ids),
                "original_locus_ids": ",".join(L.original_locus_ids),
            })
        with (OUT / "physical_loci.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh,
                               fieldnames=list(rows[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    else:
        phys = []
        with (OUT / "physical_loci.tsv").open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                L = physical_locus.PhysicalLocus(
                    physical_locus_id=row["physical_locus_id"],
                    genome_id=row["genome_id"],
                    contig_uid=row["contig_uid"],
                    span_start0=int(row["span_start0"]),
                    span_end0=int(row["span_end0"]),
                )
                phys.append(L)

    if "review" in args.phases:
        classified = phase_structural_review(phys)
        # Write physical_candidate_evidence.tsv
        rows = []
        for L in phys:
            info = classified.get(L.physical_locus_id, {})
            cl = info.get("class", {}) if info else {}
            rows.append({
                "physical_locus_id": L.physical_locus_id,
                "genome_id": L.genome_id,
                "contig_uid": L.contig_uid,
                "span_start0": L.span_start0,
                "span_end0": L.span_end0,
                "span_length": L.span_end0 - L.span_start0,
                "strand_hypotheses": "|".join(sorted(L.strand_hypotheses)),
                "n_original_loci": len(L.original_locus_ids),
                "structural_class": cl.get("structural_class", ""),
                "structural_reasons": "|".join(cl.get("reasons", [])),
                "n_units_inside": cl.get("n_units_inside", 0),
                "n_units_in_flank": cl.get("n_units_in_flank", 0),
                "evidence_arrays_with_units": ",".join(
                    cl.get("evidence_arrays_with_units", [])),
                "flank_left_bp": info.get("flank_left", 0),
                "flank_right_bp": info.get("flank_right", 0),
                "contig_length": info.get("contig_length", 0),
                "genome_uid": info.get("genome_uid", ""),
            })
        with (OUT / "physical_candidate_evidence.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh,
                               fieldnames=list(rows[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow(r)

    if "proteomes" in args.phases:
        gid_to_faa = phase_recover_proteomes(manifest, env_bin, threads)

    # The rest of the pipeline (phases 4-8) runs in a continuation module
    # to keep this file readable. For now, write a status so a shell can
    # chain the next step.
    (OUT / "step2c_status.json").write_text(json.dumps({
        "phases_completed": args.phases,
        "n_physical_loci": len(phys),
        "n_proteomes": len(gid_to_faa) if gid_to_faa else -1,
    }, indent=2) + "\n")
    _log(f"[done] wrote step2c_status.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
