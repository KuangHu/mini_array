#!/usr/bin/env python3
"""Step 2C analyze — phases 4..8 after step2c prepare has:
  - run physical dedup (physical_loci.tsv + physical_candidate_evidence.tsv)
  - recovered proteomes (step2c_pilot/proteins/*.faa)

Phases here:
  4. extract ±10 kb CDS per physical locus → candidate_all_cds.tsv
  5. cluster neighborhood CDS by blastp all-vs-all → protein_family_members.tsv
  6. BLAST family reps against all proteomes → homolog_neighborhoods.tsv
  7. intersect homologs with Step 2 anchor hits → rna_protein_associations.tsv
  8. updated candidate_cards.md + REPORT_STEP2C.md
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
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
from step2c_lib import proteins as prot_mod  # noqa: E402

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
OUT = PROJECT / "results" / "step2c_v2"
STEP2 = PROJECT / "results" / "step2_v2"

NEIGHBORHOOD_BP = 10000
MIN_PIDENT_FAM = 0.30
MIN_COV_FAM = 0.70
MIN_PIDENT_HOMO = 0.30
MIN_COV_HOMO = 0.70


def _log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def _parse_prodigal_faa(faa: Path, genome_id: str
                        ) -> tuple[dict[str, dict], dict[str, str]]:
    """Return (orfs_by_contig, protein_seqs). ORF is dict(start,end,strand,pos,
    protein_id); protein_seqs maps protein_id → sequence."""
    orfs: dict[str, list[dict]] = defaultdict(list)
    seqs: dict[str, str] = {}
    for hdr, seq in prot_mod.read_fasta(faa):
        m = prot_mod.parse_prodigal_header(hdr)
        if "start" not in m or "end" not in m:
            continue
        pid = f"{genome_id}|{m['contig']}|{m['pos']}"
        orfs[m["contig"]].append({
            "pos": m["pos"], "start1": m["start"], "end1": m["end"],
            "strand": m["strand"], "protein_id": pid,
            "length_aa": max(0, (m["end"] - m["start"] + 1) // 3 - 1),
        })
        seqs[pid] = seq.rstrip("*")
    for v in orfs.values():
        v.sort(key=lambda x: x["start1"])
    return orfs, seqs


def phase_neighborhood_cds(phys_rows: list[dict]) -> tuple[list[dict], dict[str, str]]:
    _log(f"[K4a] extracting ±{NEIGHBORHOOD_BP//1000} kb CDS for "
          f"{len(phys_rows)} physical loci")
    # Index phys by genome
    by_gid: dict[str, list[dict]] = defaultdict(list)
    for r in phys_rows:
        by_gid[r["genome_id"]].append(r)
    all_rows = []
    all_proteins: dict[str, str] = {}
    for gid, lst in by_gid.items():
        guid = gid.replace(".", "_")
        faa = OUT / "proteins" / f"{guid}.faa"
        if not faa.is_file():
            continue
        orfs, seqs = _parse_prodigal_faa(faa, gid)
        for r in lst:
            cuid = r["contig_uid"]
            s0 = int(r["span_start0"])
            e0 = int(r["span_end0"])
            ws, we = s0 - NEIGHBORHOOD_BP, e0 + NEIGHBORHOOD_BP
            for o in orfs.get(cuid, []):
                o_s0 = o["start1"] - 1
                o_e0 = o["end1"]
                if o_e0 < ws or o_s0 > we:
                    continue
                d = max(0, min(abs(o_s0 - e0), abs(s0 - o_e0)))
                all_rows.append({
                    "physical_locus_id": r["physical_locus_id"],
                    "genome_id": gid,
                    "contig_uid": cuid,
                    "protein_id": o["protein_id"],
                    "pos": o["pos"],
                    "start0": o_s0,
                    "end0": o_e0,
                    "strand": o["strand"],
                    "length_aa": o["length_aa"],
                    "distance_bp_to_candidate": d,
                })
                all_proteins[o["protein_id"]] = seqs[o["protein_id"]]
    _log(f"[K4a] CDS collected: {len(all_rows)} rows, "
          f"unique proteins: {len(all_proteins)}")
    return all_rows, all_proteins


def phase_cluster_families(proteins: dict[str, str], threads: int) -> list[fam_mod.ProteinFamily]:
    _log("[K4b] clustering neighborhood proteins via blastp all-vs-all")
    env_bin = PROJECT / "envs" / "identify" / "bin"
    faa_combined = OUT / "candidate_neighborhood_proteins.faa"
    with faa_combined.open("w") as fh:
        for pid, seq in sorted(proteins.items()):
            fh.write(f">{pid}\n{seq}\n")
    db_prefix = OUT / "nbr_db" / "nbr"
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    fam_mod.build_db(faa_combined, db_prefix, env_bin)
    out_tsv = OUT / "nbr_blastp_all_vs_all.tsv"
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        fam_mod.blastp(faa_combined, db_prefix, out_tsv, env_bin, threads=threads)
    families = fam_mod.cluster_connected_components(
        out_tsv, sorted(proteins), min_pident=MIN_PIDENT_FAM, min_cov=MIN_COV_FAM
    )
    _log(f"[K4b] protein families: {len(families)} "
          f"(singletons: {sum(1 for f in families if f.n_members == 1)})")
    return families


def phase_homologs(families: list[fam_mod.ProteinFamily],
                    proteins: dict[str, str], threads: int) -> Path:
    _log("[K5a] searching for homologs in all 1000 proteomes")
    env_bin = PROJECT / "envs" / "identify" / "bin"
    # Representative FASTA (longest member of each family)
    reps_fa = OUT / "family_representatives.faa"
    # pick longest member
    with reps_fa.open("w") as fh:
        for f in families:
            best = max(f.members, key=lambda m: len(proteins[m]))
            seq = proteins[best]
            fh.write(f">{f.family_id}|{best}\n{seq}\n")

    # Concatenate all proteomes into one DB (reuse if present)
    all_db_fna = OUT / "all_proteomes.faa"
    if not all_db_fna.is_file():
        with all_db_fna.open("w") as outfh:
            for faa in sorted((OUT / "proteins").glob("*.faa")):
                guid = faa.stem
                with faa.open() as inp:
                    for ln in inp:
                        if ln.startswith(">"):
                            # Prefix header with genome id so sseqid is unique
                            # keep the original header text after a |
                            outfh.write(f">{guid}|{ln[1:]}")
                        else:
                            outfh.write(ln)
    db_prefix = OUT / "all_db" / "all"
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    if not (db_prefix.with_suffix(".phr").is_file()
             or Path(f"{db_prefix}.phr").is_file()):
        fam_mod.build_db(all_db_fna, db_prefix, env_bin)
    out_tsv = OUT / "homologs_raw.tsv"
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        homolog_search.blastp_homologs(reps_fa, db_prefix, out_tsv,
                                         env_bin, threads=threads)
    _log(f"[K5a] homolog BLAST done ({out_tsv.stat().st_size} bytes)")
    return out_tsv


def phase_rna_protein_associations(phys_rows: list[dict],
                                    families: list[fam_mod.ProteinFamily],
                                    homolog_tsv: Path) -> dict:
    """For each family, count homologs across genomes, and intersect each
    homolog location with Step 2 anchor hits (hits_real.tsv).

    Produces:
      - homolog_neighborhoods.tsv  (one row per homolog hit)
      - rna_protein_associations.tsv (one row per family × genome summary)
    """
    _log("[K5b] parsing homolog hits + Step 2 anchor hits for intersection")
    # Load Step 2 anchor hits — only the ones that formed pairs we actually
    # care about = pairs_real.tsv (one record per anchor).
    anchors_by_contig: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    pairs_tsv = STEP2 / "pairs_real.tsv"
    if pairs_tsv.is_file():
        with pairs_tsv.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                anchors_by_contig[r["contig"]].append(
                    (int(r["anchor_start0"]), int(r["anchor_end0"]),
                     r["strand"])
                )
    _log(f"[K5b] anchors indexed: "
          f"{sum(len(v) for v in anchors_by_contig.values())}")

    # Also need genome_id for each contig (to tag which genome each anchor is in)
    contigs = {}
    with (PROJECT / "results" / "staging" / "contigs.1000.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            contigs[row["contig_uid"]] = row

    # Family lookup
    sseq_to_family: dict[str, str] = {}
    for f in families:
        for m in f.members:
            sseq_to_family[m] = f.family_id

    def _parse_sseqid(sseqid: str) -> tuple[str, str, int] | None:
        # format <guid>|<prodigal_header>  where prodigal header begins with the
        # original contig_X_N ...
        if "|" not in sseqid:
            return None
        guid, rest = sseqid.split("|", 1)
        parts = rest.split(" # ")
        first = parts[0]
        idx = first.rfind("_")
        if idx < 0:
            return None
        contig = first[:idx]
        try:
            pos = int(first[idx + 1:])
            s1 = int(parts[1]); e1 = int(parts[2]); strand = int(parts[3])
            return (guid, contig, pos, s1, e1, strand)
        except (ValueError, IndexError):
            return None

    nbr_rows = []
    assoc_counter: dict[tuple[str, str], dict] = defaultdict(
        lambda: {"n_homologs": 0, "n_with_rna_nearby": 0}
    )
    with homolog_tsv.open() as fh:
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            qid = p[0]
            # qid format: <family_id>|<protein_id>
            family_id = qid.split("|")[0]
            try:
                pident = float(p[2]); length = int(p[3])
                qlen = int(p[4]); slen = int(p[5])
                bitscore = float(p[7])
            except ValueError:
                continue
            cov = length / qlen if qlen else 0
            if pident / 100.0 < MIN_PIDENT_HOMO or cov < MIN_COV_HOMO:
                continue
            parsed = _parse_sseqid(p[1])
            if parsed is None:
                continue
            guid, contig, pos, s1, e1, strand = parsed
            # Reverse-map guid -> gid via contigs
            cmeta = contigs.get(contig)
            if cmeta is None or cmeta["genome_uid"] != guid:
                # cross-genome mismatch is possible if headers collide — skip
                continue
            gid = cmeta["genome_id"]
            # Check if an anchor hit exists on this contig within ±neighborhood_bp
            s0, e0 = s1 - 1, e1
            found_anchor = False
            min_dist = None
            for (a0, a1, strand_a) in anchors_by_contig.get(contig, []):
                dist = min(abs(a0 - e0), abs(s0 - a1))
                if dist <= NEIGHBORHOOD_BP:
                    found_anchor = True
                    if min_dist is None or dist < min_dist:
                        min_dist = dist
            nbr_rows.append({
                "family_id": family_id,
                "query_protein": qid.split("|", 1)[1] if "|" in qid else qid,
                "subject_genome_id": gid,
                "subject_contig": contig,
                "subject_protein_pos": pos,
                "subject_strand": strand,
                "subject_start0": s0,
                "subject_end0": e0,
                "pident": pident,
                "cov": round(cov, 3),
                "bitscore": bitscore,
                "rna_anchor_within_10kb": "true" if found_anchor else "false",
                "distance_to_nearest_anchor_bp": min_dist if min_dist is not None else "",
            })
            key = (family_id, gid)
            assoc_counter[key]["n_homologs"] += 1
            if found_anchor:
                assoc_counter[key]["n_with_rna_nearby"] += 1

    with (OUT / "homolog_neighborhoods.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh,
                           fieldnames=list(nbr_rows[0].keys()) if nbr_rows else
                           ["family_id", "query_protein", "subject_genome_id",
                             "subject_contig", "subject_protein_pos",
                             "subject_strand", "subject_start0", "subject_end0",
                             "pident", "cov", "bitscore",
                             "rna_anchor_within_10kb",
                             "distance_to_nearest_anchor_bp"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in nbr_rows:
            w.writerow(r)

    # Summary per family × genome
    assoc_rows = []
    for (fam, gid), info in sorted(assoc_counter.items()):
        assoc_rows.append({
            "family_id": fam,
            "genome_id": gid,
            "n_homologs": info["n_homologs"],
            "n_homologs_with_rna_anchor_within_10kb": info["n_with_rna_nearby"],
        })
    with (OUT / "rna_protein_associations.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh,
                           fieldnames=["family_id", "genome_id", "n_homologs",
                                        "n_homologs_with_rna_anchor_within_10kb"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in assoc_rows:
            w.writerow(r)
    _log(f"[K5b] homolog rows: {len(nbr_rows)}; "
          f"family×genome summaries: {len(assoc_rows)}")
    return {"n_homolog_rows": len(nbr_rows), "n_assoc_rows": len(assoc_rows)}


def main(argv: list[str]) -> int:
    global PROJECT, OUT, STEP2
    require_slurm()
    p = argparse.ArgumentParser(prog="step2c_analyze.py")
    p.add_argument("--config", type=Path, default=HERE.parent / "config/step2.yaml")
    p.add_argument("--phases", nargs="+",
                   default=["cds", "families", "homologs", "associations"],
                   help="which phases to run")
    args = p.parse_args(argv)
    from step2 import load_cfg, require_v2_run, sha256
    cfg = load_cfg(args.config)
    PROJECT = cfg["project_root"]
    OUT, STEP2 = PROJECT / cfg["review_c_dir"], PROJECT / cfg["output_dir"]
    method = require_v2_run(STEP2)["method_version"]
    status = json.loads((OUT / "step2c_status.json").read_text())
    if status.get("method_version") != method or status.get("step2_manifest_sha256") != sha256(STEP2 / "run_manifest.json"):
        raise ValueError("Rerun v2 structural review before protein analysis")

    OUT.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))

    phys_rows = list(csv.DictReader(
        (OUT / "physical_candidate_evidence.tsv").open(), delimiter="\t"))
    phys_rows = [r for r in phys_rows if r["structural_class"] != "edge_excluded"]
    _log(f"[load] physical loci (post-review) = {len(phys_rows)}")

    if "cds" in args.phases:
        cds_rows, proteins = phase_neighborhood_cds(phys_rows)
        with (OUT / "candidate_all_cds.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh,
                               fieldnames=list(cds_rows[0].keys()) if cds_rows
                               else ["physical_locus_id", "genome_id",
                                      "contig_uid", "protein_id", "pos",
                                      "start0", "end0", "strand", "length_aa",
                                      "distance_bp_to_candidate"],
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in cds_rows:
                w.writerow(r)
        # Also save protein sequences for downstream
        with (OUT / "candidate_neighborhood_proteins.faa").open("w") as fh:
            for pid, seq in sorted(proteins.items()):
                fh.write(f">{pid}\n{seq}\n")
    else:
        proteins = {}
        with (OUT / "candidate_neighborhood_proteins.faa").open() as fh:
            cur_h = None; cur = []
            for ln in fh:
                if ln.startswith(">"):
                    if cur_h is not None:
                        proteins[cur_h] = "".join(cur)
                    cur_h = ln[1:].strip(); cur = []
                else:
                    cur.append(ln.strip())
            if cur_h is not None:
                proteins[cur_h] = "".join(cur)
        _log(f"[load] cached proteins = {len(proteins)}")

    families: list[fam_mod.ProteinFamily] = []
    if "families" in args.phases:
        families = phase_cluster_families(proteins, threads)
        with (OUT / "protein_family_members.tsv").open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["family_id", "n_members", "representative",
                         "member_protein_id"])
            for f in families:
                for m in f.members:
                    w.writerow([f.family_id, f.n_members, f.representative, m])
    else:
        # reload
        members = defaultdict(list)
        with (OUT / "protein_family_members.tsv").open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                members[r["family_id"]].append(r["member_protein_id"])
        for fid, mem in sorted(members.items()):
            families.append(fam_mod.ProteinFamily(
                family_id=fid, representative=mem[0],
                members=mem, n_members=len(mem),
            ))

    if "homologs" in args.phases:
        homolog_tsv = phase_homologs(families, proteins, threads)
    else:
        homolog_tsv = OUT / "homologs_raw.tsv"

    if "associations" in args.phases:
        summary = phase_rna_protein_associations(phys_rows, families, homolog_tsv)
    else:
        summary = {}

    (OUT / "step2c_analyze_status.json").write_text(
        json.dumps({"phases": args.phases,
                     "n_physical": len(phys_rows),
                     "n_families": len(families),
                     **(summary or {})}, indent=2) + "\n")
    _log("[done]")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
