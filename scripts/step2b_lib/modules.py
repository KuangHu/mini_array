"""Group candidate neighborhoods into module families.

A `module` = set of candidates sharing the same HMM-signature in their
±10 kb neighborhood. Candidates with "no_hmm_hits" are kept as a
single bucket but are still enumerated separately.

For each module we report:
  - member locus_ids
  - unique genomes
  - unique repeat_ids used (distinct queries that landed on these loci)
  - unique "phyla" (via manifest independent_group_id)
  - near-identical-genome duplicates flagged separately: pairs of loci from
    different genomes but with identical neighborhood_hmm_signature AND
    distance between candidate start0 and ORFs below a tight tolerance
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass


@dataclass
class Module:
    module_id: str
    hmm_signature: str
    n_members: int
    member_locus_ids: list[str]
    genomes: set[str]
    phyla: set[str]
    unique_queries: set[str]
    classes: dict[str, int]


def group_by_signature(candidate_rows: list[dict]) -> list[Module]:
    by_sig: dict[str, list[dict]] = defaultdict(list)
    for row in candidate_rows:
        by_sig[row["hmm_signature"]].append(row)
    out: list[Module] = []
    for idx, sig in enumerate(sorted(by_sig)):
        members = by_sig[sig]
        locus_ids = [m["locus_id"] for m in members]
        genomes = {m["genome_id"] for m in members}
        phyla = {m.get("phylum", "NA") for m in members}
        qset = set()
        for m in members:
            if m.get("query_set"):
                for q in m["query_set"].split(","):
                    if q.strip():
                        qset.add(q.strip())
        classes = defaultdict(int)
        for m in members:
            classes[m["structural_class"]] += 1
        out.append(Module(
            module_id=f"mod_{idx+1:04d}",
            hmm_signature=sig,
            n_members=len(members),
            member_locus_ids=locus_ids,
            genomes=genomes,
            phyla=phyla,
            unique_queries=qset,
            classes=dict(classes),
        ))
    return out


MODULE_COLS = ["module_id", "hmm_signature", "n_members",
                "n_genomes", "n_phyla", "n_unique_queries",
                "structural_class_distribution",
                "representative_locus_id"]


def write_modules_tsv(modules: list[Module], out_path):
    import csv
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MODULE_COLS)
        for m in sorted(modules, key=lambda x: (-x.n_members, x.module_id)):
            dist = ",".join(f"{k}:{v}" for k, v in sorted(m.classes.items()))
            w.writerow([
                m.module_id, m.hmm_signature, m.n_members,
                len(m.genomes), len(m.phyla), len(m.unique_queries),
                dist, m.member_locus_ids[0],
            ])
