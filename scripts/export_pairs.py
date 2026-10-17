#!/usr/bin/env python3
"""Export every Step 2 pair as position + sequence.

Reads results/step2_pilot/pairs_real.tsv (15,472 rows). For each pair:
  - compute pair_id = step2_pair_<NNNNN>
  - extract genomic-forward DNA for repeat_1, spacer, repeat_2, and the
    full pair span [min(start), max(end)) directly from staged FASTA
  - classify its context bucket using the SAME context function Step 2
    used (known_array_interior / possible_soft_context / standalone /
    contig_edge / suspected_terminal_extension / standalone)
  - emit all_pairs.fna, all_pair_units.fna, all_pairs_index.tsv, README.md

Pure extraction. No RC, no re-identification, no new biology. sbatch gate
enforced.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _slurm_gate import require_slurm  # noqa: E402
from step2_lib import context as context_mod  # noqa: E402

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
STAGED_DIR = PROJECT / "input" / "staged"
STEP2 = PROJECT / "results" / "step2_pilot"
STEP1 = PROJECT / "results" / "step1_patched_min_score_0.75"
OUT = PROJECT / "results" / "sequence_catalog" / "pairs"


def _read_staged(path: Path) -> dict[str, str]:
    seqs = {}
    cur = None
    parts: list[str] = []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(parts)
                cur = ln[1:].split()[0]
                parts = []
            else:
                parts.append(ln.strip())
        if cur is not None:
            seqs[cur] = "".join(parts)
    return seqs


class _Cache:
    def __init__(self): self._c: dict[str, dict[str, str]] = {}
    def get(self, guid: str) -> dict[str, str]:
        if guid not in self._c:
            self._c[guid] = _read_staged(STAGED_DIR / f"{guid}.fna")
        return self._c[guid]


def main() -> int:
    require_slurm()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # Context lookup inputs
    contig_lookup = {}
    with (PROJECT / "results" / "staging" / "contigs.1000.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            contig_lookup[r["contig_uid"]] = r
    print(f"[load] contigs: {len(contig_lookup)}")

    hard = context_mod.load_known_array_bed(
        STEP1 / "known_array_regions.bed")
    soft = context_mod.load_soft_context(STEP1 / "arrays.all.tsv")
    contig_lengths = {k: int(v["length"]) for k, v in contig_lookup.items()}
    print(f"[load] hard mask contigs: {len(hard)}, soft mask contigs: {len(soft)}")

    pairs_tsv = STEP2 / "pairs_real.tsv"
    pairs = list(csv.DictReader(pairs_tsv.open(), delimiter="\t"))
    print(f"[load] pairs: {len(pairs)}")

    cache = _Cache()

    pairs_fa = (OUT / "all_pairs.fna").open("w")
    units_fa = (OUT / "all_pair_units.fna").open("w")
    index_rows: list[dict] = []
    bucket_counts = Counter()
    n_oob = 0

    for i, r in enumerate(pairs, start=1):
        contig = r["contig"]
        meta = contig_lookup.get(contig)
        if meta is None:
            n_oob += 1
            continue
        guid = meta["genome_uid"]
        gid = meta["genome_id"]
        seqs = cache.get(guid)
        contig_seq = seqs.get(contig)
        if contig_seq is None:
            n_oob += 1
            continue

        a0 = int(r["anchor_start0"]); a1 = int(r["anchor_end0"])
        b0 = int(r["second_start0"]); b1 = int(r["second_end0"])
        # Normalize so repeat_1 is upstream, repeat_2 is downstream
        if a0 <= b0:
            r1_s, r1_e, r1_q = a0, a1, r["anchor_qseqid"]
            r2_s, r2_e, r2_q = b0, b1, r["second_qseqid"]
            r1_role = "anchor"; r2_role = "second"
            r2_pident = r["second_pident"]; r2_cov = r["second_cov"]
        else:
            r1_s, r1_e, r1_q = b0, b1, r["second_qseqid"]
            r2_s, r2_e, r2_q = a0, a1, r["anchor_qseqid"]
            r1_role = "second"; r2_role = "anchor"
            r2_pident = r["anchor_pident"]; r2_cov = r["anchor_cov"]
        sp_s, sp_e = r1_e, r2_s
        if sp_e < sp_s:
            n_oob += 1
            continue

        pair_start = r1_s
        pair_end = r2_e

        if pair_start < 0 or pair_end > len(contig_seq) or pair_start >= pair_end:
            n_oob += 1
            continue

        pair_id = f"step2_pair_{i:06d}"
        bucket = context_mod.classify_locus(contig, pair_start, pair_end,
                                             hard, soft, contig_lengths)
        bucket_counts[bucket] += 1

        pair_seq = contig_seq[pair_start:pair_end]
        pairs_fa.write(f">{pair_id}\n{pair_seq}\n")

        r1_seq = contig_seq[r1_s:r1_e]
        sp_seq = contig_seq[sp_s:sp_e]
        r2_seq = contig_seq[r2_s:r2_e]
        units_fa.write(f">{pair_id}|repeat_1\n{r1_seq}\n")
        units_fa.write(f">{pair_id}|spacer\n{sp_seq}\n")
        units_fa.write(f">{pair_id}|repeat_2\n{r2_seq}\n")

        has_amb_pair = any(c not in "ACGT" for c in pair_seq)
        row = {
            "pair_id": pair_id,
            "genome_id": gid,
            "contig_uid": contig,
            "step2_strand": r["strand"],
            "pair_start0": pair_start,
            "pair_end0": pair_end,
            "pair_length": pair_end - pair_start,
            "pair_sequence_sha256": hashlib.sha256(pair_seq.encode()).hexdigest(),
            "has_ambiguous_bases_pair": "true" if has_amb_pair else "false",
            "step2_context_bucket": bucket,
            "repeat_1_role_in_pair": r1_role,
            "repeat_1_query_id": r1_q,
            "repeat_1_start0": r1_s, "repeat_1_end0": r1_e,
            "repeat_1_length": r1_e - r1_s,
            "spacer_start0": sp_s, "spacer_end0": sp_e,
            "spacer_length": sp_e - sp_s,
            "repeat_2_role_in_pair": r2_role,
            "repeat_2_query_id": r2_q,
            "repeat_2_start0": r2_s, "repeat_2_end0": r2_e,
            "repeat_2_length": r2_e - r2_s,
            "repeat_2_pident_vs_query": r2_pident,
            "repeat_2_cov_vs_query": r2_cov,
            "anchor_pident_vs_query": r["anchor_pident"],
            "anchor_cov_vs_query": r["anchor_cov"],
            "step2_second_is_full": r.get("second_is_full", ""),
        }
        index_rows.append(row)

        if i % 2000 == 0:
            print(f"[extract] {i}/{len(pairs)}")

    pairs_fa.close(); units_fa.close()
    print(f"[extract] done; wrote {len(index_rows)} pairs, skipped {n_oob} OOB/missing")

    idx_cols = list(index_rows[0].keys()) if index_rows else []
    with (OUT / "all_pairs_index.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=idx_cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in index_rows:
            w.writerow(r)

    # README
    readme = []
    readme.append("# Step 2 pair-level sequence export (all 15,472 pairs)\n")
    readme.append("Pure extraction. Every pair in `step2_pilot/pairs_real.tsv` "
                   "has its full span + three component sequences dumped from "
                   "the staged FASTA (genomic-forward DNA). Pair-level "
                   "correspondence is preserved: identical component "
                   "sequences across pairs are NOT deduped; each pair keeps "
                   "its own records.\n")
    readme.append("## Query library — origin confirmation\n")
    readme.append("All anchors / seconds come from the 213 NR100 HQ repeat "
                   "sequences in "
                   "`results/step1_patched_min_score_0.75/repeats.hq.nr100.fna`. "
                   "That library was built from the 652 Bona-fide + "
                   "coordinate-pass identify arrays on the same frozen "
                   "1000-genome batch, deduped (exact + RC). No external "
                   "repeat library is used.\n")
    readme.append("## Record counts\n")
    readme.append(f"- Pairs exported: **{len(index_rows)}** (of 15,472 raw "
                   f"pairs; {n_oob} dropped because contig/sequence lookup "
                   f"failed)")
    readme.append(f"- Full-span FASTA: `all_pairs.fna`  "
                   f"({len(index_rows)} records)")
    readme.append(f"- Unit FASTA: `all_pair_units.fna` "
                   f"({len(index_rows) * 3} records = pairs × 3 unit types)")
    readme.append(f"- Index: `all_pairs_index.tsv` "
                   f"({len(index_rows)} rows × {len(idx_cols)} columns)")
    readme.append("")
    readme.append("## Step 2 context bucket distribution\n")
    readme.append("| bucket | n_pairs |")
    readme.append("|---|---:|")
    for k, n in sorted(bucket_counts.items(), key=lambda kv: -kv[1]):
        readme.append(f"| `{k}` | {n} |")
    readme.append("")
    readme.append("## FASTA ID conventions\n")
    readme.append("- `all_pairs.fna`: `step2_pair_<NNNNNN>`")
    readme.append("- `all_pair_units.fna`: `step2_pair_<NNNNNN>|{repeat_1,"
                   "spacer,repeat_2}`")
    readme.append("- `repeat_1` is always the UPSTREAM repeat of the pair "
                   "(lower start0). The `repeat_*_role_in_pair` column in the "
                   "index records which was the original BLAST anchor vs "
                   "secondary.")
    readme.append("")
    readme.append("## Non-negotiables\n")
    readme.append("- Only source of sequence: `input/staged/<genome_uid>.fna`.")
    readme.append("- All sequences genomic-forward; no RC applied.")
    readme.append("- Ambiguous bases preserved as-is.")
    readme.append("- No deduplication (pair-level correspondence intact).")
    readme.append("- No new search, no boundary / threshold / parameter change.")
    readme.append("- `allow_full_run: false` preserved. sbatch-only.")
    (OUT / "README.md").write_text("\n".join(readme) + "\n")

    dt = time.time() - t0
    print(f"[done] {dt:.1f}s")
    print(f"wrote: all_pairs.fna, all_pair_units.fna, "
          f"all_pairs_index.tsv, README.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
