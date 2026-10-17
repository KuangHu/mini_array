#!/usr/bin/env python3
"""Phase H6: compute full delivery metrics for the batch-1000 Step 1 export.

Produces results/batch1000_delivery_metrics/*.tsv + a Markdown summary covering
the user's point-7 deliverables:

  (a) primary (0.75) vs strict (0.90) subsets — arrays / NR100 / families
  (b) per-gate filter loss bona-fide → HQ
  (c) family cross-genome support counts
  (d) stratification by taxonomy / source_type / assembly contiguity
  (e) per-tool wall time, failure rate, per-genome disk cost
  (f) full-DB (113 k) resource extrapolation
"""
from __future__ import annotations

import csv
import json
import os
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
PRIMARY = PROJECT / "results" / "step1_patched_min_score_0.75"
STRICT = PROJECT / "results" / "step1_patched_min_score_0.90"
MANIFEST = PROJECT / "input" / "genomes.tsv"
STAGED_DIR = PROJECT / "input" / "staged"
RAW_ID = PROJECT / "raw_patched" / "identify"
RAW_CC = PROJECT / "raw" / "cctyper"
GTDB_META = Path("/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/bac120_metadata.tsv")

OUT_DIR = PROJECT / "results" / "batch1000_delivery_metrics"


def _read(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _sanitize(gid: str) -> str:
    return gid.replace(".", "_")


def section_a_libraries() -> dict:
    pr = _read(PRIMARY / "arrays.all.tsv")
    st = _read(STRICT / "arrays.all.tsv")
    pr_hq = [r for r in pr if r["qc_tier"] == "hq"]
    st_hq = [r for r in st if r["qc_tier"] == "hq"]
    pr_nr100 = sum(1 for _ in (PRIMARY / "repeats.hq.nr100.fna").read_text().splitlines() if _.startswith(">"))
    st_nr100 = sum(1 for _ in (STRICT / "repeats.hq.nr100.fna").read_text().splitlines() if _.startswith(">"))
    pr_fam = len(_read(PRIMARY / "repeat_families.tsv"))
    st_fam = len(_read(STRICT / "repeat_families.tsv"))
    return {
        "primary_0.75": {
            "candidates_total": len(pr),
            "bona_fide_coord_pass": sum(1 for r in pr if r["category"] == "bona_fide" and r["coordinate_check_pass"] == "true"),
            "hq_arrays": len(pr_hq),
            "reserve_arrays": sum(1 for r in pr if r["qc_tier"] == "reserve"),
            "context_only_arrays": sum(1 for r in pr if r["qc_tier"] == "context_only"),
            "invalid_arrays": sum(1 for r in pr if r["qc_tier"] == "invalid"),
            "nr100_repeats_hq": pr_nr100,
            "families_hq": pr_fam,
        },
        "strict_0.90": {
            "candidates_total": len(st),
            "bona_fide_coord_pass": sum(1 for r in st if r["category"] == "bona_fide" and r["coordinate_check_pass"] == "true"),
            "hq_arrays": len(st_hq),
            "reserve_arrays": sum(1 for r in st if r["qc_tier"] == "reserve"),
            "context_only_arrays": sum(1 for r in st if r["qc_tier"] == "context_only"),
            "invalid_arrays": sum(1 for r in st if r["qc_tier"] == "invalid"),
            "nr100_repeats_hq": st_nr100,
            "families_hq": st_fam,
        },
    }


def section_b_filter_loss() -> dict:
    """Per-gate loss walking bona-fide → HQ at primary threshold."""
    rows = _read(PRIMARY / "arrays.all.tsv")
    bf = [r for r in rows if r["category"] == "bona_fide"]
    stats = {"bona_fide_total": len(bf)}
    coord_ok = [r for r in bf if r["coordinate_check_pass"] == "true"]
    stats["coord_pass"] = len(coord_ok)
    # QC reasons for the ones that are reserve (not HQ)
    reserve = [r for r in coord_ok if r["qc_tier"] == "reserve"]
    hq = [r for r in coord_ok if r["qc_tier"] == "hq"]
    # Tally first reason per reserved candidate
    reason_counts = Counter()
    for r in reserve:
        try:
            reasons = json.loads(r["qc_reasons"])
        except Exception:
            reasons = []
        if not reasons:
            reason_counts["unknown"] += 1
            continue
        # Normalize labels
        for raw in reasons:
            label = raw.split(" (")[0].split("<")[0].split("=")[0].strip()
            reason_counts[label] += 1
    stats["hq_at_primary"] = len(hq)
    stats["reserve_at_primary"] = len(reserve)
    stats["reserve_reasons"] = dict(reason_counts.most_common())
    return stats


def section_c_family_support() -> dict:
    """Per-family cross-genome support at the primary threshold."""
    members = _read(PRIMARY / "repeat_family_members.tsv")
    sources = _read(PRIMARY / "repeat_sources.tsv")
    rep_to_fam = {m["repeat_id"]: m["family_id"] for m in members}
    fam_rows = defaultdict(lambda: {
        "arrays": set(), "genomes": set(), "repeats": set(), "phyla": set()
    })
    # Map genome_id → phylum from manifest
    manifest = {m["genome_id"]: m for m in _read(MANIFEST)}
    for s in sources:
        rid = s["repeat_id"]
        fid = rep_to_fam.get(rid)
        if fid is None:
            continue
        gid = s["genome_id"]
        fam_rows[fid]["repeats"].add(rid)
        fam_rows[fid]["arrays"].add(s["array_id"])
        fam_rows[fid]["genomes"].add(gid)
        phy = manifest.get(gid, {}).get("independent_group_id", "NA")
        fam_rows[fid]["phyla"].add(phy)
    out = {
        "families_total": len(fam_rows),
        "rows": [
            {
                "family_id": fid,
                "n_repeats": len(info["repeats"]),
                "n_arrays": len(info["arrays"]),
                "n_genomes": len(info["genomes"]),
                "n_phyla": len(info["phyla"]),
            }
            for fid, info in sorted(fam_rows.items())
        ],
    }
    return out


def section_d_stratification() -> dict:
    """Bona-fide / HQ recovery by source_type and phylum."""
    manifest = _read(MANIFEST)
    gid_to_meta = {m["genome_id"]: m for m in manifest}
    rows = _read(PRIMARY / "arrays.all.tsv")

    by_source = defaultdict(lambda: {"genomes": set(), "bf": 0, "hq": 0, "any_candidate": 0})
    by_phylum = defaultdict(lambda: {"genomes": set(), "bf": 0, "hq": 0, "any_candidate": 0})
    for m in manifest:
        by_source[m["source_type"]]["genomes"].add(m["genome_id"])
        by_phylum[m["independent_group_id"]]["genomes"].add(m["genome_id"])
    for r in rows:
        meta = gid_to_meta.get(r["genome_id"])
        if not meta:
            continue
        s = meta["source_type"]
        p = meta["independent_group_id"]
        by_source[s]["any_candidate"] += 1
        by_phylum[p]["any_candidate"] += 1
        if r["category"] == "bona_fide" and r["coordinate_check_pass"] == "true":
            by_source[s]["bf"] += 1
            by_phylum[p]["bf"] += 1
            if r["qc_tier"] == "hq":
                by_source[s]["hq"] += 1
                by_phylum[p]["hq"] += 1

    def flatten(buckets):
        return [
            {
                "key": k,
                "n_genomes_in_manifest": len(v["genomes"]),
                "candidates_total": v["any_candidate"],
                "bona_fide_coord_pass": v["bf"],
                "hq_arrays": v["hq"],
            }
            for k, v in sorted(buckets.items(),
                                key=lambda kv: (-len(kv[1]["genomes"]), kv[0]))
        ]

    # Also split by assembly contiguity using contig count from staged contigs.tsv:
    staged_contigs = PROJECT / "results" / "staging" / "contigs.1000.tsv"
    contigs_per_genome = Counter()
    if staged_contigs.is_file():
        for r in _read(staged_contigs):
            contigs_per_genome[r["genome_id"]] += 1
    bucket_labels = [
        ("1-10 contigs", lambda n: n <= 10),
        ("11-100 contigs", lambda n: 10 < n <= 100),
        ("101-500 contigs", lambda n: 100 < n <= 500),
        (">500 contigs", lambda n: n > 500),
    ]
    by_contigs = defaultdict(lambda: {"genomes": set(), "bf": 0, "hq": 0, "any_candidate": 0})
    for gid, n in contigs_per_genome.items():
        for label, pred in bucket_labels:
            if pred(n):
                by_contigs[label]["genomes"].add(gid)
                break
    for r in rows:
        n = contigs_per_genome.get(r["genome_id"], 0)
        label = next((lab for lab, pred in bucket_labels if pred(n)), "?")
        by_contigs[label]["any_candidate"] += 1
        if r["category"] == "bona_fide" and r["coordinate_check_pass"] == "true":
            by_contigs[label]["bf"] += 1
            if r["qc_tier"] == "hq":
                by_contigs[label]["hq"] += 1

    return {
        "by_source_type": flatten(by_source),
        "by_phylum": flatten(by_phylum)[:25],  # top 25 phyla by #genomes
        "by_contig_count": flatten(by_contigs),
    }


def section_e_runtime_disk() -> dict:
    """Per-tool wall time, failure rate, per-genome disk cost."""
    # sacct for identify + cctyper jobs — pull the last 24 h window
    def _sacct(name: str) -> list[dict]:
        r = subprocess.run(
            ["sacct", "-X", "-u", "kuangh", "-S", "2026-10-15",
             "--format=JobID,JobName,State,ExitCode,Elapsed", "-P", "-n"],
            capture_output=True, text=True,
        )
        rows = []
        for ln in r.stdout.strip().splitlines():
            parts = ln.split("|")
            if len(parts) < 5:
                continue
            jid, jn, st, ec, el = parts[:5]
            if jn.strip() == name:
                rows.append({"jobid": jid, "state": st, "exit": ec, "elapsed": el})
        return rows

    def _secs(el: str) -> int:
        # format: HH:MM:SS or D-HH:MM:SS
        if "-" in el:
            d, rest = el.split("-")
            h, m, s = rest.split(":")
            return int(d) * 86400 + int(h) * 3600 + int(m) * 60 + int(s)
        parts = el.split(":")
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
        return 0

    def _stats(rows):
        if not rows:
            return {"n": 0}
        completed = [r for r in rows if r["state"] == "COMPLETED"]
        failed = [r for r in rows if r["state"] != "COMPLETED"]
        secs = [_secs(r["elapsed"]) for r in completed if r["elapsed"]]
        return {
            "n_total": len(rows),
            "n_completed": len(completed),
            "n_failed": len(failed),
            "wall_sec_median": statistics.median(secs) if secs else 0,
            "wall_sec_mean": round(statistics.mean(secs), 1) if secs else 0,
            "wall_sec_max": max(secs) if secs else 0,
            "wall_sec_p95": round(sorted(secs)[int(0.95 * len(secs))], 1) if len(secs) >= 20 else 0,
        }

    identify_rows = _sacct("cri_id_patch")
    cctyper_rows = _sacct("cctyper")
    identify_stats = _stats(identify_rows)
    cctyper_stats = _stats(cctyper_rows)

    def _du(path: Path) -> int:
        """Return total bytes of a directory (fast version)."""
        total = 0
        for root, _, files in os.walk(path):
            for f in files:
                try:
                    total += os.path.getsize(os.path.join(root, f))
                except OSError:
                    pass
        return total

    def _per_genome_disk(root: Path) -> dict:
        sizes = []
        for d in sorted(root.iterdir()):
            if d.is_dir():
                sizes.append(_du(d))
        if not sizes:
            return {"n": 0}
        return {
            "n": len(sizes),
            "median_bytes": int(statistics.median(sizes)),
            "mean_bytes": int(statistics.mean(sizes)),
            "max_bytes": max(sizes),
            "total_bytes": sum(sizes),
        }

    return {
        "identify": {
            **identify_stats,
            "per_genome_disk": _per_genome_disk(RAW_ID),
        },
        "cctyper": {
            **cctyper_stats,
            "per_genome_disk": _per_genome_disk(RAW_CC),
        },
        "staged_input_bytes": _du(STAGED_DIR),
    }


def section_f_extrapolation(runtime: dict, libraries: dict) -> dict:
    """Extrapolate 113k-genome full-DB estimates from this batch."""
    full_n = 107235  # GTDB r220 bac120 reps with FASTAs
    batch_n = libraries["primary_0.75"]["candidates_total"]
    # Candidates are proportional to #genomes; use raw genomes = 1000.
    scale = full_n / 1000
    id_sec = runtime["identify"].get("wall_sec_median", 0)
    cc_sec = runtime["cctyper"].get("wall_sec_median", 0)
    cpu_h_per_job = (id_sec + cc_sec) * 4 / 3600  # 4 CPUs per job
    disk_id = runtime["identify"]["per_genome_disk"].get("median_bytes", 0)
    disk_cc = runtime["cctyper"]["per_genome_disk"].get("median_bytes", 0)
    staged_bytes = runtime["staged_input_bytes"]
    return {
        "assumed_full_n": full_n,
        "cpu_hours_est": round(cpu_h_per_job * full_n, 1),
        "wall_hours_at_30_concurrent_jobs": round(
            (id_sec + cc_sec) * full_n / 30 / 3600, 1),
        "disk_identify_bytes_est": disk_id * full_n,
        "disk_cctyper_bytes_est": disk_cc * full_n,
        "disk_staged_bytes_est": int(staged_bytes * full_n / 1000),
        "notes": [
            "CPU-hour estimate assumes identify + cctyper medians per genome.",
            "Wall-hour estimate assumes 30 jobs in flight on the standard partition.",
            "Disk estimates scale per-genome medians to the full 107k bac120 reps.",
        ],
    }


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


def _write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                            lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main() -> int:
    import sys as _s
    _s.path.insert(0, "/home/kuangh/tools/CRISPR_array_related/scripts")
    from _slurm_gate import require_slurm
    require_slurm()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    a = section_a_libraries()
    b = section_b_filter_loss()
    c = section_c_family_support()
    d = section_d_stratification()
    e = section_e_runtime_disk()
    f = section_f_extrapolation(e, a)

    _write_json(OUT_DIR / "a_libraries.json", a)
    _write_json(OUT_DIR / "b_filter_loss.json", b)
    _write_json(OUT_DIR / "c_family_support.json", c)
    _write_json(OUT_DIR / "d_stratification.json", d)
    _write_json(OUT_DIR / "e_runtime_disk.json", e)
    _write_json(OUT_DIR / "f_full_db_extrapolation.json", f)
    _write_tsv(OUT_DIR / "c_family_support.tsv",
               ["family_id", "n_repeats", "n_arrays", "n_genomes", "n_phyla"],
               c["rows"])

    # Produce a one-page Markdown summary
    md = ["# Batch-1000 delivery metrics", ""]
    md.append("## (a) Library sizes")
    md.append("")
    md.append("| level | primary (0.75) | strict (0.90) |")
    md.append("|---|---:|---:|")
    for k in ("candidates_total", "bona_fide_coord_pass", "hq_arrays",
              "reserve_arrays", "context_only_arrays", "invalid_arrays",
              "nr100_repeats_hq", "families_hq"):
        md.append(f"| {k} | {a['primary_0.75'][k]} | {a['strict_0.90'][k]} |")
    md.append("")

    md.append("## (b) Filter loss bona-fide → HQ (primary 0.75)")
    md.append("")
    md.append(f"- bona-fide calls: **{b['bona_fide_total']}**")
    md.append(f"- coord-check pass: **{b['coord_pass']}**")
    md.append(f"- HQ at primary: **{b['hq_at_primary']}**")
    md.append(f"- reserve at primary: **{b['reserve_at_primary']}**")
    if b["reserve_reasons"]:
        md.append("")
        md.append("| gate | count |")
        md.append("|---|---:|")
        for k, v in b["reserve_reasons"].items():
            md.append(f"| `{k}` | {v} |")
    md.append("")

    md.append(f"## (c) Family cross-genome support (primary, {c['families_total']} families)")
    md.append("")
    md.append("Top 10 families by #genomes supporting:")
    md.append("")
    md.append("| family_id | n_repeats | n_arrays | n_genomes | n_phyla |")
    md.append("|---|---:|---:|---:|---:|")
    for row in sorted(c["rows"], key=lambda r: -r["n_genomes"])[:10]:
        md.append(f"| `{row['family_id'][:30]}` | {row['n_repeats']} | {row['n_arrays']} | {row['n_genomes']} | {row['n_phyla']} |")
    md.append("")

    md.append("## (d) Stratified recovery")
    md.append("")
    md.append("### by source_type")
    md.append("")
    md.append("| key | n_genomes | candidates | bona-fide | hq |")
    md.append("|---|---:|---:|---:|---:|")
    for r in d["by_source_type"]:
        md.append(f"| `{r['key']}` | {r['n_genomes_in_manifest']} | {r['candidates_total']} | {r['bona_fide_coord_pass']} | {r['hq_arrays']} |")
    md.append("")
    md.append("### by assembly contiguity")
    md.append("")
    md.append("| key | n_genomes | candidates | bona-fide | hq |")
    md.append("|---|---:|---:|---:|---:|")
    for r in d["by_contig_count"]:
        md.append(f"| `{r['key']}` | {r['n_genomes_in_manifest']} | {r['candidates_total']} | {r['bona_fide_coord_pass']} | {r['hq_arrays']} |")
    md.append("")
    md.append("### by phylum (top 25 by #genomes)")
    md.append("")
    md.append("| phylum | n_genomes | candidates | bona-fide | hq |")
    md.append("|---|---:|---:|---:|---:|")
    for r in d["by_phylum"]:
        md.append(f"| {r['key']} | {r['n_genomes_in_manifest']} | {r['candidates_total']} | {r['bona_fide_coord_pass']} | {r['hq_arrays']} |")
    md.append("")

    md.append("## (e) Runtime and disk")
    md.append("")
    md.append(f"| tool | n_total | n_completed | n_failed | wall_median_s | wall_mean_s | wall_max_s | disk_median_bytes |")
    md.append(f"|---|---:|---:|---:|---:|---:|---:|---:|")
    for name in ("identify", "cctyper"):
        s = e[name]
        md.append(f"| {name} | {s.get('n_total',0)} | {s.get('n_completed',0)} | "
                   f"{s.get('n_failed',0)} | {s.get('wall_sec_median',0)} | "
                   f"{s.get('wall_sec_mean',0)} | {s.get('wall_sec_max',0)} | "
                   f"{s['per_genome_disk'].get('median_bytes',0)} |")
    md.append("")
    md.append(f"- staged input total: {e['staged_input_bytes'] / 1e9:.2f} GB")
    md.append("")

    md.append("## (f) Full-DB (107k) resource extrapolation")
    md.append("")
    md.append(f"- CPU-hours est: **{f['cpu_hours_est']:.1f}**")
    md.append(f"- Wall-hours @ 30 concurrent jobs: **{f['wall_hours_at_30_concurrent_jobs']:.1f}**")
    md.append(f"- identify output est: {f['disk_identify_bytes_est'] / 1e9:.2f} GB")
    md.append(f"- cctyper output est: {f['disk_cctyper_bytes_est'] / 1e9:.2f} GB")
    md.append(f"- staged input est: {f['disk_staged_bytes_est'] / 1e9:.2f} GB")
    md.append("")

    (OUT_DIR / "SUMMARY.md").write_text("\n".join(md) + "\n")
    print(f"wrote {OUT_DIR}/{{a..f}}_*.json + c_family_support.tsv + SUMMARY.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
