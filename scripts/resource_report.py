#!/usr/bin/env python3
"""Phase I1: rebuild batch-1000 resource accounting with explicit units.

Reads sacct for identify + cctyper jobs, computes:
  - wall time (per-tool total + mean + median + p95 + max)
  - CPU time (AllocCPUS × wall), same stats
  - combined (identify + cctyper) per-genome aggregates
  - full-DB (107,235 bac120 reps) extrapolation using means (long tail safe)
  - wall-time extrapolation at explicit concurrency assumptions

Disk accounting split:
  - fixed overhead (envs, software, DBs — doesn't scale with input)
  - variable per-genome (staged + raw outputs, scales linearly)
  - peak temporary during run (work_patched/, best-effort estimate)

Writes results/batch1000_delivery_metrics/e_runtime_disk_v2.json
      results/batch1000_delivery_metrics/f_full_db_extrapolation_v2.json
      results/batch1000_delivery_metrics/SUMMARY_v2.md (sections e + f only)
"""
from __future__ import annotations

import csv
import json
import os
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
OUT = PROJECT / "results" / "batch1000_delivery_metrics"
FULL_DB_N = 107_235


def _elapsed_secs(elapsed: str) -> int:
    if "-" in elapsed:
        d, rest = elapsed.split("-")
        parts = rest.split(":")
        return int(d) * 86400 + int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
    parts = elapsed.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
    return 0


def _sacct_rows(name: str) -> list[dict]:
    r = subprocess.run(
        ["sacct", "-X", "-u", os.environ["USER"], "-S", "2026-10-15",
         "--format=JobID,JobName,State,Elapsed,CPUTimeRAW,AllocCPUS",
         "-P", "-n"],
        capture_output=True, text=True,
    )
    out = []
    for ln in r.stdout.strip().splitlines():
        p = ln.split("|")
        if len(p) < 6:
            continue
        jid, jn, st, el, cpu_raw, alloc = p[:6]
        if jn.strip() != name:
            continue
        if st != "COMPLETED":
            continue
        try:
            wall = _elapsed_secs(el)
            cpu = int(cpu_raw)
            ac = int(alloc)
        except ValueError:
            continue
        out.append({"jobid": jid, "wall_sec": wall, "cpu_sec": cpu, "alloc_cpus": ac})
    return out


def _stats(vals: list[float]) -> dict:
    if not vals:
        return {"n": 0}
    vs = sorted(vals)
    def p(q):
        i = min(len(vs) - 1, int(round((len(vs) - 1) * q)))
        return vs[i]
    return {
        "n": len(vals),
        "sum": sum(vs),
        "mean": round(statistics.mean(vs), 2),
        "median": statistics.median(vs),
        "p95": p(0.95),
        "max": max(vs),
    }


def _du_bytes(path: Path) -> int:
    total = 0
    try:
        for root, _, files in os.walk(path):
            for f in files:
                try:
                    total += os.path.getsize(os.path.join(root, f))
                except OSError:
                    pass
    except OSError:
        pass
    return total


def main() -> int:
    import sys as _s
    _s.path.insert(0, "/home/kuangh/tools/CRISPR_array_related/scripts")
    from _slurm_gate import require_slurm
    require_slurm()
    # ---- (e) runtime ----
    id_rows = _sacct_rows("cri_id_patch")
    cc_rows = _sacct_rows("cctyper")

    def _tool_block(rows):
        wall = [r["wall_sec"] for r in rows]
        cpu = [r["cpu_sec"] for r in rows]  # AllocCPUS × wall for each job
        return {
            "n_completed": len(rows),
            "wall_sec": _stats(wall),
            "cpu_sec": _stats(cpu),
            "alloc_cpus_per_job": sorted({r["alloc_cpus"] for r in rows}),
        }

    identify = _tool_block(id_rows)
    cctyper = _tool_block(cc_rows)

    # Combined per-genome: mean(identify_cpu) + mean(cctyper_cpu)
    mean_cpu_per_genome_sec = (identify["cpu_sec"].get("mean", 0)
                                + cctyper["cpu_sec"].get("mean", 0))
    mean_wall_per_genome_sec = (identify["wall_sec"].get("mean", 0)
                                 + cctyper["wall_sec"].get("mean", 0))
    median_cpu_per_genome_sec = (identify["cpu_sec"].get("median", 0)
                                  + cctyper["cpu_sec"].get("median", 0))
    median_wall_per_genome_sec = (identify["wall_sec"].get("median", 0)
                                   + cctyper["wall_sec"].get("median", 0))

    runtime = {
        "window": "sacct from 2026-10-15 00:00, filter --state=COMPLETED",
        "note": "CPUTimeRAW is reported by SLURM; it equals AllocCPUS × wall for each job.",
        "identify": identify,
        "cctyper": cctyper,
        "combined_per_genome": {
            "cpu_sec_mean": round(mean_cpu_per_genome_sec, 2),
            "cpu_sec_median": median_cpu_per_genome_sec,
            "cpu_hours_mean": round(mean_cpu_per_genome_sec / 3600, 4),
            "cpu_hours_median": round(median_cpu_per_genome_sec / 3600, 4),
            "wall_sec_mean": round(mean_wall_per_genome_sec, 2),
            "wall_sec_median": median_wall_per_genome_sec,
        },
    }

    # ---- disk accounting split ----
    def _pg(path: Path, n: int) -> dict:
        total = _du_bytes(path)
        return {"total_bytes": total, "n_items": n,
                "per_item_bytes": round(total / n, 1) if n else 0}

    staged_dir = PROJECT / "input" / "staged"
    raw_id = PROJECT / "raw_patched" / "identify"
    raw_cc = PROJECT / "raw" / "cctyper"
    work_patched = PROJECT / "work_patched"
    results_dir = PROJECT / "results"
    envs_dir = PROJECT / "envs"
    sw_dir = PROJECT / "software"
    prov_dir = PROJECT / "provenance"

    disk = {
        "fixed_overhead": {
            "envs": {"path": str(envs_dir), "bytes": _du_bytes(envs_dir)},
            "software": {"path": str(sw_dir), "bytes": _du_bytes(sw_dir)},
            "provenance": {"path": str(prov_dir), "bytes": _du_bytes(prov_dir)},
        },
        "variable_per_genome": {
            "staged_input": _pg(staged_dir, 1000),
            "raw_patched_identify": _pg(raw_id, 1000),
            "raw_cctyper": _pg(raw_cc, 1000),
        },
        "peak_temporary_during_run": {
            "work_patched": {"path": str(work_patched), "bytes": _du_bytes(work_patched),
                              "note": "steady-state after cleanup — true peak was higher"},
        },
        "results_cumulative": {"path": str(results_dir), "bytes": _du_bytes(results_dir),
                                "note": "4 variant trees + reports + metrics"},
    }
    fixed_total = sum(v["bytes"] for v in disk["fixed_overhead"].values())
    variable_per_genome = sum(v["per_item_bytes"]
                               for v in disk["variable_per_genome"].values())
    disk["summary"] = {
        "fixed_overhead_bytes": fixed_total,
        "variable_per_genome_bytes": round(variable_per_genome, 1),
    }

    # ---- (f) full-DB extrapolation (uses MEAN per-genome, not median) ----
    full_db = FULL_DB_N
    extrap = {
        "assumed_full_db_n_genomes": full_db,
        "cpu_hours_total_at_mean": {
            "identify": round(identify["cpu_sec"].get("mean", 0) * full_db / 3600, 1),
            "cctyper": round(cctyper["cpu_sec"].get("mean", 0) * full_db / 3600, 1),
            "combined": round(mean_cpu_per_genome_sec * full_db / 3600, 1),
        },
        "wall_hours_total_at_30_concurrent": {
            "note": ("Assumes 30 jobs run concurrently throughout. Each job uses "
                     "4 CPUs. identify + cctyper can overlap because they use "
                     "different envs; the longer tool dominates wall time."),
            "identify": round(identify["wall_sec"].get("mean", 0) * full_db / 30 / 3600, 1),
            "cctyper": round(cctyper["wall_sec"].get("mean", 0) * full_db / 30 / 3600, 1),
            "overlapping_pipelines": round(
                max(identify["wall_sec"].get("mean", 0),
                    cctyper["wall_sec"].get("mean", 0)) * full_db / 30 / 3600, 1),
        },
        "disk_bytes_at_full_db": {
            "fixed_overhead": fixed_total,
            "variable_per_genome_times_fullN": round(variable_per_genome * full_db, 0),
            "total_steady_state": round(fixed_total + variable_per_genome * full_db, 0),
            "note": ("Peak temporary during run depends on concurrent-jobs setting; "
                     "work_patched sampled at steady state (post-cleanup) at "
                     f"{disk['peak_temporary_during_run']['work_patched']['bytes']} B "
                     "for 1000-genome batch — not representative of in-flight peak."),
        },
        "sources_of_discrepancy_vs_prior_report": [
            "Prior report wrote '476 s × 4 cores = 0.53 cpu-h' — incorrect units; "
            "476 cpu-sec = 0.132 cpu-h.",
            "Prior report used medians for extrapolation, which underestimates "
            "because of a long tail (identify max was 1209 s vs median 96 s).",
            "Means give a slightly larger total; both are reported here.",
        ],
    }

    # ---- write outputs ----
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "e_runtime_disk_v2.json").write_text(json.dumps(runtime | {"disk": disk}, indent=2) + "\n")
    (OUT / "f_full_db_extrapolation_v2.json").write_text(json.dumps(extrap, indent=2) + "\n")

    md = []
    md.append("# Resource accounting (v2 — correct units)\n")
    md.append("## (e) Runtime per tool (sacct COMPLETED jobs)\n")
    md.append("| tool | n | wall_mean_s | wall_median_s | wall_p95_s | wall_max_s | cpu_mean_cpu·s | cpu_median_cpu·s | cpu_h_mean/job |")
    md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name, b in [("identify", identify), ("cctyper", cctyper)]:
        ws = b["wall_sec"]; cs = b["cpu_sec"]
        md.append(f"| {name} | {b['n_completed']} | {ws.get('mean',0)} | {ws.get('median',0)} | "
                   f"{ws.get('p95',0)} | {ws.get('max',0)} | {cs.get('mean',0)} | {cs.get('median',0)} | "
                   f"{round(cs.get('mean',0)/3600, 4)} |")
    md.append("")
    md.append("- Each job allocates **4 CPUs**; `cpu_sec = alloc_cpus × wall_sec` (sacct `CPUTimeRAW`).")
    md.append(f"- **Per-genome combined cpu-hours (mean)**: {runtime['combined_per_genome']['cpu_hours_mean']}")
    md.append(f"- **Per-genome combined cpu-hours (median)**: {runtime['combined_per_genome']['cpu_hours_median']}")
    md.append(f"- Per-genome combined wall-sec (mean): {runtime['combined_per_genome']['wall_sec_mean']}")
    md.append("")
    md.append("## Disk accounting (split)\n")
    md.append("### Fixed overhead (does not scale with input)")
    md.append("| component | bytes | GB |")
    md.append("|---|---:|---:|")
    for k, v in disk["fixed_overhead"].items():
        md.append(f"| {k} | {v['bytes']} | {v['bytes']/1e9:.2f} |")
    md.append(f"| **sum** | {fixed_total} | **{fixed_total/1e9:.2f}** |")
    md.append("")
    md.append("### Variable per-genome (scales linearly)")
    md.append("| component | total_bytes (1000 g) | bytes/genome | MB/genome |")
    md.append("|---|---:|---:|---:|")
    for k, v in disk["variable_per_genome"].items():
        md.append(f"| {k} | {v['total_bytes']} | {int(v['per_item_bytes'])} | {v['per_item_bytes']/1e6:.2f} |")
    md.append(f"| **sum** | — | {int(variable_per_genome)} | **{variable_per_genome/1e6:.2f}** |")
    md.append("")
    md.append("### Peak temporary during run")
    md.append(f"- `work_patched/` sampled post-cleanup: {disk['peak_temporary_during_run']['work_patched']['bytes']} B "
               "(not representative of mid-run peak).")
    md.append("")

    md.append("## (f) Full-DB extrapolation (107,235 bac120 reps)\n")
    md.append("### CPU time (uses MEAN per-genome; medians underestimate the long tail)")
    md.append("| tool | cpu-hours |")
    md.append("|---|---:|")
    for k, v in extrap["cpu_hours_total_at_mean"].items():
        md.append(f"| {k} | {v} |")
    md.append("")
    md.append("### Wall time @ 30 concurrent 4-CPU jobs")
    md.append("| scenario | wall-hours |")
    md.append("|---|---:|")
    for k, v in extrap["wall_hours_total_at_30_concurrent"].items():
        if isinstance(v, (int, float)):
            md.append(f"| {k} | {v} |")
    md.append(f"\nNote: {extrap['wall_hours_total_at_30_concurrent']['note']}\n")
    md.append("### Disk @ full-DB (steady state after cleanup)")
    md.append("| component | bytes | GB |")
    md.append("|---|---:|---:|")
    for k, v in extrap["disk_bytes_at_full_db"].items():
        if isinstance(v, (int, float)):
            md.append(f"| {k} | {v} | {v/1e9:.2f} |")
    md.append(f"\n**Peak temporary disk during run depends on concurrency**; measured here only at post-cleanup steady state. Do NOT quote the full-DB total as a validated peak.\n")

    md.append("## Reconciliation with prior report")
    for n in extrap["sources_of_discrepancy_vs_prior_report"]:
        md.append(f"- {n}")
    md.append("")

    (OUT / "SUMMARY_v2.md").write_text("\n".join(md) + "\n")
    print(f"wrote {OUT}/e_runtime_disk_v2.json, f_full_db_extrapolation_v2.json, SUMMARY_v2.md")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
