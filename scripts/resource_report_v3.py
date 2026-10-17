#!/usr/bin/env python3
"""Phase J1: clarify the per-tool vs combined wall time story + record
that SLURM did not capture actual TotalCPU (process-level CPU time)."""
from __future__ import annotations

import csv
import json
import subprocess
import statistics
from pathlib import Path

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
OUT = PROJECT / "results" / "batch1000_delivery_metrics"


def _elapsed_s(s: str) -> int:
    if "-" in s:
        d, rest = s.split("-"); parts = rest.split(":")
        return int(d) * 86400 + int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
    parts = s.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
    return 0


def _pull(name: str) -> list[dict]:
    r = subprocess.run(
        ["sacct", "-X", "-u", "kuangh", "-S", "2026-10-15",
         "--format=JobID,JobName,State,Elapsed,CPUTimeRAW,TotalCPU,AllocCPUS",
         "-P", "-n"],
        capture_output=True, text=True,
    )
    out = []
    for ln in r.stdout.strip().splitlines():
        p = ln.split("|")
        if len(p) < 7 or p[1].strip() != name or p[2] != "COMPLETED":
            continue
        try:
            out.append({
                "jobid": p[0],
                "wall_s": _elapsed_s(p[3]),
                "cputime_raw_s": int(p[4]),
                "totalcpu_s": _elapsed_s(p[5]),
                "alloc_cpus": int(p[6]),
            })
        except ValueError:
            continue
    return out


def _mean(vs): return round(statistics.mean(vs), 2) if vs else 0


def main() -> int:
    import sys as _s
    _s.path.insert(0, "/home/kuangh/tools/CRISPR_array_related/scripts")
    from _slurm_gate import require_slurm
    require_slurm()
    identify = _pull("cri_id_patch")
    cctyper = _pull("cctyper")

    totalcpu_captured_identify = sum(1 for r in identify if r["totalcpu_s"] > 0)
    totalcpu_captured_cctyper = sum(1 for r in cctyper if r["totalcpu_s"] > 0)

    md = []
    md.append("# Resource accounting (v3 — pool clarification + TotalCPU note)\n")
    md.append("## What the 'wall @ 30 concurrent' numbers mean\n")
    md.append("Each tool's extrapolated wall time @ 30 concurrent 4-CPU jobs is "
               "computed as `mean_wall_per_job × full_db_N / 30`. This gives:\n")
    id_mean_wall = _mean([r["wall_s"] for r in identify])
    cc_mean_wall = _mean([r["wall_s"] for r in cctyper])
    full_db = 107_235
    md.append(f"- identify alone: `{id_mean_wall}s × {full_db} / 30 / 3600` = "
               f"**{round(id_mean_wall * full_db / 30 / 3600, 1)} h**")
    md.append(f"- cctyper alone: `{cc_mean_wall}s × {full_db} / 30 / 3600` = "
               f"**{round(cc_mean_wall * full_db / 30 / 3600, 1)} h**")
    md.append("")
    md.append("The two pools run in **different conda envs** so they can share "
               "the SLURM scheduler but NOT a single 30-slot pool. The combined "
               "wall-time outcome depends on how you allocate the 30 slots:\n")
    md.append("| scenario | wall-hours (combined) |")
    md.append("|---|---:|")
    md.append(f"| Sequential: 30 slots do identify first then cctyper | "
               f"**{round((id_mean_wall + cc_mean_wall) * full_db / 30 / 3600, 1)}** |")
    md.append(f"| Two separate 30-slot pools (identify + cctyper in parallel) | "
               f"**{round(max(id_mean_wall, cc_mean_wall) * full_db / 30 / 3600, 1)}** "
               f"(60 concurrent CPUs × 4 = 240 cores at peak) |")
    md.append(f"| One 30-slot pool, 50/50 split | "
               f"**{round((id_mean_wall + cc_mean_wall) * full_db / 15 / 3600, 1)}** |")
    md.append("")
    md.append("The 'overlapping_pipelines' number in the v2 report (129.3 h) was "
               "the second row above — two separate 30-slot pools, which doubles "
               "the actual CPU footprint. If a single 30-slot cap is enforced, "
               "the real wall time is the first row (~154 h).\n")

    md.append("## Allocated CPU vs actual TotalCPU\n")
    md.append("SLURM reported `TotalCPU=00:00:00` for every one of the "
               f"{len(identify) + len(cctyper)} completed jobs "
               "(identify + cctyper) — i.e. process-level CPU time was not "
               "captured for this cluster/qos. The only CPU metric we can "
               "report accurately is **AllocCPUS × Elapsed = CPUTimeRAW**, "
               "i.e. the *reserved* CPU-time, not the actually-consumed CPU time.\n")
    md.append(f"- identify jobs with non-zero TotalCPU: "
               f"**{totalcpu_captured_identify} / {len(identify)}**")
    md.append(f"- cctyper jobs with non-zero TotalCPU: "
               f"**{totalcpu_captured_cctyper} / {len(cctyper)}**")
    md.append("")
    md.append("The real process-level CPU is almost certainly **lower** than "
               "the allocated figure (external tools and prodigal are not perfectly "
               "parallel). The 18,533 cpu-hour extrapolation is an upper bound. "
               "To get a true TotalCPU, re-run would be required with cgroup "
               "accounting enabled — explicitly out-of-scope for this round.\n")

    (OUT / "SUMMARY_v3.md").write_text("\n".join(md) + "\n")
    (OUT / "e_runtime_disk_v3.json").write_text(json.dumps({
        "identify_mean_wall_s": id_mean_wall,
        "cctyper_mean_wall_s": cc_mean_wall,
        "wall_hours_at_30_concurrent": {
            "sequential_single_pool": round(
                (id_mean_wall + cc_mean_wall) * full_db / 30 / 3600, 1),
            "two_parallel_pools_30_each": round(
                max(id_mean_wall, cc_mean_wall) * full_db / 30 / 3600, 1),
            "single_pool_split_15_15": round(
                (id_mean_wall + cc_mean_wall) * full_db / 15 / 3600, 1),
        },
        "totalcpu_capture_rate": {
            "identify": f"{totalcpu_captured_identify}/{len(identify)}",
            "cctyper": f"{totalcpu_captured_cctyper}/{len(cctyper)}",
        },
        "note": ("CPUTimeRAW = AllocCPUS × Elapsed = *allocated* CPU-sec. "
                 "TotalCPU = *actually consumed* CPU-sec. Only the former "
                 "is non-zero in this cluster's sacct output; the real CPU "
                 "footprint is almost certainly smaller than the extrapolation."),
    }, indent=2) + "\n")
    print(f"wrote {OUT}/SUMMARY_v3.md and e_runtime_disk_v3.json")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
