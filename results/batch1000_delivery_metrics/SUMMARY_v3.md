# Resource accounting (v3 — pool clarification + TotalCPU note)

## What the 'wall @ 30 concurrent' numbers mean

Each tool's extrapolated wall time @ 30 concurrent 4-CPU jobs is computed as `mean_wall_per_job × full_db_N / 30`. This gives:

- identify alone: `130.26s × 107235 / 30 / 3600` = **129.3 h**
- cctyper alone: `25.28s × 107235 / 30 / 3600` = **25.1 h**

The two pools run in **different conda envs** so they can share the SLURM scheduler but NOT a single 30-slot pool. The combined wall-time outcome depends on how you allocate the 30 slots:

| scenario | wall-hours (combined) |
|---|---:|
| Sequential: 30 slots do identify first then cctyper | **154.4** |
| Two separate 30-slot pools (identify + cctyper in parallel) | **129.3** (60 concurrent CPUs × 4 = 240 cores at peak) |
| One 30-slot pool, 50/50 split | **308.9** |

The 'overlapping_pipelines' number in the v2 report (129.3 h) was the second row above — two separate 30-slot pools, which doubles the actual CPU footprint. If a single 30-slot cap is enforced, the real wall time is the first row (~154 h).

## Allocated CPU vs actual TotalCPU

SLURM reported `TotalCPU=00:00:00` for every one of the 2119 completed jobs (identify + cctyper) — i.e. process-level CPU time was not captured for this cluster/qos. The only CPU metric we can report accurately is **AllocCPUS × Elapsed = CPUTimeRAW**, i.e. the *reserved* CPU-time, not the actually-consumed CPU time.

- identify jobs with non-zero TotalCPU: **0 / 1126**
- cctyper jobs with non-zero TotalCPU: **0 / 993**

The real process-level CPU is almost certainly **lower** than the allocated figure (external tools and prodigal are not perfectly parallel). The 18,533 cpu-hour extrapolation is an upper bound. To get a true TotalCPU, re-run would be required with cgroup accounting enabled — explicitly out-of-scope for this round.

