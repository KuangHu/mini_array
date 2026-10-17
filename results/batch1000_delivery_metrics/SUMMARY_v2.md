# Resource accounting (v2 — correct units)

## (e) Runtime per tool (sacct COMPLETED jobs)

| tool | n | wall_mean_s | wall_median_s | wall_p95_s | wall_max_s | cpu_mean_cpu·s | cpu_median_cpu·s | cpu_h_mean/job |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| identify | 1126 | 130.26 | 96.0 | 343 | 1209 | 521.04 | 384.0 | 0.1447 |
| cctyper | 993 | 25.28 | 23 | 45 | 84 | 101.12 | 92 | 0.0281 |

- Each job allocates **4 CPUs**; `cpu_sec = alloc_cpus × wall_sec` (sacct `CPUTimeRAW`).
- **Per-genome combined cpu-hours (mean)**: 0.1728
- **Per-genome combined cpu-hours (median)**: 0.1322
- Per-genome combined wall-sec (mean): 155.54

## Disk accounting (split)

### Fixed overhead (does not scale with input)
| component | bytes | GB |
|---|---:|---:|
| envs | 5508699838 | 5.51 |
| software | 982948577 | 0.98 |
| provenance | 66112 | 0.00 |
| **sum** | 6491714527 | **6.49** |

### Variable per-genome (scales linearly)
| component | total_bytes (1000 g) | bytes/genome | MB/genome |
|---|---:|---:|---:|
| staged_input | 2899725668 | 2899725 | 2.90 |
| raw_patched_identify | 293523608 | 293523 | 0.29 |
| raw_cctyper | 154064924 | 154064 | 0.15 |
| **sum** | — | 3347314 | **3.35** |

### Peak temporary during run
- `work_patched/` sampled post-cleanup: 141533175 B (not representative of mid-run peak).

## (f) Full-DB extrapolation (107,235 bac120 reps)

### CPU time (uses MEAN per-genome; medians underestimate the long tail)
| tool | cpu-hours |
|---|---:|
| identify | 15520.5 |
| cctyper | 3012.1 |
| combined | 18532.6 |

### Wall time @ 30 concurrent 4-CPU jobs
| scenario | wall-hours |
|---|---:|
| identify | 129.3 |
| cctyper | 25.1 |
| overlapping_pipelines | 129.3 |

Note: Assumes 30 jobs run concurrently throughout. Each job uses 4 CPUs. identify + cctyper can overlap because they use different envs; the longer tool dominates wall time.

### Disk @ full-DB (steady state after cleanup)
| component | bytes | GB |
|---|---:|---:|
| fixed_overhead | 6491714527 | 6.49 |
| variable_per_genome_times_fullN | 358949238237.0 | 358.95 |
| total_steady_state | 365440952764.0 | 365.44 |

**Peak temporary disk during run depends on concurrency**; measured here only at post-cleanup steady state. Do NOT quote the full-DB total as a validated peak.

## Reconciliation with prior report
- Prior report wrote '476 s × 4 cores = 0.53 cpu-h' — incorrect units; 476 cpu-sec = 0.132 cpu-h.
- Prior report used medians for extrapolation, which underestimates because of a long tail (identify max was 1209 s vs median 96 s).
- Means give a slightly larger total; both are reported here.

