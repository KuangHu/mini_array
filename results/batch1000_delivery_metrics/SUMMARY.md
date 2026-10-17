# Batch-1000 delivery metrics

## (a) Library sizes

| level | primary (0.75) | strict (0.90) |
|---|---:|---:|
| candidates_total | 69278 | 69278 |
| bona_fide_coord_pass | 652 | 652 |
| hq_arrays | 239 | 75 |
| reserve_arrays | 1318 | 1482 |
| context_only_arrays | 67621 | 67621 |
| invalid_arrays | 100 | 100 |
| nr100_repeats_hq | 213 | 65 |
| families_hq | 155 | 39 |

## (b) Filter loss bona-fide → HQ (primary 0.75)

- bona-fide calls: **658**
- coord-check pass: **652**
- HQ at primary: **239**
- reserve at primary: **413**

| gate | count |
|---|---:|
| `flank` | 401 |
| `ambiguous_bases_in_context` | 22 |
| `n_full_repeats` | 1 |

## (c) Family cross-genome support (primary, 155 families)

Top 10 families by #genomes supporting:

| family_id | n_repeats | n_arrays | n_genomes | n_phyla |
|---|---:|---:|---:|---:|
| `fam_7a805f541ac1b6ee72cf` | 7 | 9 | 9 | 9 |
| `fam_6f393e61ef58ac27c10c` | 7 | 11 | 7 | 5 |
| `fam_32bff8536f7e4b0e82bd` | 5 | 5 | 5 | 5 |
| `fam_18747d6a513562247094` | 5 | 8 | 4 | 3 |
| `fam_c746fdff9d0f000b08a8` | 3 | 4 | 4 | 3 |
| `fam_88b45a39a9cd4a86b359` | 3 | 4 | 3 | 3 |
| `fam_0cf6ecb62bb5aafb0abe` | 2 | 2 | 2 | 2 |
| `fam_12de37b63356e8a5bf00` | 2 | 2 | 2 | 2 |
| `fam_1b23be31875941438c5a` | 2 | 2 | 2 | 2 |
| `fam_26dc06856c7f90ece5a6` | 2 | 3 | 2 | 1 |

## (d) Stratified recovery

### by source_type

| key | n_genomes | candidates | bona-fide | hq |
|---|---:|---:|---:|---:|
| `MAG` | 933 | 64681 | 570 | 165 |
| `isolate` | 67 | 4597 | 82 | 74 |

### by assembly contiguity

| key | n_genomes | candidates | bona-fide | hq |
|---|---:|---:|---:|---:|
| `101-500 contigs` | 567 | 37925 | 337 | 88 |
| `11-100 contigs` | 285 | 17053 | 203 | 82 |
| `>500 contigs` | 85 | 8299 | 48 | 12 |
| `1-10 contigs` | 63 | 6001 | 64 | 57 |

### by phylum (top 25 by #genomes)

| phylum | n_genomes | candidates | bona-fide | hq |
|---|---:|---:|---:|---:|
| p__4484-113 | 7 | 571 | 1 | 1 |
| p__AUK180 | 7 | 601 | 1 | 1 |
| p__Acidobacteriota | 7 | 402 | 2 | 1 |
| p__Actinomycetota | 7 | 1146 | 2 | 2 |
| p__Aerophobota | 7 | 280 | 2 | 0 |
| p__Aquificota | 7 | 86 | 4 | 1 |
| p__Armatimonadota | 7 | 453 | 5 | 1 |
| p__Atribacterota | 7 | 230 | 1 | 0 |
| p__BMS3Abin14 | 7 | 599 | 0 | 0 |
| p__Bacillota | 7 | 203 | 2 | 2 |
| p__Bacillota_A | 7 | 660 | 2 | 2 |
| p__Bacillota_B | 7 | 493 | 16 | 15 |
| p__Bacillota_C | 7 | 568 | 4 | 3 |
| p__Bacillota_D | 7 | 312 | 6 | 2 |
| p__Bacillota_E | 7 | 250 | 4 | 2 |
| p__Bacillota_F | 7 | 205 | 8 | 8 |
| p__Bacillota_G | 7 | 127 | 7 | 1 |
| p__Bacillota_H | 7 | 140 | 5 | 3 |
| p__Bacillota_I | 7 | 84 | 0 | 0 |
| p__Bacteroidota | 7 | 516 | 0 | 0 |
| p__Bdellovibrionota | 7 | 219 | 0 | 0 |
| p__Bipolaricaulota | 7 | 102 | 2 | 0 |
| p__CAIJMQ01 | 7 | 131 | 0 | 0 |
| p__CG03 | 7 | 200 | 3 | 1 |
| p__CLD3 | 7 | 276 | 8 | 1 |

## (e) Runtime and disk

| tool | n_total | n_completed | n_failed | wall_median_s | wall_mean_s | wall_max_s | disk_median_bytes |
|---|---:|---:|---:|---:|---:|---:|---:|
| identify | 1126 | 1126 | 0 | 96.0 | 130.3 | 1209 | 194494 |
| cctyper | 1310 | 993 | 317 | 23 | 25.3 | 84 | 138821 |

- staged input total: 2.90 GB

## (f) Full-DB (107k) resource extrapolation

- CPU-hours est: **14178.8**
- Wall-hours @ 30 concurrent jobs: **118.2**
- identify output est: 20.86 GB
- cctyper output est: 14.89 GB
- staged input est: 310.95 GB

