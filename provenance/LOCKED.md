# Step 1 — frozen implementation spec

**Date locked:** 2026-10-15T21:33:13-07:00
**Approval ref:** User review 2026-10-15, "验收决定：按当前报告，受控 patched Step 1 pilot 通过".

Downstream work must use these exact artifacts. Do NOT change classifier
model files, feature definitions, search parameters, or thresholds without
re-approval.

## Pinned software

| Component | Source | Commit | Notes |
|---|---|---|---|
| CRISPRidentify (upstream) | https://github.com/BackofenLab/CRISPRidentify.git | `474896ce3c9b89a609e69ac9d2fcdf78e4e2aad4` | reference commit; use CRISPRidentify.patched/ in practice |
| CRISPRidentify (patched fork) | software/CRISPRidentify.patched/ | patch set below | always invoked instead of upstream |
| CRISPRCasTyper | https://github.com/Russel88/CRISPRCasTyper.git | `a2d51831d4802533bb113bac23cf9df2853063b0` | version 1.9.0 |
| CRISPRloci (resource source) | https://github.com/BackofenLab/CRISPRloci.git | `d25287b4f8a1c7dc1f885aeda25e61f85b1aecdc` | sole source for models_tandem.hmm |

## Pinned resources

| Resource | SHA256 | Size |
|---|---|---|
| models_tandem.hmm | `818363ec57f0b3b664ce4267ff15314090046a8b17fd709ac7b5926a15a67336` | 13723938 bytes |
| Verified_repeats_dataset1.fa (md5) | `5d2677c20ebdee53f7945eefb5001738` | identical to CRISPRloci copy |
| Verified_repeats_dataset2.fa (md5) | `e936ddec90de10911a14fc8e53c58ec2` | identical to CRISPRloci copy |

## Pinned patches

| File | SHA256 |
|---|---|
| `provenance/patches/components_evaluation.py.patch` | `655f8999bbd4bae1d569f2d295c976313cc02c24f4353622eab703b56c34c246` |
| `provenance/patches/CRISPRidentify.py.patch` | `bb9904e78d90ebf016b502620a67a096ca8a9970b57584f33e250b5299e4afb7` |
| `provenance/patches/_resources.py.new` | `b7061a66efb5e1e5a16fe581b49ff18b4722d01e15769ca3aa77942ac2938b78` |

## Pinned environments

| Env | Explicit lockfile | SHA256 |
|---|---|---|
| envs/identify (Py 3.7.6, TF 2.3, sklearn 0.22) | provenance/identify.conda.explicit.txt | `8094fefa714b45a4adaeec6463a144ebe7c7eb8702e7e9377b3df3800652ef13` |
| envs/driver (Py 3.11, cctyper 1.9.0, edlib, biopython) | provenance/driver.conda.explicit.txt | `62167c08e87e8765f03b64311219b29606db8d425bca0568f16e074dd9bb18cb` |

## Pinned classifier models

| Dir | Hash of file list |
|---|---|
| software/CRISPRidentify.patched/trained_models/ | `b693f9d2f7829a152df9ea9c4fbbe46ab631e8242b9c1c5d8a134066d7ad4d18` |

## Operating QC configuration

| Parameter | Primary (0.75 standard) | Strict comparator (0.90) |
|---|---|---|
| seed_qc.min_score | 0.75 | 0.90 |
| seed_qc.min_full_repeats | 3 | 3 |
| seed_qc.min_unique_spacers | 2 | 2 |
| seed_qc.context_bp | 100 | 100 |
| seed_qc.max_ambiguous_bases_in_context | 0 | 0 |
| seed_qc.min_copy_normalized_edit_similarity | 0.90 | 0.90 |
| seed_qc.min_copy_length_ratio | 0.90 | 0.90 |
| seed_qc.require_cas | false | false |
| cas_annotation.max_local_gap_bp | 10000 | 10000 |
| families.min_normalized_edit_similarity | 0.90 | 0.90 |
| families.min_length_ratio | 0.90 | 0.90 |

**Non-negotiables:**
- Native primary `Bona-fide` category required for HQ; Possible stays in reserve.
- Flank: 100 bp required on both sides, measured against actual contig length.
- Cas absence is NOT a disqualifying condition.
- CRISPRCasTyper disagreement is NOT a disqualifying condition.
- Coordinate check (per step1.md §4) must pass; failures stay `invalid`.
