# Reproducing the catalog

This repository ships the frozen coordinate + sequence catalog and the
code that produced it. The underlying raw outputs (staged FASTAs, conda
envs, software clones, cctyper per-genome output, patched identify
per-genome output) are NOT shipped — they are large and derivable.

## Prerequisites

- SLURM cluster with `sbatch` reachable from the submit node.
- Conda / mamba.
- GTDB r220 bac120 representative genomes available locally.

## Steps

1. Install envs from the locked explicit files in `provenance/`:
   ```
   mamba create -p envs/identify --file provenance/identify.conda.explicit.txt
   mamba create -p envs/driver   --file provenance/driver.conda.explicit.txt
   ```

2. Pin the three upstream tools to the recorded commits:
   ```
   git clone --recurse-submodules https://github.com/BackofenLab/CRISPRidentify.git software/CRISPRidentify
   git -C software/CRISPRidentify checkout $(cat provenance/CRISPRidentify.commit)

   git clone https://github.com/Russel88/CRISPRCasTyper.git software/CRISPRCasTyper
   git -C software/CRISPRCasTyper checkout $(cat provenance/CRISPRCasTyper.commit)

   git clone --recurse-submodules https://github.com/BackofenLab/CRISPRloci.git software/CRISPRloci
   git -C software/CRISPRloci checkout $(cat provenance/CRISPRloci.commit)
   ```

3. Create the patched CRISPRidentify fork:
   ```
   rsync -a --exclude='.git' software/CRISPRidentify/ software/CRISPRidentify.patched/
   mkdir -p software/CRISPRidentify.patched/tools/hmm_search
   cp software/CRISPRloci/CRISPRidentify/tools/hmm_search/models_tandem.hmm \
      software/CRISPRidentify.patched/tools/hmm_search/
   patch -p0 -d software/CRISPRidentify.patched/ < provenance/patches/CRISPRidentify.py.patch
   patch -p0 -d software/CRISPRidentify.patched/components/ < provenance/patches/components_evaluation.py.patch
   cp provenance/patches/_resources.py.new software/CRISPRidentify.patched/components/_resources.py
   ```

4. Build the pilot manifest + stage FASTAs:
   ```
   python scripts/build_pilot_manifest.py --config config/step1.yaml --n 1000
   python scripts/step1.py prepare --config config/step1.yaml
   ```

5. Submit Step 1 + Step 2 via `scripts/sbatch/*.sh` wrappers (see
   `step1.md` for the full procedure).

6. Rebuild the coordinate + sequence catalogs:
   ```
   sbatch scripts/sbatch/coordinate_catalog.sh
   sbatch scripts/sbatch/finalize_coordinate_catalog.sh
   sbatch scripts/sbatch/sequence_catalog.sh
   sbatch scripts/sbatch/export_pairs.sh
   ```

Every entry point enforces a `SLURM_JOB_ID` gate; direct invocation on a
login node will exit with code 78.
