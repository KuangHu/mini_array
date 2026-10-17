# Draft upstream issue — CRISPRidentify feature-extraction defects

Draft only — NOT auto-submitted. Review before posting.

Target: https://github.com/BackofenLab/CRISPRidentify/issues

Repro inputs below use CRISPRidentify's own `TestInput/NC_006513.fa` so no
private project paths or data are shared.

---

## Title

Silent feature-extraction failures depress all Certainty Scores (`hmmr_score`,
`blast_score_{1,2}` permanently 0.0) at commit `474896ce`

## Observed behavior

Running CRISPRidentify at commit `474896ce3c9b89a609e69ac9d2fcdf78e4e2aad4`
produces output where the ML classifier's `hmmr_score`, `blast_score_1`, and
`blast_score_2` features are **always** 0.0 regardless of the input, even on
very large, high-identity arrays that other tools (CRISPRCasTyper, CRT) call
`Trusted`. This systematically caps the final Certainty Score below the
bona-fide threshold.

## Minimal reproducer

```bash
# Fresh clone + env per the shipped environment.yml
git clone --recurse-submodules https://github.com/BackofenLab/CRISPRidentify.git
cd CRISPRidentify
git rev-parse HEAD    # 474896ce...

# Any working directory other than the clone root reproduces the bug.
mkdir /tmp/cri_run && cd /tmp/cri_run

python <path-to>/CRISPRidentify/CRISPRidentify.py \
  --file <path-to>/CRISPRidentify/TestInput/NC_006513.fa \
  --result_folder /tmp/cri_result \
  --model ALL --cas False --is_element False
```

Inspect any `*_Candidates.txt` under `/tmp/cri_result/*/`. The `Array features`
block will contain `blast_score_1: 0.0`, `blast_score_2: 0.0`, `hmmr_score: 0.0`
for every candidate, even for strong arrays with hundreds of near-identical
repeats.

## Root causes (3 distinct defects)

### 1. `tools/hmm_search/models_tandem.hmm` is missing from the repo

`components_evaluation.py:303, 872` reference `tools/hmm_search/models_tandem.hmm`
but the file is not shipped — `ls tools/` shows only `blasting/`,
`CRISPRcasIdentifier/`, `strand_prediction/`. The file is still present in the
sibling `CRISPRloci` repo at
`CRISPRidentify/tools/hmm_search/models_tandem.hmm` (commit
`d25287b4...`, SHA256 `818363ec57f0b3...`, 79 HMM profiles).

Suggested fix: ship `models_tandem.hmm` with CRISPRidentify (same file as in
CRISPRloci) or document that it must be fetched separately.

### 2. BLAST invocations use cwd-dependent relative paths

`components_evaluation.py:406, 419`:

```python
no_binary_cmd = f"blastn -query file_with_all_consensus.fa -db tools/blasting/{db_file} ..."
```

The `-db tools/blasting/...` path only resolves when the process cwd is
CRISPRidentify's own source directory. Any wrapper that `cd`s elsewhere
breaks BLAST silently (`subprocess.Popen(...).communicate()` discards stderr).

Suggested fix: resolve `tools/blasting/` via `os.path.dirname(__file__)` (an
absolute path), or document that the current working directory MUST be
CRISPRidentify's root.

### 3. `hmm_model` f-string typo on line 900

`components_evaluation.py:900`:

```python
no_binary_cmd = f"tools/hmm_search/hmmsearch --tblout result_hmm_{self.index}.out hmm_model protein_{self.index}.fa"
```

The literal string `hmm_model` is passed to hmmsearch instead of the
`hmm_model` parameter. hmmsearch tries to open a file named `hmm_model` and
fails; `subprocess.Popen(...).communicate()` discards the error.

Suggested fix: `... {hmm_model} ...` (f-string substitution).

## Impact

Affects every real-world run that is not executed with the CRISPRidentify
source directory as its cwd. Scoring behavior converges toward "every array
is Possible" because 3 of ~13 input features are permanently zero. We verified
this on 100 GTDB r220 representative genomes: zero bona-fide calls, even for
arrays other callers classify as Trusted I-C/III-D/… systems with hundreds of
spacers.

## Workaround we are using

Local fork at
<redacted>/CRISPRidentify.patched/ that:
- Resolves `tools/blasting/*` and `tools/hmm_search/models_tandem.hmm` via
  absolute paths relative to the installation.
- Fixes the `hmm_model` f-string.
- Wraps subprocess calls so non-zero exit codes become loud errors.
- Pre-flight validates all required resources + external executables at
  Pipeline init.

Patch diffs available on request.
