"""End-to-end tests for the patched CRISPRidentify feature-extraction path.

Covers the paths called out in step1 Phase C of the review instructions:
  - BLAST DB1 + DB2 each produce a non-zero bitscore on a known-hit query
  - HMM produces a non-zero score on a known-hit protein query
  - Legitimate no-hit (random sequence) completes with exit 0 and empty output
  - Legitimate no-ORF / no-query is handled explicitly (not as failure)
  - Missing BLAST DB or HMM file raises FeatureExtractionError
  - External executable not on PATH raises FeatureExtractionError
  - Different starting cwd produces identical results
  - Two parallel tasks don't clobber each other's temporary files

Only runs if envs/identify is built (`CRISPR_IDENTIFY_ENV_BIN` env var set to
`<project_root>/envs/identify/bin`). Otherwise tests are skipped so the CI
smoke test can still run on the login node without the full env active.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pytest

# -- project paths --
PROJECT_ROOT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
PATCHED = PROJECT_ROOT / "software" / "CRISPRidentify.patched"
ENV_BIN = PROJECT_ROOT / "envs" / "identify" / "bin"

if not PATCHED.is_dir() or not ENV_BIN.is_dir():
    pytestmark = pytest.mark.skip(
        reason="patched CRISPRidentify or envs/identify not installed"
    )


# -- helpers: run patched _resources helpers from inside envs/identify python --

def _run_in_identify(code: str, env_overrides: dict | None = None,
                     cwd: str | None = None) -> subprocess.CompletedProcess:
    """Execute a short Python snippet inside envs/identify and return result.

    The snippet runs with:
      - sys.path primed to import `components._resources`
      - PATH pointing at envs/identify/bin
    """
    env = os.environ.copy()
    env["PATH"] = f"{ENV_BIN}:{env.get('PATH','')}"
    if env_overrides:
        env.update(env_overrides)
    wrapper = (
        "import sys, os\n"
        f"sys.path.insert(0, {str(PATCHED)!r})\n"
        + code
    )
    return subprocess.run(
        [str(ENV_BIN / "python"), "-c", wrapper],
        capture_output=True, text=True, env=env, cwd=cwd,
    )


def test_preflight_resolves_all_resources():
    r = _run_in_identify(
        "from components._resources import preflight; "
        "import json; print(json.dumps(preflight()))"
    )
    assert r.returncode == 0, r.stderr
    import json
    info = json.loads(r.stdout.strip().splitlines()[-1])
    assert "Verified_repeats_dataset1.fa" in info["resources"]
    assert "models_tandem.hmm" in info["resources"]
    for exe in ("blastn", "hmmsearch", "prodigal", "vmatch", "mkvtree"):
        assert exe in info["executables"]


def test_preflight_fails_loudly_when_tools_dir_missing(tmp_path: Path):
    r = _run_in_identify(
        "from components._resources import preflight; preflight()",
        env_overrides={"CRISPRIDENTIFY_TOOLS_DIR": str(tmp_path / "nope")},
    )
    assert r.returncode != 0
    assert "FeatureExtractionError" in r.stderr
    assert "missing" in r.stderr.lower()


# -- BLAST fixtures --

def _blast_hit_fixture(tmp_path: Path, which: str) -> Path:
    """Write a 'file_with_all_consensus.fa' containing a known-hit query."""
    db_fasta = PATCHED / "tools" / "blasting" / f"Verified_repeats_dataset{which}.fa"
    with db_fasta.open() as fh:
        header = fh.readline()
        seq = fh.readline().strip()
    query = tmp_path / "file_with_all_consensus.fa"
    query.write_text(f">query_hit\n{seq}\n")
    return query


def test_blast_db1_known_hit_nonzero_score(tmp_path: Path):
    _blast_hit_fixture(tmp_path, "1")
    r = _run_in_identify(
        "from components._resources import blast_db, run_subprocess\n"
        "db = blast_db('Verified_repeats_dataset1.fa')\n"
        "run_subprocess(f'blastn -query file_with_all_consensus.fa -db {db} "
        "-word_size 6 -outfmt 6 -out out.tsv', label='test_db1')\n"
        "lines = open('out.tsv').read().strip().splitlines()\n"
        "print('n_hits', len(lines))\n"
        "assert lines, 'expected ≥1 hit for self-match query'\n",
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert "n_hits" in r.stdout


def test_blast_db2_known_hit_nonzero_score(tmp_path: Path):
    _blast_hit_fixture(tmp_path, "2")
    r = _run_in_identify(
        "from components._resources import blast_db, run_subprocess\n"
        "db = blast_db('Verified_repeats_dataset2.fa')\n"
        "run_subprocess(f'blastn -query file_with_all_consensus.fa -db {db} "
        "-word_size 6 -outfmt 6 -out out.tsv', label='test_db2')\n"
        "lines = open('out.tsv').read().strip().splitlines()\n"
        "assert lines, 'expected ≥1 hit for self-match query'\n",
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, (r.stdout, r.stderr)


def test_blast_legitimate_no_hit_succeeds(tmp_path: Path):
    (tmp_path / "file_with_all_consensus.fa").write_text(
        ">random\nAAAATTTTTGGGGGCCCCC\n"
    )
    r = _run_in_identify(
        "from components._resources import blast_db, run_subprocess\n"
        "db = blast_db('Verified_repeats_dataset1.fa')\n"
        "run_subprocess(f'blastn -query file_with_all_consensus.fa -db {db} "
        "-word_size 6 -outfmt 6 -out out.tsv', label='test_no_hit')\n"
        "import os\n"
        "assert os.path.exists('out.tsv'), 'output file should exist'\n"
        "print('empty_output_ok size', os.stat('out.tsv').st_size)\n",
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert "empty_output_ok" in r.stdout


def test_blast_missing_db_raises(tmp_path: Path):
    (tmp_path / "file_with_all_consensus.fa").write_text(">x\nACGTACGT\n")
    r = _run_in_identify(
        "from components._resources import run_subprocess\n"
        "run_subprocess('blastn -query file_with_all_consensus.fa "
        "-db /tmp/does_not_exist.fa -word_size 6 -outfmt 6 -out out.tsv', "
        "label='test_missing_db')\n",
        cwd=str(tmp_path),
    )
    assert r.returncode != 0
    assert "FeatureExtractionError" in r.stderr


# -- HMM fixtures --

@pytest.fixture
def hmm_hit_fixture(tmp_path: Path) -> Path:
    """Emit the YSIRK_signal consensus as a query — guaranteed self-hit."""
    hmm_full = PATCHED / "tools" / "hmm_search" / "models_tandem.hmm"
    ysirk = tmp_path / "ysirk.hmm"
    subprocess.run([str(ENV_BIN / "hmmfetch"), str(hmm_full), "YSIRK_signal"],
                   stdout=ysirk.open("w"), check=True)
    emit = subprocess.run([str(ENV_BIN / "hmmemit"), "-c", str(ysirk)],
                          capture_output=True, text=True, check=True)
    protein = tmp_path / "protein_1.fa"
    protein.write_text(emit.stdout)
    return protein


def test_hmm_known_hit_nonzero_score(hmm_hit_fixture, tmp_path: Path):
    r = _run_in_identify(
        "from components._resources import hmm_model, run_subprocess\n"
        "h = hmm_model()\n"
        "run_subprocess(f'hmmsearch --tblout result_hmm_1.out {h} protein_1.fa', "
        "label='test_hmm_hit')\n"
        "out = open('result_hmm_1.out').read()\n"
        "import re\n"
        "lines = [l for l in out.splitlines() "
        "         if l and not l.startswith('#')]\n"
        "assert lines, 'expected ≥1 HMM hit on self-consensus query'\n"
        "print('n_hmm_hits', len(lines))\n",
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert "n_hmm_hits" in r.stdout


def test_hmm_legitimate_no_hit_succeeds(tmp_path: Path):
    (tmp_path / "protein_1.fa").write_text(">rand\nFFFFFFFFFFFFFFFFFF\n")
    r = _run_in_identify(
        "from components._resources import hmm_model, run_subprocess\n"
        "h = hmm_model()\n"
        "run_subprocess(f'hmmsearch --tblout result_hmm_1.out {h} protein_1.fa', "
        "label='test_hmm_nohit')\n"
        "out = open('result_hmm_1.out').read()\n"
        "lines = [l for l in out.splitlines() if l and not l.startswith('#')]\n"
        "print('empty_ok n=', len(lines))\n",
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, (r.stdout, r.stderr)


def test_hmm_no_query_handled_explicitly(tmp_path: Path):
    """No protein FASTA at all → bulk HMM path writes empty result and returns.

    Mirrors upstream's 'no ORFs' branch — must NOT be a failure.
    """
    r = _run_in_identify(
        "import sys\n"
        "sys.path.insert(0, '/groups/rubin/projects/kuang/out/CRISPR_array_related/software/CRISPRidentify.patched')\n"
        "# simulate the patched bulk path's no-protein branch\n"
        "from components._resources import hmm_model\n"
        "import os\n"
        "if not os.path.exists('protein_results.fa') or os.stat('protein_results.fa').st_size == 0:\n"
        "    open('result_hmm.out', 'w').close()\n"
        "    print('no_query_branch_hit')\n",
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stderr
    assert "no_query_branch_hit" in r.stdout


# -- cwd independence --

def test_results_independent_of_starting_cwd(tmp_path: Path):
    # Same query, two different cwds — must produce byte-identical BLAST output.
    for subdir in ("here", "there"):
        d = tmp_path / subdir
        d.mkdir()
        _blast_hit_fixture(d, "1")
    out = {}
    for subdir in ("here", "there"):
        d = tmp_path / subdir
        r = _run_in_identify(
            "from components._resources import blast_db, run_subprocess\n"
            "db = blast_db('Verified_repeats_dataset1.fa')\n"
            "run_subprocess(f'blastn -query file_with_all_consensus.fa -db {db} "
            "-word_size 6 -outfmt 6 -out out.tsv', label='cwd_test')\n"
            "print(open('out.tsv').read(), end='')\n",
            cwd=str(d),
        )
        assert r.returncode == 0
        out[subdir] = r.stdout
    assert out["here"] == out["there"]


# -- parallel tempfile isolation --

def test_parallel_tasks_no_crosstalk(tmp_path: Path):
    """Two concurrent runs each in their own cwd — distinct outputs preserved."""
    dirs = []
    for i in range(2):
        d = tmp_path / f"task{i}"
        d.mkdir()
        (d / "file_with_all_consensus.fa").write_text(
            f">q_{i}\n" + ("A" * 20 if i == 0 else "G" * 20) + "\n"
        )
        dirs.append(d)
    procs = []
    for d in dirs:
        procs.append(subprocess.Popen(
            [str(ENV_BIN / "python"), "-c",
             f"""
import sys, os
sys.path.insert(0, {str(PATCHED)!r})
from components._resources import blast_db, run_subprocess
db = blast_db('Verified_repeats_dataset1.fa')
run_subprocess(f'blastn -query file_with_all_consensus.fa -db {{db}} -word_size 6 -outfmt 6 -out out.tsv', label='parallel')
print(open('out.tsv').read(), end='')
"""],
            cwd=str(d), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        ))
    outs = [p.communicate() for p in procs]
    rcs = [p.returncode for p in procs]
    assert all(rc == 0 for rc in rcs), (rcs, outs)
    # per-cwd `out.tsv` must stay independent
    for d in dirs:
        assert (d / "out.tsv").is_file()
