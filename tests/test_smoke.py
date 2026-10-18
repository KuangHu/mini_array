"""Smoke tests: config loads, directory layout exists."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path, PurePosixPath

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config" / "step1.yaml"


def test_config_loads():
    cfg = yaml.safe_load(CONFIG.read_text())
    assert cfg["schema_version"] == 1
    assert PurePosixPath(cfg["project_root"]).is_absolute()
    assert PurePosixPath(cfg["code_root"]).is_absolute()


def test_project_root_tree_present():
    cfg = yaml.safe_load(CONFIG.read_text())
    p = Path(cfg["project_root"])
    if not p.is_dir():
        pytest.skip("Cluster integration check: configured project root is not mounted")
    # results/step1 has been split into variant subdirs (patched/broken etc.);
    # check only the stable skeleton here.
    for sub in ("software", "envs", "input/staged", "work", "raw/identify",
                "raw/cctyper", "results", "provenance", "logs"):
        assert (p / sub).is_dir(), f"missing {p / sub}"


def test_doctor_runs():
    """`step1.py doctor` should at least exit without crashing (status may be nonzero
    while envs/manifest aren't fully populated)."""
    result = subprocess.run(
        [sys.executable, str(REPO / "scripts/step1.py"), "doctor",
         "--config", str(CONFIG)],
        capture_output=True, text=True,
    )
    assert "checks passed" in result.stdout, (result.stdout, result.stderr)
