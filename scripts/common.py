"""File provenance and the small LSST input contract shared by the runners."""

import hashlib
import json
import subprocess
from pathlib import Path


def sha256(path):
    """Hash file bytes so reused inputs can be checked before calculation."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def revision(path):
    """Record the Git commit and all local changes, without changing Git."""
    command = ["git", "-C", str(path)]
    return {
        "commit": subprocess.check_output(
            command + ["rev-parse", "HEAD"], text=True).strip(),
        "status": subprocess.check_output(
            command + ["status", "--porcelain"], text=True).strip(),
    }


def write_json(path, values):
    """Save a human-readable record; reject NaN rather than hiding failure."""
    Path(path).write_text(json.dumps(values, indent=2, allow_nan=False) + "\n")


def load_bundle(path, schema):
    """Read a manifest only after checking every named scientific input."""
    path = Path(path)
    record = json.loads((path / "manifest.json").read_text())
    if record["schema"] != schema:
        raise ValueError(f"Expected {schema}, got {record['schema']}")
    for name, expected in record["files"].items():
        if sha256(path / name) != expected:
            raise ValueError(f"Input changed: {path / name}")
    return record


def require_thread_environment():
    """Require an explicit OpenMP allocation and the single-thread BLAS rule."""
    import os

    count = int(os.environ.get("OMP_NUM_THREADS", "0"))
    if count < 1:
        raise ValueError("Set OMP_NUM_THREADS before running the comparison")
    if os.environ.get("OPENBLAS_NUM_THREADS") != "1":
        raise ValueError("Set OPENBLAS_NUM_THREADS=1 before this calculation")
    return count
