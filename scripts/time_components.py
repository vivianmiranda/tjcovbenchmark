"""Time checked SSC/trispectrum pilots in fresh, sequential processes.

Run from the activated CoCoA environment on an otherwise quiet machine.
Supply an existing accuracy campaign and the separate TJPCov interpreter.
Every numerical array must reproduce its accuracy archive exactly before
the measured construction times are accepted. No full-covariance estimate
is derived from these component pilots.
"""

import argparse
import json
import os
import platform
import re
import signal
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np

from common import require_thread_environment, sha256, write_json
from refresh_comparison import SCRIPTS, plan, verify_record


CASES = ("cocoa_ssc_low_i2", "ssc_low_a16",
         "cocoa_ssc_high_i2", "ssc_high_a16",
         "trispectrum_cocoa_i2", "trispectrum_native_m255")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--tjpcov-python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    threads = require_thread_environment()
    if threads > 8 or args.repeats < 2 or args.output.exists():
        parser.error("Use at most eight threads, at least two repeats, "
                     "and a fresh output directory")
    archive, output = args.archive.resolve(), args.output.resolve()
    stages = {s["name"]: s for s in plan(
        archive, args.cocoa.resolve(), args.tjpcov.resolve())}
    python = {"cocoa": sys.executable,
              "tjpcov": str(args.tjpcov_python.absolute())}
    output.mkdir(parents=True)
    record = {
        "schema": "tjpcov-component-timings-v1", "status": "running",
        "threads": threads, "repeats": args.repeats,
        "platform": platform.platform(), "runner_sha256": sha256(__file__),
        "policy": "Sequential fresh processes, first-use tables included; "
                  "no disk-cache hits; benchmark exports and plotting excluded; "
                  "the native TJPCov SSC block-cache write is included",
        "environment": {k: os.environ.get(k) for k in (
            "OMP_NUM_THREADS", "OMP_PROC_BIND", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
        "input_manifest_sha256": sha256(archive / "lsst_y1/manifest.json"),
        "runs": [],
    }
    journal = output / "timing.json"
    write_json(journal, record)
    try:
        for repeat in range(args.repeats):
            # Alternate code order between rounds to expose systematic drift.
            names = CASES if repeat % 2 == 0 else tuple(reversed(CASES))
            for name in names:
                stage = stages[name]
                folder = output / f"{name}_r{repeat+1}"
                arguments = stage["arguments"].copy()
                arguments[arguments.index("--output")+1] = str(folder)
                script = SCRIPTS / stage["script"]
                command = [python[stage["environment"]], str(script),
                           *arguments, "--timing"]
                environment = dict(os.environ, PYTHONUNBUFFERED="1")
                if stage["environment"] == "tjpcov":
                    for key in ("PYTHONPATH", "PYTHONHOME", "LD_LIBRARY_PATH",
                                "DYLD_LIBRARY_PATH", "VIRTUAL_ENV", "ROOTDIR",
                                "CONDA_PREFIX"):
                        environment.pop(key, None)
                    environment["PATH"] = (
                        str(Path(python["tjpcov"]).parent) + os.pathsep
                        + environment.get("PATH", ""))
                log = output / f"{folder.name}.log"
                run = {"case": name, "repeat": repeat+1, "command": command,
                       "log": str(log), "status": "running",
                       "script_sha256": sha256(script)}
                record["runs"].append(run)
                write_json(journal, record)
                print(f"Starting {name}, repeat {repeat+1}", flush=True)
                # The external timer records peak child memory. Its wall time
                # includes imports/verification; only the in-script compute
                # timers are used in the public component comparison.
                timed = ["/usr/bin/time", "-l" if sys.platform == "darwin"
                         else "-v", *command]
                with log.open("x") as stream:
                    process = subprocess.Popen(
                        timed, stdout=stream, stderr=subprocess.STDOUT,
                        env=environment, cwd=SCRIPTS.parent,
                        start_new_session=True)
                    try:
                        run["returncode"] = process.wait(timeout=600)
                    except (subprocess.TimeoutExpired, KeyboardInterrupt):
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                        raise
                if run["returncode"]:
                    raise RuntimeError(f"Failed {name}; inspect {log}")
                manifest_path = folder / "manifest.json"
                run["manifest_sha256"] = verify_record(manifest_path)
                run["archive_manifest_sha256"] = verify_record(
                    archive / name / "manifest.json")
                manifest = json.loads(manifest_path.read_text())
                old = json.loads((archive / name / "manifest.json").read_text())
                if manifest["input_manifest_sha256"] != old["input_manifest_sha256"]:
                    raise ValueError("Input fingerprint differs from accuracy run")
                checks = {}
                for filename in manifest["files"]:
                    if not filename.endswith(".npz"):
                        continue
                    with np.load(folder / filename) as new, np.load(
                            archive / name / filename) as saved:
                        if set(new.files) != set(saved.files):
                            raise ValueError("Array keys changed")
                        for key in saved.files:
                            np.testing.assert_array_equal(new[key], saved[key])
                        checks[filename] = {"all_arrays_exact": True,
                                            "arrays": saved.files}
                text = log.read_text()
                match = re.search(r"(\d+)\s+maximum resident set size", text)
                if sys.platform != "darwin":
                    match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
                run["peak_memory_bytes"] = (int(match[1]) * (
                    1 if sys.platform == "darwin" else 1024) if match else None)
                run.update(status="completed", timing=manifest["timing"],
                           accuracy_checks=checks)
                write_json(journal, record)
        record["summary"] = {}
        for name in CASES:
            runs = [r for r in record["runs"] if r["case"] == name]
            record["summary"][name] = {
                field: {"seconds": [r["timing"][field] for r in runs],
                        "mean": statistics.mean(r["timing"][field] for r in runs),
                        "sample_stdev": statistics.stdev(r["timing"][field] for r in runs)}
                for field in ("setup_seconds", "construction_seconds")}
        record["status"] = "completed"
    except (Exception, KeyboardInterrupt) as error:
        record.update(status="stopped", error=str(error))
        raise
    finally:
        write_json(journal, record)
    print(f"Verified timing record: {journal}", flush=True)


if __name__ == "__main__":
    main()
