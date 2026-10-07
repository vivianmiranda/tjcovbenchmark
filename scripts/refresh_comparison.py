"""Run the existing small TJPCov comparison stages sequentially.

Launch in the activated CoCoA environment. Supply the separate TJPCov
Python executable without resolving its virtual-environment symlink.
The default campaign writes only to a fresh work directory. --list prints
the exact commands without importing either code or starting a job.

Gaussian means the complete 5x5/30x30 Gaussian matrices, with signal and
noise pieces. SSC means complete 5x5 point-sampled shear SSC. The separated
halo trispectra are not projected cNG matrices. No G+SSC+cNG total or native
real-space SSC/cNG is manufactured when that comparison is not implemented.
"""

import argparse
import json
import os
import shlex
import signal
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

from common import require_thread_environment, sha256, write_json


GROUPS = ("inputs", "gaussian", "ssc", "halo", "trispectrum")
SCRIPTS = Path(__file__).resolve().parent


def plan(work, cocoa, tjpcov):
    """List the existing native exports, matching diagnostics and figures."""
    steps = []
    inputs = work / "lsst_y1"

    def add(group, name, environment, script, arguments, record=None):
        steps.append({
            "group": group, "name": name, "environment": environment,
            "script": script, "arguments": list(map(str, arguments)),
            "record": str(record) if record is not None else None,
        })

    def export(group, name, environment, script, arguments):
        add(group, name, environment, script,
            [*arguments, "--output", work / name],
            work / name / "manifest.json")

    export("inputs", "lsst_y1", "cocoa", "export_lsst_y1.py",
           ["--cocoa", cocoa])
    for kind in ("shear", "3x2"):
        for band, limits in (("low", [30, 150]), ("high", [1500, 1620])):
            name = f"gaussian_{kind}_{band}"
            arguments = [inputs, "--tjpcov", tjpcov, "--ell-range", *limits]
            if kind == "3x2":
                arguments.append("--include-lenses")
            export("gaussian", name, "tjpcov", "run_gaussian.py", arguments)
            comparison = work / f"comparison_{kind}_{band}"
            add("gaussian", comparison.name, "cocoa", "compare_gaussian.py",
                [work / name, "--cocoa", cocoa, "--output", comparison],
                comparison / "comparison.json")
            add("gaussian", f"plot_{kind}_{band}", "tjpcov", "plot_gaussian.py",
                [comparison, "--output", work / "figures" / name])

    # Native CCL time/window refinement uses one fixed supplied power grid.
    # Its N_K test refines its own response table, not the CAMB bundle.
    native_cases = [
        ("low_qag", "low", 1, 1, "qag_quad"),
        *[(f"low_a{n}", "low", n, 1, "qag_quad") for n in (2, 4, 8, 16)],
        ("low_k2", "low", 1, 2, "qag_quad"),
        ("low_a8_spline", "low", 8, 1, "spline"),
        ("high_qag", "high", 1, 1, "qag_quad"),
        ("high_a8", "high", 8, 1, "qag_quad"),
        ("high_a16", "high", 16, 1, "qag_quad"),
    ]
    for label, band, a_refinement, k_refinement, method in native_cases:
        export("ssc", f"ssc_{label}", "tjpcov", "run_ssc.py", [
            inputs, work / f"gaussian_shear_{band}", "--tjpcov", tjpcov,
            "--match-power-a-range", "--a-refinement", a_refinement,
            "--k-refinement", k_refinement, "--integration-method", method,
        ])
    for band in ("low", "high"):
        for level in range(3):
            export("ssc", f"cocoa_ssc_{band}_i{level}", "cocoa",
                   "run_cocoa_ssc.py", [
                       inputs, work / f"gaussian_shear_{band}",
                       "--cocoa", cocoa, "--accuracy-boost", 1,
                       "--integration-accuracy", level,
                   ])
        for label, level, native in (("base", 0, "qag"), ("fine", 2, "a16")):
            name = f"ssc_comparison_{band}_{label}"
            add("ssc", name, "tjpcov", "compare_ssc.py", [
                work / f"cocoa_ssc_{band}_i{level}", work / f"ssc_{band}_{native}",
                "--output", work / name, "--figures", work / "figures" / name,
            ], work / name / "comparison.json")
        export("ssc", f"ssc_models_{band}_ccl", "tjpcov",
               "diagnose_ssc_models.py", [
                   "ccl", inputs, work / f"gaussian_shear_{band}",
                   work / f"ssc_{band}_a16", work / f"cocoa_ssc_{band}_i2",
                   "--tjpcov", tjpcov,
               ])
        export("ssc", f"ssc_models_{band}", "cocoa",
               "diagnose_ssc_models.py", [
                   "cocoa", work / f"ssc_models_{band}_ccl",
                   work / f"cocoa_ssc_{band}_i2", "--cocoa", cocoa,
               ])
    name = "ssc_sampling_a2"
    add("ssc", name, "tjpcov", "diagnose_ssc_sampling.py", [
        inputs, work / "gaussian_shear_low", work / "ssc_low_qag",
        work / "ssc_low_a2", "--tjpcov", tjpcov, "--output", work / name,
    ], work / name / "report.json")
    add("ssc", "collect_ssc", "tjpcov", "collect_ssc_refresh.py", [
        work, "--output", work / "results/ssc",
    ], work / "results/ssc/comparison.json")
    add("ssc", "plot_ssc_sampling", "tjpcov", "plot_ssc_sampling.py", [
        work / "results/ssc/comparison.json", "--output", work / "figures/ssc",
    ])
    add("ssc", "plot_ssc_models", "tjpcov", "plot_ssc_models.py", [
        "--low", work / "ssc_models_low", "--high", work / "ssc_models_high",
        "--low-ccl", work / "ssc_models_low_ccl",
        "--high-ccl", work / "ssc_models_high_ccl",
        "--results", work / "results/ssc_models",
        "--figures", work / "figures/ssc_models",
    ], work / "results/ssc_models/comparison.json")
    add("ssc", "galaxy_bias_placement", "tjpcov", "diagnose_galaxy_bias.py", [
        inputs, work / "gaussian_3x2_low", "--tjpcov", tjpcov,
        "--output", work / "galaxy_bias_placement",
    ], work / "galaxy_bias_placement/report.json")

    for count, refinement in ((128, 1), (255, 2)):
        export("halo", f"halo_native_m{count}", "tjpcov", "run_halo.py",
               [inputs, "--tjpcov", tjpcov, "--mass-refinement", refinement])
    for level in range(3):
        export("halo", f"halo_cocoa_i{level}", "cocoa", "export_cocoa_halo.py",
               [inputs, work / "halo_native_m255", "--cocoa", cocoa,
                "--integration-accuracy", level])
        name = f"halo_comparison_i{level}"
        add("halo", name, "tjpcov", "compare_halo.py", [
            work / "halo_native_m255", work / f"halo_cocoa_i{level}",
            "--output", work / name, "--figures", work / "figures" / name,
        ], work / name / "comparison.json")

    for level in range(3):
        export("trispectrum", f"trispectrum_cocoa_i{level}", "cocoa",
               "export_cocoa_trispectrum.py", [
                   inputs, "--cocoa", cocoa, "--integration-accuracy", level,
               ])
    native_cases = [
        ("native_m128", 1, []), ("native_m255", 2, []),
        ("direct_m255", 2, ["--direct-growth"]),
        ("bhattacharya_m255", 2, ["--concentration", "bhattacharya13"]),
        ("tinker10_m255", 2, ["--mass-function", "tinker10"]),
        ("both_m255", 2, ["--concentration", "bhattacharya13",
                          "--mass-function", "tinker10"]),
    ]
    for label, refinement, options in native_cases:
        export("trispectrum", f"trispectrum_{label}", "tjpcov",
               "run_trispectrum.py", [
                   inputs, work / "trispectrum_cocoa_i0", "--tjpcov", tjpcov,
                   "--mass-refinement", refinement, *options,
               ])
    for label in ("native", "direct"):
        export("trispectrum", f"trispectrum_comparison_{label}", "cocoa",
               "compare_trispectrum.py", [
                   work / "trispectrum_cocoa_i0",
                   work / f"trispectrum_{label}_m255",
               ])
    name = "trispectrum_power_refinement"
    add("trispectrum", name, "cocoa", "diagnose_trispectrum_power.py", [
        inputs, work / "trispectrum_cocoa_i0", "--cocoa", cocoa,
        "--factors", 1, 2, "--output", work / name,
    ], work / name / "report.json")
    add("trispectrum", "plot_trispectrum", "tjpcov", "plot_trispectrum.py", [
        work, "--output", work / "results/trispectrum",
        "--figures", work / "figures/trispectrum",
    ], work / "results/trispectrum/comparison.json")
    return steps


def verify_record(path):
    """Accept completion only after checking the saved report and file hashes."""
    record = json.loads(path.read_text())
    if (record.get("status", "completed") != "completed"
            or not record.get("passed", True)):
        raise RuntimeError(f"Scientific validation did not pass: {path}")
    for name, expected in record.get("files", {}).items():
        if sha256(path.parent / name) != expected:
            raise RuntimeError(f"Saved output hash mismatch: {path.parent / name}")
    return sha256(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--tjpcov-python", type=Path, required=True)
    parser.add_argument("--cocoa-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--groups", nargs="+", choices=GROUPS, default=GROUPS)
    parser.add_argument("--start-at", help="explicit first remaining stage name")
    parser.add_argument("--timeout", type=int, default=600,
                        help="hard per-stage limit, at most 600 seconds")
    parser.add_argument("--list", action="store_true", help="print commands only")
    args = parser.parse_args()
    if not 0 < args.timeout <= 600:
        parser.error("Per-stage timeout must be 1..600 seconds")
    work = args.output.resolve()
    steps = [s for s in plan(work, args.cocoa.resolve(), args.tjpcov.resolve())
             if s["group"] in args.groups]
    if args.start_at:
        names = [step["name"] for step in steps]
        if args.start_at not in names:
            parser.error("--start-at must name a stage in the selected groups")
        steps = steps[names.index(args.start_at):]
    # absolute(), unlike resolve(), retains .local/bin/python and therefore
    # the correct virtual environment instead of its base interpreter.
    python = {"cocoa": args.cocoa_python.absolute(),
              "tjpcov": args.tjpcov_python.absolute()}
    for step in steps:
        step["command"] = [str(python[step["environment"]]),
                           str(SCRIPTS / step["script"]), *step["arguments"]]
        step["script_sha256"] = sha256(SCRIPTS / step["script"])
    if args.list:
        for step in steps:
            print(f"[{step['group']}] {step['name']}")
            print(shlex.join(step["command"]))
        return
    threads = require_thread_environment()
    if threads > 6:
        parser.error("This laptop correctness campaign permits at most 6 threads")
    if not all(path.is_file() for path in python.values()):
        parser.error("Both existing Python executables are required")
    work.mkdir(parents=True, exist_ok=True)
    (work / "logs").mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    journal = work / f"refresh_{stamp}.json"
    record = {"schema": "tjpcov-production-power-refresh-v1", "status": "running",
              "threads": threads, "timeout_seconds_per_stage": args.timeout,
              "timing_policy": "Accuracy refresh; stage times are not benchmarks",
              "runner_sha256": sha256(__file__), "stages": steps}
    write_json(journal, record)

    # Each fresh child owns its caches and power tables. Waiting for that
    # process before spawning the next keeps numerical work sequential.
    # Never reuse a TJPCov disk-cache folder after changing power inputs.
    for step in steps:
        log = work / "logs" / f"{step['name']}_{stamp}.log"
        environment = dict(os.environ, PYTHONUNBUFFERED="1")
        if step["environment"] == "tjpcov":
            for key in ("PYTHONPATH", "PYTHONHOME", "LD_LIBRARY_PATH",
                        "DYLD_LIBRARY_PATH", "VIRTUAL_ENV", "ROOTDIR",
                        "CONDA_PREFIX"):
                environment.pop(key, None)
            environment["PATH"] = (str(python["tjpcov"].parent) + os.pathsep
                                   + environment.get("PATH", ""))
        step.update(status="running", log=str(log))
        write_json(journal, record)
        print(f"Starting {step['name']} (limit {args.timeout} s)", flush=True)
        started = perf_counter()
        try:
            if sha256(SCRIPTS / step["script"]) != step["script_sha256"]:
                raise RuntimeError("Stage script changed after the plan was saved")
            with log.open("x") as stream:
                process = subprocess.Popen(
                    step["command"], stdout=stream, stderr=subprocess.STDOUT,
                    env=environment, cwd=SCRIPTS.parent, start_new_session=True,
                )
                try:
                    returncode = process.wait(timeout=args.timeout)
                except (subprocess.TimeoutExpired, KeyboardInterrupt):
                    if process.poll() is None:
                        os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    raise
            if returncode != 0:
                raise RuntimeError(f"Exit {returncode}; inspect {log}")
            if step["record"]:
                step["record_sha256"] = verify_record(Path(step["record"]))
            step["status"] = "completed"
        except (Exception, KeyboardInterrupt) as error:
            step.update(status="failed", error=str(error))
            record["status"] = "stopped"
            raise
        finally:
            step["wall_seconds_diagnostic_only"] = perf_counter()-started
            write_json(journal, record)
    record["status"] = "completed"
    write_json(journal, record)
    print(f"Selected stages finished; review their scientific results: {journal}")


if __name__ == "__main__":
    main()
