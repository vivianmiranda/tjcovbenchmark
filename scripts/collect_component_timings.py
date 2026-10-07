"""Check a completed timing campaign and write a portable publication record."""

import argparse
import json
from pathlib import Path

import numpy as np

from common import sha256, write_json
from time_components import CASES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a fresh publication record")
    run_path = args.run / "timing.json"
    record = json.loads(run_path.read_text())
    if record["status"] != "completed":
        raise ValueError("Timing campaign has not completed")
    if len(record["runs"]) != len(CASES)*record["repeats"]:
        raise ValueError("Incomplete timing campaign")
    record["journal_sha256"] = sha256(run_path)
    record["policy"] = (
        "Sequential fresh processes, first-use tables included; no disk-cache "
        "hits; benchmark exports and plotting excluded; the native TJPCov "
        "SSC block-cache write is included")
    record["scope_audit"] = (
        "The original journal's generic file-writing exclusion was too broad. "
        "TJPCov writes a small SSC cache inside get_covariance_block; its full "
        "call is timed unchanged. This corrects the description, not the timings.")
    preflight = args.run / "preflight.json"
    if preflight.exists():
        record["preflight"] = json.loads(preflight.read_text())
        record["preflight_sha256"] = sha256(preflight)
    record["cases"] = {}
    allowed_changes = {"core", "lsst_y1", "threads", "script_sha256",
                       "timing_scope", "timing",
                       "native_block_seconds_diagnostic_only"}
    for run in record["runs"]:
        name = run["case"]
        folder = args.run / f"{name}_r{run['repeat']}"
        saved = args.archive / name
        if run["status"] != "completed" or run["returncode"] != 0:
            raise ValueError("A timing process failed")
        new_path, old_path = folder/"manifest.json", saved/"manifest.json"
        if (sha256(new_path) != run["manifest_sha256"]
                or sha256(old_path) != run["archive_manifest_sha256"]):
            raise ValueError("A manifest changed after measurement")
        new, old = json.loads(new_path.read_text()), json.loads(old_path.read_text())
        for key in set(new) | set(old):
            if key not in allowed_changes and new.get(key) != old.get(key):
                raise ValueError(f"Scientific setting or fingerprint changed: {name}/{key}")
        script = Path(__file__).parent / Path(run["command"][1]).name
        if sha256(script) != run["script_sha256"]:
            raise ValueError(f"Timing script changed: {script}")
        for filename in new["files"]:
            if (sha256(folder/filename) != new["files"][filename]
                    or sha256(saved/filename) != old["files"][filename]):
                raise ValueError("Saved scientific output changed")
            if not filename.endswith(".npz"):
                continue
            with np.load(folder/filename) as a, np.load(saved/filename) as b:
                if set(a.files) != set(b.files):
                    raise ValueError("Array keys changed")
                for key in b.files:
                    if (a[key].dtype != b[key].dtype or a[key].shape != b[key].shape
                            or a[key].tobytes() != b[key].tobytes()):
                        raise ValueError(f"Bitwise discrepancy: {name}/{filename}/{key}")
        run["all_arrays_bitwise_equal_to_accuracy_archive"] = True
        run["scientific_files"] = new["files"]
        if name not in record["cases"]:
            # Keep each case's full settings/source provenance once, rather
            # than repeating it for every timing sample.
            record["cases"][name] = {
                key: value for key, value in new.items()
                if key not in ("timing", "files", "native_block_seconds_diagnostic_only")}
    record["ratios_tjpcov_over_cocoa"] = {}
    for a, b in zip(CASES[::2], CASES[1::2]):
        record["ratios_tjpcov_over_cocoa"][a] = (
            record["summary"][b]["construction_seconds"]["mean"]
            / record["summary"][a]["construction_seconds"]["mean"])
    # Replace the local checkout root in commands/provenance with a portable
    # placeholder. Hashes still refer to the exact original files.
    root = str(Path(__file__).resolve().parents[2])
    record = json.loads(json.dumps(record).replace(root, "<workspace>"))
    write_json(args.output, record)
    print(f"Verified {len(record['runs'])} processes; all scientific arrays bitwise equal")


if __name__ == "__main__":
    main()
