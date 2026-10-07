"""Archive the refreshed small SSC campaign and its numerical refinements.

Each saved matrix has five effective multipoles and all 25 entries.
This collector runs neither covariance code. It keeps native table
refinement separate from the difference between the physical models.
"""

import argparse
from pathlib import Path

from common import load_bundle, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a fresh archive directory")

    import numpy as np

    records, arrays, matrices = {}, {}, {}
    names = ["ssc_low_qag", "ssc_low_a2", "ssc_low_a4", "ssc_low_a8",
             "ssc_low_a16", "ssc_low_k2", "ssc_low_a8_spline",
             "ssc_high_qag", "ssc_high_a8", "ssc_high_a16"]
    names += [f"cocoa_ssc_{band}_i{level}"
              for band in ("low", "high") for level in range(3)]
    input_hash = None
    for name in names:
        folder = args.work / name
        code = "cocoa" if name.startswith("cocoa") else "tjpcov"
        record = load_bundle(folder, f"{code}-ssc-shear-v1")
        if record["status"] != "completed":
            raise ValueError(f"Incomplete native SSC: {name}")
        if input_hash is None:
            input_hash = record["input_manifest_sha256"]
        if record["input_manifest_sha256"] != input_hash:
            raise ValueError("SSC cases do not share one power bundle")
        with np.load(folder / "ssc.npz", allow_pickle=False) as data:
            for key in data.files:
                if not np.all(np.isfinite(data[key])):
                    raise ValueError(f"Nonfinite {name}/{key}")
                arrays[f"{name}__{key}"] = data[key]
            matrices[name] = data["ssc"]
        if (matrices[name].shape != (5, 5)
                or np.any(np.diag(matrices[name]) <= 0)):
            raise ValueError(f"Invalid five-multipole SSC: {name}")
        if code == "cocoa" and record["settings"]["accuracy_boost"] != 1:
            raise ValueError("This campaign holds CoCoA global boost at one")
        records[name] = {
            "manifest_sha256": sha256(folder / "manifest.json"),
            "manifest": record,
        }

    # Normalize a refinement by the rms products of its finer matrix.
    # This remains useful for small off-diagonal entries and retains the
    # raw arrays, including any negative modes, for subsequent inspection.
    comparisons = []
    pairs = [(f"ssc_low_{a}", f"ssc_low_{b}") for a, b in
             (("qag", "a2"), ("a2", "a4"), ("a4", "a8"),
              ("a8", "a16"), ("qag", "k2"), ("a8", "a8_spline"))]
    pairs += [("ssc_high_qag", "ssc_high_a8"),
              ("ssc_high_a8", "ssc_high_a16")]
    pairs += [(f"cocoa_ssc_{band}_i{level}",
               f"cocoa_ssc_{band}_i{level+1}")
              for band in ("low", "high") for level in range(2)]
    for coarse, fine in pairs:
        np.testing.assert_array_equal(arrays[f"{coarse}__ell"],
                                      arrays[f"{fine}__ell"])
        if (records[coarse]["manifest"]["gaussian_manifest_sha256"]
                != records[fine]["manifest"]["gaussian_manifest_sha256"]):
            raise ValueError("Refinement changed the effective multipoles")
        rms = np.sqrt(np.diag(matrices[fine]))
        scaled = (matrices[coarse]-matrices[fine])/rms[:, None]/rms[None, :]
        comparisons.append({
            "coarse": coarse, "fine": fine,
            "max_abs_difference_over_fine_rms_product": float(
                np.max(np.abs(scaled))),
            "diagonal_fractional_change": (
                np.diag(matrices[coarse])/np.diag(matrices[fine])-1).tolist(),
        })

    args.output.mkdir(parents=True)
    np.savez_compressed(args.output / "ssc_native.npz", **arrays)
    write_json(args.output / "comparison.json", {
        "schema": "cocoa-tjpcov-native-ssc-campaign-v1",
        "scope": "Native five-multipole shear SSC; no G/cNG/total matrix",
        "input_manifest_sha256": input_hash,
        "cocoa_accuracy": "AB=1; independent integration levels 0,1,2",
        "timing_policy": "Accuracy campaign; elapsed logs are not benchmarks",
        "refinements": comparisons, "runs": records,
        "script_sha256": sha256(__file__),
        "files": {"ssc_native.npz": sha256(args.output / "ssc_native.npz")},
    })
    print(args.output / "comparison.json")


if __name__ == "__main__":
    main()
