"""Feed TJPCov's saved spectra and band operator to CoCoA's real C kernels.

Run in the Cocoa environment. No survey initialization or CAMB call is
needed: both implementations now receive the identical signal and noise.
Every entry of the five-band Gaussian matrix is retained. There are no
SSC/cNG arrays here because those components have not been calculated.
"""

import argparse
import json
import sys
from pathlib import Path

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error("Choose a new output directory to preserve comparisons")
    threads = require_thread_environment()
    run = args.run.resolve()
    schema = json.loads((run / "manifest.json").read_text())["schema"]
    if schema not in ("tjpcov-gaussian-shear-v1", "tjpcov-gaussian-fields-v1"):
        raise ValueError("Expected a native shear or three-field Gaussian pilot")
    record = load_bundle(run, schema)
    if record["status"] != "completed":
        raise ValueError("TJPCov run did not complete")

    import numpy as np
    from scipy.linalg import eigvalsh

    data = np.load(run / "gaussian.npz", allow_pickle=False)
    ell, spectra, noise = data["ell"], data["spectra"], data["noise_power"]
    operators, pairs = data["operators"], data["pairs"]
    expected_pairs = [[0, 0]]
    nfield = 1
    if schema == "tjpcov-gaussian-fields-v1":
        nfield = 3
        expected_pairs = [[0, 0], [0, 1], [1, 1], [0, 2], [1, 2], [2, 2]]
    if (ell.ndim != 1 or not np.all(np.diff(ell) == 1)
            or not np.all(ell == ell.astype(int)) or ell[0] < 2
            or spectra.shape != (ell.size, nfield, nfield)
            or noise.shape != (nfield,) or operators.shape != (5, ell.size)
            or not np.array_equal(pairs, expected_pairs)):
        raise ValueError("Malformed ell, field, observable-pair or band layout")
    for values in (spectra, noise, operators):
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("Invalid input spectrum, noise or operator")
    np.testing.assert_allclose(operators.sum(axis=1), 1, rtol=0, atol=1e-14)

    cocoa = args.cocoa.resolve()
    sys.path.insert(0, str(cocoa / "projects/lsst_y1"))
    import cosmolike_lsst_y1_interface as interface

    if not getattr(interface, "has_covariance", False):
        raise RuntimeError("Compile LSST Y1 with covariance support first")
    compute = interface.covariance.covariance_gaussian_fourier
    arguments = dict(pairs=pairs, operators=operators, ell_min=int(ell[0]),
                     area_sr=record["area_sr"])
    cocoa_components = {
        "total": compute(spectra=spectra, noise=noise, **arguments),
        "sample_variance": compute(spectra=spectra,
                                   noise=np.zeros_like(noise), **arguments),
        "noise": compute(spectra=np.zeros_like(spectra), noise=noise,
                         **arguments),
    }
    cocoa_components["mixed"] = (cocoa_components["total"]
                                 - cocoa_components["sample_variance"]
                                 - cocoa_components["noise"])
    tjpcov_components = {key: data[key] for key in cocoa_components}
    rms = np.sqrt(np.diag(cocoa_components["total"]))
    if not np.all(np.isfinite(rms)) or np.any(rms <= 0):
        raise ValueError("CoCoA total has invalid variances")

    # Use the same CoCoA total rms product for every component. This keeps
    # tiny/zero cross entries in the comparison without division by them.
    # Fractional errors on nonzero entries are a separate diagnostic.
    tolerance = 1e-11
    size = len(pairs) * operators.shape[0]
    report = {"schema": "cocoa-tjpcov-gaussian-comparison-v1",
              "scope": record["scope"], "native_run": record,
              "native_manifest_sha256": sha256(run / "manifest.json"),
              "component_diagnostics": {}, "passed": True,
              "variance_scaled_tolerance": tolerance,
              "shape": [size, size], "threads": threads,
              "core": revision(cocoa / "external_modules/code/cosmolike_core"),
              "lsst_y1": revision(cocoa / "projects/lsst_y1"),
              "interface_sha256": sha256(interface.__file__),
              "script_sha256": sha256(__file__)}
    for name, reference in cocoa_components.items():
        candidate = tjpcov_components[name]
        if candidate.shape != reference.shape or not np.all(np.isfinite(candidate)):
            raise ValueError(f"Invalid TJPCov component: {name}")
        np.testing.assert_allclose(candidate, candidate.T, rtol=1e-13, atol=0)
        difference = candidate - reference
        scaled = difference / rms[:, None] / rms[None, :]
        nonzero = reference != 0
        maximum = float(np.max(np.abs(scaled)))
        report["component_diagnostics"][name] = {
            "max_variance_scaled_difference": maximum,
            "max_fractional_nonzero_difference": float(np.max(
                np.abs(difference[nonzero] / reference[nonzero])))
                if np.any(nonzero) else None,
        }
        report["passed"] &= maximum < tolerance

    positive = True
    for name, components in (("cocoa", cocoa_components),
                              ("tjpcov", tjpcov_components)):
        scaled = components["total"] / rms[:, None] / rms[None, :]
        minimum = float(eigvalsh(scaled)[0])
        report[f"{name}_min_scaled_eigenvalue"] = minimum
        positive &= minimum > 0
    report["passed"] &= positive
    if positive:
        ratios = eigvalsh(tjpcov_components["total"], cocoa_components["total"])
        report["generalized_variance_ratio_range"] = [
            float(ratios[0]), float(ratios[-1])]
        report["max_generalized_variance_change"] = float(
            np.max(np.abs(ratios - 1)))
    else:
        report["generalized_variance_ratio_range"] = None

    output.mkdir(parents=True)
    arrays = {"ell_edges": data["edges"], "effective_ell": data["effective_ell"],
              "pairs": pairs}
    for label, components in (("cocoa", cocoa_components),
                               ("tjpcov", tjpcov_components)):
        arrays.update({f"{label}_{key}": value for key, value in components.items()})
    np.savez_compressed(output / "matrices.npz", **arrays)
    report["files"] = {"matrices.npz": sha256(output / "matrices.npz")}
    write_json(output / "comparison.json", report)
    print(f"Shared-spectrum Gaussian assembly passed: {report['passed']}")
    if not report["passed"]:
        raise SystemExit(f"Comparison failed; inspect {output / 'comparison.json'}")


if __name__ == "__main__":
    main()
