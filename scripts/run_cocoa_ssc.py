"""Compute CoCoA's native shear SSC at the TJPCov effective multipoles.

This is an unbinned SSC model comparison, separate from the band-averaged
Gaussian assembly test. It uses the production C interface and existing
CoCoA geometry/halo helpers. No TJPCov response or mask is substituted.
Run in the Cocoa environment, with an explicit OpenMP thread allocation.
"""

import argparse
import sys
from pathlib import Path

from common import load_bundle, require_thread_environment, revision
from common import sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("gaussian_run", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--integration-accuracy", type=int, choices=range(5),
                        default=0)
    parser.add_argument("--accuracy-boost", type=int, choices=(1, 2, 4, 8),
                        default=1)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a fresh output directory")
    threads = require_thread_environment()
    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    native = load_bundle(args.gaussian_run, "tjpcov-gaussian-shear-v1")
    if native["input_manifest_sha256"] != sha256(args.inputs / "manifest.json"):
        raise ValueError("The TJPCov bands belong to a different input bundle")

    import numpy as np

    cocoa = args.cocoa.resolve()
    project = cocoa / "projects/lsst_y1"
    core = cocoa / "external_modules/code/cosmolike_core"
    sys.path[:0] = [str(project), str(project / "covariance"), str(core)]
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.gaussian import limber_spectra
    from cosmolike_notebook_utils.covariance.geometry import cap_mask
    from cosmolike_notebook_utils.covariance.halo import halo_power_response
    from cosmolike_notebook_utils.covariance.forecast import _json_array
    from cosmolike_notebook_utils.covariance.accuracy import covariance_accuracy

    settings = configuration(
        accuracy_boost=args.accuracy_boost,
        integration_accuracy=args.integration_accuracy,
        gaussian={"nonlimber": False, "ia": "none"},
    )
    # The existing response helper resolves its derivative step and mass
    # rule from the global controls. Do not silently ignore a new survey
    # override if these defaults are changed in a later project release.
    helper = covariance_accuracy(accuracy_boost=args.accuracy_boost,
                                 integration_accuracy=args.integration_accuracy)
    for key in ("response_step", "halo_mass_nquad"):
        if settings[key] != helper[key]:
            raise ValueError(f"The response helper must honor the survey {key}")
    if settings["cosmology"] != inputs["cosmology"]:
        raise ValueError("CoCoA cosmology changed since input export")
    if sha256(project / settings["source_file"]) != inputs["source_nz_sha256"]:
        raise ValueError("CoCoA source distribution changed since export")
    if settings["area_deg2"] != inputs["area_deg2"]:
        raise ValueError("CoCoA and TJPCov survey areas differ")
    tables = initialize(interface=ci, settings=settings)
    original = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    np.testing.assert_array_equal(tables["z_2D"], original["z"])
    np.testing.assert_array_equal(10.0**tables["log10k_2D"], original["k_h_mpc"])
    for name in ("linear", "nonlinear", "linear_cb"):
        key = f"lnP_{name}"
        power = np.exp(tables[key].reshape(original[f"p_{name}"].shape,
                                          order="F"))
        np.testing.assert_array_equal(power, original[f"p_{name}"])

    ell = np.load(args.gaussian_run / "gaussian.npz")["effective_ell"]
    backend = ci.covariance
    snapshot = limber_spectra(
        interface=backend, ell=ell, a_edges=settings["a_edges"],
        nquad=settings["radial_nquad"], nwindow=settings["nwindow"],
        include_ia=False, include_rsd=False, linear=False,
    )
    # Foreground matter lenses the selected source population. Its kernel
    # W enters twice in each shear spectrum; a shear field has no galaxy
    # count normalization to subtract from the background response.
    geometry = snapshot["geometry"]
    a, distance, dchi = geometry[0], geometry[2], geometry[3]
    field = snapshot["nlens"] + inputs["source_bin_1based"] - 1
    window = snapshot["windows"][1, field]
    wave = (ell[None, :] + 0.5) / distance[:, None]
    response = halo_power_response(
        interface=backend, a=a, k=wave, lnm_edges=settings["lnm_edges"],
        accuracy_boost=args.accuracy_boost, mnu=settings["cosmology"]["mnu"],
        integration_accuracy=args.integration_accuracy,
    )
    # Each shear leg contributes the Fourier spin-2 transfer. Multiplying
    # the matter response by its square gives the response of a shear
    # spectrum; the later outer product supplies all four shear legs.
    spin = np.sqrt((ell-1)*ell*(ell+1)*(ell+2)) / (ell+0.5)**2
    pair = np.ascontiguousarray(np.broadcast_to(window**2, response.T.shape))
    shell = backend.covariance_ssc_shell_response(
        distance=np.ascontiguousarray(distance), signal=np.zeros(len(ell)),
        pair_window=pair, mean_window=np.zeros_like(pair),
        power_response=np.ascontiguousarray(response.T * spin[:, None]**2),
    )

    # The native survey example uses the raw spherical-cap mask power.
    # Its long-mode sum includes the monopole. TJPCov's circular-disc
    # variance is deliberately left as a distinct model choice.
    area = settings["area_deg2"] * (np.pi/180)**2
    mask = cap_mask(area, settings["mask_ell_max"])
    long_power = np.empty((len(a), len(mask)))
    for row, value in enumerate(a):
        long_power[row] = backend.covariance_power(
            a=float(value), k=(np.arange(len(mask))+0.5)/distance[row],
            linear=True,
        )
    variance = backend.covariance_ssc_mask_variance(
        mask_cl=mask, area_sr=area, distance=np.ascontiguousarray(distance),
        power=long_power,
    )
    covariance = backend.covariance_project(
        left=shell, right=shell, weight=np.ascontiguousarray(dchi*variance),
    )
    if covariance.shape != (5, 5) or not np.all(np.isfinite(covariance)):
        raise ValueError("Invalid CoCoA SSC matrix")
    rms = np.sqrt(np.diag(covariance))
    if not np.all(np.isfinite(rms)) or np.any(rms <= 0):
        raise ValueError("Invalid CoCoA SSC diagonal")
    correlation = covariance / rms[:, None] / rms[None, :]

    import json
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output / "ssc.npz", ell=ell, ssc=covariance,
                        geometry=geometry, window=window, response=response,
                        variance=variance, shell=shell)
    record = {
        "schema": "cocoa-ssc-shear-v1", "status": "completed",
        "scope": "Native shear SSC at effective ell; no band averaging",
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "gaussian_manifest_sha256": sha256(args.gaussian_run / "manifest.json"),
        "core": revision(core), "lsst_y1": revision(project),
        "interface_sha256": sha256(ci.__file__),
        "script_sha256": sha256(__file__), "threads": threads,
        "settings": json.loads(json.dumps(settings, default=_json_array)),
        "minimum_correlation_eigenvalue": float(np.linalg.eigvalsh(
            (correlation+correlation.T)/2)[0]),
        "max_scaled_asymmetry": float(np.max(np.abs(correlation-correlation.T))),
        "timing_scope": "Accuracy only; concurrent validation, no benchmark",
        "files": {"ssc.npz": sha256(args.output / "ssc.npz")},
    }
    write_json(args.output / "manifest.json", record)
    print(f"Saved CoCoA native five-multipole SSC to {args.output}")


if __name__ == "__main__":
    main()
