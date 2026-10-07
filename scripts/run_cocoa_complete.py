"""Compute CoCoA's complete LSST source-3 Fourier covariance pilot.

The 100 bands have integer edges 15,45,...,3015 and nominal centres
30,60,...,3000. Gaussian covariance uses TJPCov's discrete ell weights,
excluding each upper edge. SSC and cNG are evaluated at the centres.
This deliberately matches the native TJPCov estimator convention; the
sum is not a uniformly band-averaged covariance. All spectra and halo
quantities remain CoCoA's native predictions from the shared CAMB inputs.

Run in the CoCoA environment. The C production interface performs every
physical integral and covariance contraction. A one-source radial cNG
contraction uses covariance_project to avoid allocating the unused three
probe slots required by covariance_project_connected. No production code,
quadrature rule, matrix entry or model prescription is replaced.
"""

import argparse
import json
import signal
import sys
from pathlib import Path
from time import perf_counter

from common import load_bundle, require_thread_environment, revision
from common import sha256, write_json


def band_operators(ell, edges):
    """Return TJPCov's ell-weighted discrete average, upper edge excluded."""
    import numpy as np

    weights = np.zeros((len(edges)-1, len(ell)))
    for row, (lower, upper) in enumerate(zip(edges[:-1], edges[1:])):
        selected = (ell >= lower) & (ell < upper)
        weights[row, selected] = ell[selected] / np.sum(ell[selected])
    np.testing.assert_allclose(weights.sum(axis=1), 1.0, rtol=0, atol=1.e-14)
    return weights


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--integration-accuracy", type=int, choices=range(3),
                        default=0)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--timing", action="store_true",
                        help="Publish first-use stage times; run only when quiet")
    parser.add_argument("--export-projection-indices", default="0,1,3,7,15,31,63,99",
                        help="Comma-separated band indices for the small matter archive")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a fresh output directory")
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    threads = require_thread_environment()
    if threads > 8:
        parser.error("The laptop comparison permits at most eight OpenMP threads")
    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    if inputs["source_bin_1based"] != 3:
        parser.error("This pilot requires LSST source bin 3")
    if inputs.get("rsd") or inputs.get("magnification"):
        parser.error("This pilot requires zero RSD and magnification")
    if inputs["power_preparation"]["installed_k_nodes"] != 11993:
        parser.error("Use the current 11,993-node shared-power archive")

    import numpy as np

    try:
        selected = np.array([int(value) for value in
                             args.export_projection_indices.split(",")], dtype=int)
    except ValueError:
        parser.error("projection indices must be comma-separated integers")
    if (len(selected) == 0 or np.any(selected < 0) or np.any(selected >= 100)
            or np.any(np.diff(selected) <= 0)):
        parser.error("projection indices must increase strictly within 0..99")

    cocoa = args.cocoa.resolve()
    project = cocoa / "projects/lsst_y1"
    core = cocoa / "external_modules/code/cosmolike_core"
    sys.path[:0] = [str(project), str(project / "interface"),
                   str(project / "covariance"), str(core)]
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.forecast import _json_array
    from cosmolike_notebook_utils.covariance.gaussian import limber_spectra
    from cosmolike_notebook_utils.covariance.geometry import cap_mask
    from cosmolike_notebook_utils.covariance.survey import _matter_covariance_tables

    settings = configuration(
        accuracy_boost=1, integration_accuracy=args.integration_accuracy,
        gaussian={"nonlimber": False, "ia": "none"},
    )
    source = inputs["source_bin_1based"]-1
    expected = (
        (settings["cosmology"], inputs["cosmology"]),
        (settings["gaussian"], inputs["gaussian"]),
        (settings["area_deg2"], inputs["area_deg2"]),
        (settings["source_density_arcmin2"][source],
         inputs["source_density_arcmin2"]),
        (settings["sigma_e_component"][source], inputs["sigma_e_component"]),
        (sha256(project / settings["source_file"]), inputs["source_nz_sha256"]),
    )
    if any(actual != reference for actual, reference in expected):
        raise ValueError("CoCoA survey, Gaussian model or cosmology changed since export")

    # A new directory and running record preserve evidence if the bounded
    # numerical job times out. No pre-existing result is overwritten.
    args.output.mkdir(parents=True)
    signal.alarm(args.timeout)
    write_json(args.output / "status.json", {"status": "running"})
    setup_started = perf_counter()
    tables = initialize(interface=ci, settings=settings)
    setup_seconds = perf_counter()-setup_started
    with np.load(args.inputs / "camb_tables.npz", allow_pickle=False) as saved:
        if set(tables) != set(saved.files):
            raise ValueError("Installed CAMB table keys differ from the archive")
        for name, values in tables.items():
            np.testing.assert_array_equal(values, saved[name],
                                          err_msg=f"Installed table changed: {name}")
    with np.load(args.inputs / "inputs.npz", allow_pickle=False) as saved:
        np.testing.assert_array_equal(tables["z_2D"], saved["z"])
        np.testing.assert_array_equal(10.0**tables["log10k_2D"], saved["k_h_mpc"])
        for name in ("linear", "nonlinear", "linear_cb"):
            values = np.exp(tables[f"lnP_{name}"].reshape(
                saved[f"p_{name}"].shape, order="F"))
            np.testing.assert_array_equal(values, saved[f"p_{name}"])
    if len(tables["log10k_2D"]) != 11993:
        raise ValueError("The installed power grid is not the current baseline")

    backend = ci.covariance
    edges = np.arange(15, 3016, 30, dtype=np.int32)
    centres = np.arange(30, 3001, 30, dtype=float)
    ell = np.arange(edges[0], edges[-1]+1, dtype=float)
    operators = band_operators(ell, edges)
    area = settings["area_deg2"]*(np.pi/180)**2
    stages = {}
    construction_started = perf_counter()

    def checkpoint(name, started):
        stages[name] = perf_counter()-started
        print(f"{name}: {perf_counter()-construction_started:.2f} s", flush=True)

    tick = perf_counter()
    snapshot = limber_spectra(
        interface=backend, ell=ell, a_edges=settings["a_edges"],
        nquad=settings["radial_nquad"], nwindow=settings["nwindow"],
        include_ia=False, include_rsd=False, linear=False,
    )
    field = snapshot["nlens"]+source
    # Retain only this field in the Gaussian assembly. The native spectrum
    # snapshot still includes all crossed fields used by production.
    spectra = np.ascontiguousarray(snapshot["spectra"][:, field:field+1,
                                                       field:field+1])
    noise = np.array([inputs["sigma_e_component"]**2 /
                      (inputs["source_density_arcmin2"]*(60*180/np.pi)**2)])
    geometry = snapshot["geometry"]
    window = np.ascontiguousarray(snapshot["windows"][1, field])
    checkpoint("native_limber_spectra", tick)

    tick = perf_counter()
    arguments = dict(pairs=np.array([[0, 0]], dtype=np.int32),
                     operators=operators, ell_min=int(ell[0]), area_sr=area)
    gaussian = backend.covariance_gaussian_fourier(
        spectra=spectra, noise=noise, **arguments)
    signal_band = backend.covariance_project(
        left=np.ascontiguousarray(spectra[:, 0, 0][None, :]),
        right=operators, weight=np.ones(len(ell)),
    )[0]
    signal_centres = spectra[(centres-ell[0]).astype(int), 0, 0]
    checkpoint("gaussian_and_signal", tick)

    tick = perf_counter()
    mask = cap_mask(area, settings["mask_ell_max"])
    spin = np.sqrt((centres-1)*centres*(centres+1)*(centres+2))/(centres+0.5)**2
    # Each shear spectrum has two shear legs. A diagonal transform applies
    # that factor at the exact requested ell; no trispectrum interpolation
    # or broad-band averaging enters the point-sampled SSC/cNG components.
    transform = np.diag(spin**2)
    matter_settings = dict(settings, mnu=settings["cosmology"]["mnu"])

    def progress(completed, total):
        print(f"Matter shell {completed}/{total}: "
              f"{perf_counter()-construction_started:.2f} s", flush=True)

    matter = _matter_covariance_tables(
        interface=backend, settings=matter_settings, geometry=geometry,
        coarse_ell=centres, transform=transform, mask_nell=len(mask),
        progress=progress,
    )
    checkpoint("shared_halo_response_and_trispectrum", tick)

    tick = perf_counter()
    distance = np.ascontiguousarray(geometry[2])
    dchi = np.ascontiguousarray(geometry[3])
    pair_window = np.ascontiguousarray(np.broadcast_to(
        window**2, matter["response"].shape))
    shell = backend.covariance_ssc_shell_response(
        distance=distance, signal=signal_centres,
        pair_window=pair_window, mean_window=np.zeros_like(pair_window),
        power_response=np.ascontiguousarray(matter["response"]),
    )
    variance = backend.covariance_ssc_mask_variance(
        mask_cl=mask, area_sr=area, distance=distance,
        power=np.ascontiguousarray(matter["long_power"]),
    )
    ssc = backend.covariance_project(
        left=shell, right=shell, weight=np.ascontiguousarray(dchi*variance),
    )
    # Every source-3 covariance entry has the same W^4 radial weight.
    # The C contraction integrates all point-pair rows directly. Its input
    # is a view of the one-source table (~54 MB for the default 672 nodes),
    # avoiding a four-probe tensor (~860 MB) with fifteen unused blocks.
    nband, nradial = len(centres), len(distance)
    cng = backend.covariance_project(
        left=np.ascontiguousarray(matter["projected"].reshape(nband**2, nradial)),
        right=np.ascontiguousarray(window[None, :]**4),
        weight=np.ascontiguousarray(dchi/(area*distance**6)),
    ).reshape(nband, nband)
    total = gaussian+ssc+cng
    checkpoint("ssc_and_connected_projection", tick)
    stages["total"] = perf_counter()-construction_started

    # Preserve the scientific arrays before any acceptance diagnostic. The
    # optional matter subset is [ell_left,ell_right,radial_node]. Its k grid
    # varies with radial distance, so it is explicitly not a ready Tk3D.
    matrices = dict(gaussian=gaussian, ssc=ssc, cng=cng, total=total,
                    ell=centres, edges=edges, signal=signal_band,
                    signal_centres=signal_centres, gaussian_ell=ell,
                    spectra=spectra, noise_power=noise, operators=operators)
    np.savez_compressed(args.output / "covariance.npz", **matrices)
    point_matter = matter["projected"][np.ix_(selected, selected,
                                             np.arange(nradial))].copy()
    point_matter /= spin[selected, None, None]**2
    point_matter /= spin[None, selected, None]**2
    np.savez_compressed(
        args.output / "projection_inputs.npz", geometry=geometry, window=window,
        ell=centres[selected], k=(centres[selected, None]+0.5)/distance[None, :],
        matter_trispectrum=point_matter, response=matter["response"],
        shell=shell, variance=variance, mask_cl=mask,
    )
    rms = np.sqrt(np.diag(total))
    if np.any(~np.isfinite(rms)) or np.any(rms <= 0):
        raise ValueError("The total has invalid diagonal variances; outputs retained")
    diagnostics = {}
    for name in ("gaussian", "ssc", "cng", "total"):
        values = matrices[name]
        if values.shape != (100, 100) or not np.all(np.isfinite(values)):
            raise ValueError(f"Invalid {name} matrix; outputs retained")
        scaled = values/rms[:, None]/rms[None, :]
        diagnostics[name] = {
            "max_total_rms_scaled_asymmetry": float(np.max(np.abs(scaled-scaled.T))),
            "minimum_diagonal": float(np.min(np.diag(values))),
            "minimum_total_rms_scaled_symmetric_eigenvalue": float(
                np.linalg.eigvalsh((scaled+scaled.T)/2)[0]),
        }
    positive = diagnostics["total"]["minimum_total_rms_scaled_symmetric_eigenvalue"] > 0
    record = {
        "schema": "cocoa-complete-shear-v1",
        "status": "completed" if positive else "nonpositive-total",
        "scope": "Native CoCoA, LSST source bin 3, complete Fourier pilot",
        "estimator": {"gaussian": "Discrete ell-weighted bands, upper edge excluded",
                      "ssc": "Point evaluation at nominal band centres",
                      "cng": "Point evaluation at nominal band centres",
                      "total": "Sum of those native TJPCov estimator conventions"},
        "source_bin_1based": 3, "shape": [100, 100], "area_sr": area,
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "input_files_sha256": inputs["files"],
        "archived_power_preparation": inputs["power_preparation"],
        "installed_tables_equal_archive": True,
        "core": revision(core), "lsst_y1": revision(project),
        "interface_sha256": sha256(ci.__file__), "script_sha256": sha256(__file__),
        "helper_sha256": {str(path.relative_to(core)): sha256(path) for path in (
            core / "cosmolike_notebook_utils/covariance/survey.py",
            core / "cosmolike_notebook_utils/covariance/power.py",
            core / "cosmolike_notebook_utils/covariance/gaussian.py")},
        "threads": threads,
        "settings": json.loads(json.dumps(settings, default=_json_array)),
        "models": {"footprint": "spherical-cap mask; native monopole retained",
                   "halo": "Native CoCoA Wynn I11, finite higher moments",
                   "gaussian": "Limber, zero IA/RSD/magnification"},
        "projection_archive": {
            "indices": selected.tolist(),
            "geometry_rows": ["a", "chi", "f_K", "dchi"],
            "distance_units": "c/H0", "window_units": "(c/H0)^-1",
            "k_units": "(c/H0)^-1", "matter_trispectrum_units": "(c/H0)^9",
            "matter_axes": ["ell_left", "ell_right", "radial_node"],
            "spin_factors_removed": True,
            "note": "Physical k varies with a; not a fixed-k Tk3D table",
        },
        "projected_matter_bytes": int(matter["projected"].nbytes),
        "diagnostics": diagnostics, "total_positive": positive,
        "timing_scope": "First-use native construction" if args.timing else "Accuracy only",
        "timing": dict(setup_seconds=setup_seconds, stages_seconds=stages)
                  if args.timing else None,
        "files": {name: sha256(args.output / name)
                  for name in ("covariance.npz", "projection_inputs.npz")},
    }
    write_json(args.output / "manifest.json", record)
    write_json(args.output / "status.json", {"status": record["status"]})
    signal.alarm(0)
    if not positive:
        raise ValueError("Total covariance is not positive; no eigenvalues were repaired")
    print(f"Saved all four 100x100 matrices to {args.output}", flush=True)


if __name__ == "__main__":
    main()
