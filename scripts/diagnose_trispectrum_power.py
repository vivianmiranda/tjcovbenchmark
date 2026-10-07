"""Isolate power-table interpolation in the unequal-k tree trispectrum.

This changes supplied tables only in this diagnostic process. A natural
cubic interpolant of the SAME installed log-power samples fills denser k
grids. The input bundle records how CAMB's native samples were prepared.
The actual CoCoA reader remains linear in log power. Halo moments, z nodes
and angular quadrature stay fixed. Direct cubic evaluation is a smooth
input control, not an independent trispectrum code or new CAMB solution.
"""

import argparse
import signal
import sys
from pathlib import Path

from common import load_bundle, require_thread_environment, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("cocoa_run", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--factors", type=int, nargs="+", default=[1, 2],
                        choices=(1, 2, 4, 8, 16),
                        help="Refine the installed grid, not raw CAMB; "
                             "default 11993 then 23985 nodes")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if (args.factors[0] != 1
            or any(b <= a for a, b in zip(args.factors, args.factors[1:]))
            or args.timeout <= 0):
        parser.error("Use increasing factors starting at 1 and positive timeout")
    threads = require_thread_environment()
    if threads > 6 or args.output.exists():
        parser.error("Use at most six threads and a new directory")
    meta = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    run = load_bundle(args.cocoa_run, "cocoa-trispectrum-v1")
    if run["input_manifest_sha256"] != sha256(args.inputs / "manifest.json"):
        raise ValueError("Different original inputs")
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)

    import numpy as np
    from scipy.interpolate import CubicSpline

    cocoa = args.cocoa.resolve()
    core = cocoa / "external_modules/code/cosmolike_core"
    sys.path[:0] = [str(core), str(cocoa / "projects/lsst_y1/covariance")]
    import cosmolike_lsst_y1_interface as ci
    if sha256(ci.__file__) != run["interface_sha256"]:
        raise ValueError("CoCoA library changed since the native export")
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.power import refine_power_tables
    settings = configuration(gaussian={"nonlimber": False, "ia": "none"},
                             integration_accuracy=run["integration_accuracy"])
    tables = initialize(interface=ci, settings=settings)
    original = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    grid = np.load(args.cocoa_run / "trispectrum.npz", allow_pickle=False)
    np.testing.assert_array_equal(10**tables["log10k_2D"], original["k_h_mpc"])
    x, z_nodes = tables["log10k_2D"], tables["z_2D"]
    np.testing.assert_array_equal(z_nodes, original["z"])
    for name in ("linear", "nonlinear", "linear_cb"):
        values = tables[f"lnP_{name}"].reshape((len(z_nodes), len(x)), order="F")
        np.testing.assert_array_equal(np.exp(values), original[f"p_{name}"])
    first, second = grid["first"], grid["second"]
    k = grid["k"]
    length = 2997.92458
    pairs = np.array([k[first], k[second]])

    def assemble(power, internal, index):
        pk = power[[first, second]]
        angular = ci.covariance.covariance_tree_averages(
            k=pairs, pk=pk, corner=grid["corner"], weight=grid["weight"],
            ps=np.ascontiguousarray(internal))
        return ci.covariance.covariance_halo_trispectrum(
            pk=pk, tree=angular, i11=np.ascontiguousarray(grid["i11"][index]),
            moments=np.ascontiguousarray(grid["moments"][index]))

    # Direct smooth control: interpolate in redshift exactly as CoCoA does,
    # while retaining cubic log-k interpolation between the original nodes.
    # Below/above the supplied k range retain the original edge power law.
    smooth = []
    smooth_power = []
    for index, z in enumerate(grid["redshift"]):
        right = np.searchsorted(z_nodes, z, side="right")
        left = right-1
        f = (z-z_nodes[left])/(z_nodes[right]-z_nodes[left])
        values = ((1-f)*np.log(original["p_linear"][left])
                  +f*np.log(original["p_linear"][right]))
        curve = CubicSpline(x, values, bc_type="natural")

        def evaluate(wave):
            q = np.log10(wave)
            logp = curve(q)
            for edge, neighbor, mask in ((0, 1, q < x[0]),
                                          (-1, -2, q > x[-1])):
                slope = (values[neighbor]-values[edge])/(x[neighbor]-x[edge])
                logp[mask] = values[edge]+slope*(q[mask]-x[edge])
            return np.exp(logp)

        power, internal = evaluate(k), evaluate(grid["internal_k"])
        smooth.append(assemble(power, internal, index))
        smooth_power.append(power)

    factors = np.array(args.factors)
    cases, powers = [], []
    for factor in factors:
        if factor > 1:
            # Use the production natural-spline helper for all three power
            # tables. Only their k sampling changes; moments stay frozen.
            new_tables = refine_power_tables(tables, refinement=int(factor))
            np.testing.assert_array_equal(new_tables["log10k_2D"][::factor], x)
            ci.set_cosmology(omegam=meta["cosmology"]["omegam"],
                             omegab=meta["cosmology"]["omegab"],
                             H0=meta["cosmology"]["H0"], **new_tables)
        rows, power_rows = [], []
        for index, z in enumerate(grid["redshift"]):
            a = float(1/(1+z))
            power = ci.covariance.covariance_power(
                a=a, k=k*length, linear=True)*length**3
            internal = ci.covariance.covariance_power(
                a=a, k=grid["internal_k"]*length, linear=True)*length**3
            rows.append(assemble(power, internal, index))
            power_rows.append(power)
        cases.append(rows)
        powers.append(power_rows)
        print(f"Power grid factor {factor}: {factor*(len(x)-1)+1} nodes", flush=True)
    cases, smooth = np.array(cases), np.array(smooth)
    np.testing.assert_allclose(cases[0], grid["terms"], rtol=1e-12, atol=0)
    if not np.isfinite(cases).all() or not np.isfinite(smooth).all():
        raise ValueError("Non-finite power interpolation diagnostic")
    # Equal pairs can probe beyond the input range near opposite angles.
    # Densifying the end intervals changes extrapolation secants, whereas
    # the smooth control retains the original slopes. Isolate pairs whose
    # every internal sample stays inside the supplied domain as well.
    interior = np.all((grid["internal_k"] >= 10**x[0])
                      & (grid["internal_k"] <= 10**x[-1]), axis=1)
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output / "diagnostic.npz", factors=factors,
                        node_counts=factors*(len(x)-1)+1, terms=cases,
                        smooth=smooth, k=k, redshift=grid["redshift"],
                        first=first, second=second, power=np.array(powers),
                        smooth_power=smooth_power, interior_pairs=interior)
    write_json(args.output / "report.json", {
        "schema": "trispectrum-power-diagnostic-v1", "status": "completed",
        "scope": "Same installed power samples; natural cubic fill, "
                 "production linear lookup; fixed halo moments",
        "spline_boundary": "natural",
        "refinement_base": "installed input-bundle grid",
        "power_preparation": meta.get("power_preparation"),
        "helper_sha256": sha256(
            core / "cosmolike_notebook_utils/covariance/power.py"),
        "timeout_seconds": args.timeout,
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "cocoa_manifest_sha256": sha256(args.cocoa_run / "manifest.json"),
        "factors": factors.tolist(), "node_counts": (factors*(len(x)-1)+1).tolist(),
        "max_fractional_vs_smooth_by_halo_order": np.max(
            np.abs(cases/smooth-1), axis=(1, 3)).tolist(),
        "interior_pair_count": int(interior.sum()),
        "interior_max_fractional_vs_smooth_by_halo_order": np.max(
            np.abs(cases[..., interior]/smooth[..., interior]-1),
            axis=(1, 3)).tolist(),
        "boundary_scope": "Dense edge secants change; smooth control keeps original extrapolation. Interior metric excludes those samples.",
        "script_sha256": sha256(__file__), "interface_sha256": sha256(ci.__file__),
        "files": {"diagnostic.npz": sha256(args.output / "diagnostic.npz")},
    })


if __name__ == "__main__":
    main()
