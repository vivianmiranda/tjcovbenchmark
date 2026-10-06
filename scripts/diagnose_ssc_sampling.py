"""Separate TJPCov SSC response-a and disc-variance-a sampling effects.

Use two completed native SSC runs differing only in their a-grid density.
The four calculations below use TJPCov's selected public CCL halo model,
its actual SACC tracer, shared CAMB power and effective multipoles. They
are lower-level CCL diagnostics, not modified native TJPCov calculations.
CAMB input k/P units are h/Mpc and (Mpc/h)^3; CCL uses Mpc^-1 and Mpc^3.
The transverse disc variance has units Mpc. CCL handles every projection
distance and shear-spin factor; this script inserts no extra factors.
"""

import argparse
import inspect
import signal
from pathlib import Path
from time import perf_counter

from common import load_bundle, require_thread_environment, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, help="verified LSST/CAMB bundle")
    parser.add_argument("gaussian_run", type=Path, help="matching SACC pilot")
    parser.add_argument("coarse_run", type=Path, help="completed native SSC")
    parser.add_argument("fine_run", type=Path, help="a-only refined native SSC")
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600,
                        help="Unix wall-time deadline in seconds")
    args = parser.parse_args()
    threads = require_thread_environment()
    if args.timeout <= 0 or threads > 6:
        parser.error("Need a positive deadline and at most six OMP threads")
    if args.output.exists():
        parser.error("Choose a fresh output directory")

    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    gaussian = load_bundle(args.gaussian_run, "tjpcov-gaussian-shear-v1")
    runs = [load_bundle(path, "tjpcov-ssc-shear-v1")
            for path in (args.coarse_run, args.fine_run)]
    for record in runs:
        if (record["status"] != "completed"
                or record["input_manifest_sha256"]
                != sha256(args.inputs / "manifest.json")
                or record["gaussian_manifest_sha256"]
                != sha256(args.gaussian_run / "manifest.json")):
            raise ValueError("SSC runs must use these completed input bundles")
    for key in ("integration_method", "source_name", "area_sr", "halo_model",
                "ccl_binary_sha256", "ccl_ssc_source_sha256", "installed_source",
                "versions", "effective_ell"):
        if runs[0][key] != runs[1][key]:
            raise ValueError(f"Native runs differ in more than a sampling: {key}")
    for key in ("A_SPLINE_MINLOG_PK", "A_SPLINE_MIN_PK", "input_a_min"):
        if runs[0]["a_domain_control"][key] != runs[1]["a_domain_control"][key]:
            raise ValueError(f"Native scale-factor endpoints differ: {key}")
    if (inputs["cosmology"]["mnu"] != 0
            or inputs["gaussian"]["ia"] != "none"
            or inputs["gaussian"]["nonlimber"]
            or inputs["rsd"] or inputs["magnification"]):
        raise ValueError("This diagnostic is for massless, Limber, zero-IA shear")

    args.output.mkdir(parents=True)
    report = {
        "schema": "tjpcov-ssc-sampling-diagnostic-v1", "status": "started",
        "timeout_seconds": args.timeout, "threads": threads,
        "inputs": {str(path.resolve()): sha256(path / "manifest.json")
                   for path in (args.inputs, args.gaussian_run,
                                args.coarse_run, args.fine_run)},
        "scope": "public CCL 2x2 a-sampling diagnostic at effective ell",
        "script_sha256": sha256(__file__), "cases": {},
    }
    write_json(args.output / "report.json", report)
    # The OS action also stops a long compiled call. Partial records and
    # matrices survive so the supervising runner can inspect a timeout.
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)

    import numpy as np
    import pyccl as ccl
    import sacc
    import tjpcov
    from tjpcov.covariance_fourier_ssc_fsky import FourierSSCHaloModelFsky

    package = Path(tjpcov.__file__).resolve().parent
    for name, expected in runs[0]["installed_source"].items():
        if (sha256(package / name) != expected
                or sha256(args.tjpcov / "tjpcov" / name) != expected):
            raise ValueError(f"TJPCov source differs from saved SSC: {name}")
    if (sha256(ccl._ccllib.__file__) != runs[0]["ccl_binary_sha256"]
            or sha256(inspect.getfile(ccl.halos.halomod_Tk3D_SSC_linear_bias))
            != runs[0]["ccl_ssc_source_sha256"]):
        raise ValueError("CCL changed since the native SSC calculations")

    saved = [np.load(path / "ssc.npz", allow_pickle=False)
             for path in (args.coarse_run, args.fine_run)]
    response_a = [entry["native_response_a"] for entry in saved]
    window_a = [entry["window_a"] for entry in saved]
    if response_a[1].size <= response_a[0].size:
        raise ValueError("The fine run must have more response-a samples")

    # Fix the cosmology's internal controls to the coarse native run.
    # Only explicit response/window arrays will change between projections.
    controls = runs[0]["a_domain_control"]
    for name in ("A_SPLINE_MINLOG_PK", "A_SPLINE_MIN_PK",
                 "A_SPLINE_NA_PK", "A_SPLINE_NLOG_PK"):
        setattr(ccl.spline_params, name, controls[name])
    lk = ccl.get_pk_spline_lk()
    for entry, record in zip(saved, runs):
        # The first saved baseline predates the explicit k-grid export.
        # Its unchanged CCL default is tested by the full baseline gate.
        if record.get("k_sampling", {}).get("k_refinement", 1) != 1:
            raise ValueError("Use an a-only refinement, not a k refinement")
        if "native_response_lk" in entry:
            np.testing.assert_array_equal(lk, entry["native_response_lk"])
    np.testing.assert_array_equal(ccl.get_pk_spline_a(), response_a[0])

    data = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    cosmos = inputs["cosmology"]
    h = cosmos["H0"] / 100
    a = 1 / (1 + data["z"][::-1])
    if (inputs["units"]["k"] != "h/Mpc"
            or inputs["units"]["power"] != "(Mpc/h)^3"
            or inputs["units"]["power_axes"] != ["redshift", "wavenumber"]
            or response_a[0][0] < a[0] or response_a[1][0] < a[0]):
        raise ValueError("Unexpected CAMB units or insufficient time support")
    powers = {
        name: {"a": a, "k": data["k_h_mpc"] * h,
               "delta_matter:delta_matter": data[f"p_{name}"][::-1] / h**3}
        for name in ("linear", "nonlinear")}
    cosmo = ccl.CosmologyCalculator(
        Omega_c=cosmos["omegam"] - cosmos["omegab"],
        Omega_b=cosmos["omegab"], h=h, n_s=cosmos["ns"],
        A_s=cosmos["As_1e9"] * 1e-9, m_nu=0,
        w0=cosmos["w"], wa=cosmos["w0pwa"] - cosmos["w"],
        pk_linear=powers["linear"], pk_nonlin=powers["nonlinear"])

    catalog = sacc.Sacc.load_fits(str(args.gaussian_run / "source_bins.fits"))
    source = gaussian["source_name"]
    if (len(catalog.tracers) != 1
            or catalog.tracers[source].quantity != "galaxy_shear"
            or "lens" in source):
        raise ValueError("Expected one source tracer and no number-count legs")
    fsky = gaussian["area_sr"] / (4 * np.pi)
    native = FourierSSCHaloModelFsky({"tjpcov": {
        "cosmo": cosmo, "sacc_file": catalog, "IA": None, "use_mpi": False,
        "outdir": str(args.output / "tracer_setup"), "fsky": fsky,
        f"Ngal_{source}": inputs["source_density_arcmin2"],
        f"sigma_e_{source}": inputs["sigma_e_component"]}})
    tracers, _ = native.get_tracer_info()
    tracer = tracers[source]
    ell = native.get_ell_eff()
    for entry in saved:
        np.testing.assert_array_equal(ell, entry["ell"])
    z_max = float(np.max(catalog.tracers[source].z))
    for nodes, selected in zip(response_a, window_a):
        keep = 1 / nodes < z_max + 1
        keep[np.sum(~keep) - 1] = True
        np.testing.assert_array_equal(nodes[keep], selected)

    # These are precisely TJPCov's native halo choices and default mass
    # integrator. Four source legs have bias 1 and no local-count subtraction.
    mass_def = ccl.halos.MassDef200m
    hmc = ccl.halos.HMCalculator(
        mass_function=ccl.halos.MassFuncTinker08(mass_def=mass_def),
        halo_bias=ccl.halos.HaloBiasTinker10(mass_def=mass_def),
        mass_def=mass_def)
    profile = ccl.halos.HaloProfileNFW(
        mass_def=mass_def, fourier_analytic=True,
        concentration=ccl.halos.ConcentrationDuffy08(mass_def=mass_def))
    responses = {}
    variances = {}
    arrays = {"ell": ell, "lk": lk, "coarse_response_a": response_a[0],
              "fine_response_a": response_a[1], "coarse_window_a": window_a[0],
              "fine_window_a": window_a[1]}
    baseline_rms = np.sqrt(np.diag(saved[0]["ssc"]))
    if np.any(baseline_rms <= 0) or not np.all(np.isfinite(baseline_rms)):
        raise ValueError("Native baseline must have a finite positive diagonal")

    def difference(matrix, reference):
        delta = matrix - reference
        return {
            "max_abs_difference_over_coarse_diagonal_rms": float(np.max(
                np.abs(delta / baseline_rms[:, None] / baseline_rms[None, :]))),
            "max_fractional_diagonal_change": float(np.max(np.abs(
                np.diag(delta) / np.diag(reference))))}

    # First reproduce the native baseline. Only after that gate succeeds
    # can either independent refinement be interpreted as a sampling effect.
    for label, response_index, window_index in (
            ("coarse_coarse", 0, 0), ("fine_coarse", 1, 0),
            ("coarse_fine", 0, 1), ("fine_fine", 1, 1)):
        print(f"Computing {label}: response grid, then variance grid", flush=True)
        started = perf_counter()
        if response_index not in responses:
            responses[response_index] = ccl.halos.halomod_Tk3D_SSC_linear_bias(
                cosmo=cosmo, hmc=hmc, prof=profile,
                a_arr=response_a[response_index], lk_arr=lk)
        if window_index not in variances:
            variances[window_index] = ccl.sigma2_B_disc(
                cosmo, a_arr=window_a[window_index], fsky=fsky)
        matrix = ccl.covariances.angular_cl_cov_SSC(
            cosmo, tracer1=tracer, tracer2=tracer, tracer3=tracer,
            tracer4=tracer, ell=ell, t_of_kk_a=responses[response_index],
            sigma2_B=(window_a[window_index], variances[window_index]),
            integration_method=runs[0]["integration_method"])
        if matrix.shape != (5, 5) or not np.all(np.isfinite(matrix)):
            raise ValueError(f"Invalid projected SSC: {label}")
        arrays[label] = matrix
        arrays[f"disc_variance_{window_index}"] = variances[window_index]
        report["cases"][label] = {
            "response_a_count": int(response_a[response_index].size),
            "window_a_count": int(window_a[window_index].size),
            "elapsed_seconds_diagnostic_only": perf_counter() - started,
            "versus_coarse_native": difference(matrix, saved[0]["ssc"])}
        np.savez_compressed(args.output / "sampling.npz", **arrays)
        write_json(args.output / "report.json", report)
        if label == "coarse_coarse":
            error = report["cases"][label]["versus_coarse_native"][
                "max_abs_difference_over_coarse_diagonal_rms"]
            if error > 1e-8:
                report["status"] = "baseline_reproduction_failed"
                write_json(args.output / "report.json", report)
                raise ValueError(f"Native baseline gate failed: {error:.8g}")

    report.update({
        "status": "completed", "integration_method": runs[0]["integration_method"],
        "halo_model": runs[0]["halo_model"], "disc_variance_units": "Mpc",
        "fixed_cosmology_a_controls": controls,
        "baseline_reproduction_tolerance": 1e-8,
        "fine_fine_versus_fine_native": difference(
            arrays["fine_fine"], saved[1]["ssc"]),
        "ccl_binary_sha256": runs[0]["ccl_binary_sha256"],
        "files": {"sampling.npz": sha256(args.output / "sampling.npz")},
    })
    write_json(args.output / "report.json", report)
    signal.alarm(0)
    print(f"Saved four SSC sampling cases to {args.output}")


if __name__ == "__main__":
    main()
