"""Compute native TJPCov SSC for the existing five-point shear pilot.

Run in the TJPCov environment after the Gaussian pilot. Reuse its exact
SACC catalog and CAMB bundle. This records a native SSC component, not a
cross-code accuracy result or a band-averaged Gaussian-plus-SSC total.
TJPCov evaluates SSC at the SACC effective multipoles; its Gaussian
calculator instead averages each band. Keep that distinction explicit.
"""

import argparse
import importlib.metadata
import inspect
import platform
import signal
from pathlib import Path
from time import perf_counter

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, help="actual LSST/CAMB export")
    parser.add_argument("gaussian_run", type=Path,
                        help="completed five-band native Gaussian run")
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--integration-method", default="qag_quad",
                        choices=("qag_quad", "spline"))
    parser.add_argument("--match-power-a-range", action="store_true",
                        help="explicitly set CCL A_SPLINE_MINLOG_PK to the "
                             "supplied CAMB table's minimum scale factor")
    parser.add_argument("--a-refinement", type=int, choices=(1, 2, 4, 8, 16),
                        default=1,
                        help="nested refinement of both CCL a-grid segments")
    parser.add_argument("--k-refinement", type=int, choices=(1, 2), default=1,
                        help="multiply CCL N_K density; nodes need not nest")
    parser.add_argument("--timeout", type=int, default=600,
                        help="Unix wall-time limit in seconds (default: 600)")
    parser.add_argument("--timing", action="store_true",
                        help="record first-use compute time on a quiet machine")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    threads = require_thread_environment()
    if threads > (8 if args.timing else 6):
        parser.error("Use at most six accuracy workers or eight timing workers")
    output = args.output.resolve()
    if output.exists():
        parser.error("Choose a fresh output directory; never reuse SSC caches")

    inputs = args.inputs.resolve()
    gaussian_run = args.gaussian_run.resolve()
    manifest = load_bundle(inputs, "lsst-y1-tjpcov-inputs-v1")
    gaussian = load_bundle(gaussian_run, "tjpcov-gaussian-shear-v1")
    if (gaussian["status"] != "completed"
            or gaussian["input_manifest_sha256"]
            != sha256(inputs / "manifest.json")):
        raise ValueError("Gaussian pilot used different or incomplete inputs")
    if (manifest["cosmology"]["mnu"] != 0
            or manifest["gaussian"]["ia"] != "none"
            or manifest["gaussian"]["nonlimber"]
            or manifest["rsd"] or manifest["magnification"]):
        raise ValueError("Need massless neutrinos, Limber and zero IA/RSD/mag")

    output.mkdir(parents=True)
    write_json(output / "run_settings.json", {
        "timeout_seconds": args.timeout, "threads": threads,
        "integration_method": args.integration_method,
        "match_power_a_range": args.match_power_a_range,
        "a_refinement": args.a_refinement, "k_refinement": args.k_refinement,
        "input_manifest_sha256": sha256(inputs / "manifest.json"),
        "gaussian_manifest_sha256": sha256(gaussian_run / "manifest.json"),
    })
    # Use the OS default action so the deadline also stops compiled code,
    # which may delay a Python signal handler. No output files are removed;
    # the supervising runner can identify SIGALRM and preserve its log.
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)

    import numpy as np
    import pyccl as ccl
    import sacc
    import tjpcov
    from tjpcov.covariance_fourier_ssc_fsky import FourierSSCHaloModelFsky

    # Identify the installed code, not only a nearby Git checkout.
    # Native TJPCov caches SSC by tracer names, so every calculation below
    # also gets a new output directory with no existing covariance files.
    package = Path(tjpcov.__file__).resolve().parent
    source_hashes = {}
    for name in ("covariance_fourier_ssc.py", "covariance_fourier_ssc_fsky.py",
                 "covariance_builder.py", "covariance_io.py"):
        installed = sha256(package / name)
        if installed != sha256(args.tjpcov / "tjpcov" / name):
            raise ValueError(f"Installed TJPCov differs from checkout: {name}")
        source_hashes[name] = installed
    if sha256(ccl._ccllib.__file__) != gaussian["ccl_binary_sha256"]:
        raise ValueError("CCL binary changed since the Gaussian input export")

    setup_started = perf_counter()
    data = np.load(inputs / "inputs.npz", allow_pickle=False)
    cosmos = manifest["cosmology"]
    h = cosmos["H0"] / 100
    a = 1 / (1 + data["z"][::-1])
    k = data["k_h_mpc"] * h
    if (manifest["units"]["k"] != "h/Mpc"
            or manifest["units"]["power"] != "(Mpc/h)^3"
            or manifest["units"]["power_axes"] != ["redshift", "wavenumber"]
            or np.any(np.diff(a) <= 0) or a[-1] != 1
            or np.any(np.diff(k) <= 0) or k[0] <= 0):
        raise ValueError("Unexpected matter-power units, axes or grid")

    # This is the same table conversion as in the Gaussian pilot:
    # CCL uses increasing a, k in Mpc^-1 and P in Mpc^3. Its distance
    # and growth calculations remain native CCL, rather than supplied CAMB.
    powers = {}
    for name in ("linear", "nonlinear"):
        values = data[f"p_{name}"]
        if (values.shape != (a.size, k.size)
                or not np.all(np.isfinite(values)) or np.any(values <= 0)):
            raise ValueError(f"Malformed {name} CAMB power table")
        powers[name] = {"a": a, "k": k,
                        "delta_matter:delta_matter": values[::-1] / h**3}

    # This documented CCL control must be set before creating the
    # cosmology. Only the low-a endpoint changes; native grid densities
    # and the logarithmic-to-linear transition remain unchanged.
    original_min_a = ccl.spline_params.A_SPLINE_MINLOG_PK
    if args.match_power_a_range:
        if a[0] >= ccl.spline_params.A_SPLINE_MIN_PK:
            raise ValueError("CAMB does not span CCL's low-a grid segment")
        ccl.spline_params.A_SPLINE_MINLOG_PK = float(a[0])
        print(f"Explicit CCL minimum-a control: {original_min_a} -> {a[0]}",
              flush=True)

    # Increase intervals rather than points on both pieces of the a
    # grid. With r=2, one new node sits between each two old nodes.
    # N_K is instead a points-per-decade control: doubling it increases
    # k resolution but does not in general retain the old interior nodes.
    base_response_a = ccl.get_pk_spline_a()
    base_response_lk = ccl.get_pk_spline_lk()
    original_nk_density = ccl.spline_params.N_K
    for name in ("A_SPLINE_NA_PK", "A_SPLINE_NLOG_PK"):
        old_count = getattr(ccl.spline_params, name)
        setattr(ccl.spline_params, name,
                args.a_refinement * (old_count - 1) + 1)
    ccl.spline_params.N_K = args.k_refinement * original_nk_density
    cosmo = ccl.CosmologyCalculator(
        Omega_c=cosmos["omegam"] - cosmos["omegab"],
        Omega_b=cosmos["omegab"], h=h, n_s=cosmos["ns"],
        A_s=cosmos["As_1e9"] * 1e-9, m_nu=0,
        w0=cosmos["w"], wa=cosmos["w0pwa"] - cosmos["w"],
        pk_linear=powers["linear"], pk_nonlin=powers["nonlinear"],
    )

    # TJPCov builds its response on all default CCL scale-factor nodes.
    # Only the background-variance grid is truncated to the source range.
    # CCL's response call forbids time extrapolation by default: detect a
    # too-short CAMB input table before starting that expensive calculation.
    native_a = cosmo.get_pk_spline_a()
    native_lk = cosmo.get_pk_spline_lk()
    np.testing.assert_allclose(native_a[::args.a_refinement], base_response_a,
                               rtol=5e-14, atol=0)
    if native_a[0] < a[0] or native_a[-1] > a[-1]:
        raise ValueError(
            f"Native SSC response needs a=[{native_a[0]}, {native_a[-1]}], "
            f"but CAMB covers [{a[0]}, {a[-1]}]. Review the input support "
            "before extending it; this script does not extrapolate silently."
        )

    catalog = sacc.Sacc.load_fits(str(gaussian_run / "source_bins.fits"))
    source = gaussian["source_name"]
    ell, _ = catalog.get_ell_cl("cl_ee", source, source)
    saved = np.load(gaussian_run / "gaussian.npz", allow_pickle=False)
    if (len(catalog.tracers) != 1 or ell.shape != (5,)
            or not np.array_equal(ell, saved["effective_ell"])):
        raise ValueError("Expected the same five-point, one-source SACC pilot")
    configuration = {
        "tjpcov": {
            "cosmo": cosmo, "sacc_file": catalog, "IA": None,
            "use_mpi": False, "outdir": str(output / "native"),
            "fsky": gaussian["area_sr"] / (4 * np.pi),
            f"Ngal_{source}": manifest["source_density_arcmin2"],
            f"sigma_e_{source}": manifest["sigma_e_component"],
        },
        "SSC": {"integration_method": args.integration_method},
    }
    native = FourierSSCHaloModelFsky(configuration)
    print(f"Computing five-point SSC with {args.integration_method}", flush=True)
    setup_seconds = perf_counter() - setup_started
    started = perf_counter()
    matrix = native.get_covariance_block(
        (source, source), (source, source), include_b_modes=False)
    elapsed = perf_counter() - started
    if matrix.shape != (5, 5) or not np.all(np.isfinite(matrix)):
        raise ValueError("Native SSC returned an invalid five-point matrix")

    # Retain the unaltered matrix, including any negative eigenvalues.
    # Positivity of an SSC component is a diagnostic; this is not a total
    # covariance and no G/cNG terms or scale cuts are inferred here.
    diagonal = np.diag(matrix)
    positive_diagonal = bool(np.all(diagonal > 0))
    peak = float(np.max(np.abs(matrix)))
    asymmetry = float(np.max(np.abs(matrix - matrix.T)))
    symmetric = asymmetry <= 1e-10 * peak
    smallest = None
    if positive_diagonal and symmetric:
        rms = np.sqrt(diagonal)
        correlation = matrix / rms[:, None] / rms[None, :]
        smallest = float(np.linalg.eigvalsh(correlation)[0])

    # Save the disc window on the same a nodes chosen inside TJPCov.
    # This second public CCL call is a diagnostic export, outside the
    # recorded native block duration. It is not a reused cached result.
    z_max = float(np.max(catalog.tracers[source].z))
    keep = 1 / native_a < z_max + 1
    keep[np.sum(~keep) - 1] = True
    window_a = native_a[keep]
    variance = ccl.sigma2_B_disc(
        cosmo, a_arr=window_a, fsky=configuration["tjpcov"]["fsky"])
    if not np.all(np.isfinite(variance)) or np.any(variance < 0):
        raise ValueError("Invalid circular-disc background variance")
    np.savez_compressed(output / "ssc.npz", ell=ell, ssc=matrix,
                        window_a=window_a, disc_variance=variance,
                        native_response_a=native_a, native_response_lk=native_lk,
                        base_response_a=base_response_a,
                        base_response_lk=base_response_lk)
    record = {
        "schema": "tjpcov-ssc-shear-v1", "status": "completed",
        "scope": "native SSC at five effective multipoles; not band averaged",
        "input_manifest_sha256": sha256(inputs / "manifest.json"),
        "gaussian_manifest_sha256": sha256(gaussian_run / "manifest.json"),
        "input_manifest": manifest, "source_name": source,
        "area_sr": gaussian["area_sr"], "shape": [5, 5],
        "effective_ell": ell.tolist(),
        "integration_method": args.integration_method,
        "a_domain_control": {
            "match_power_a_range": args.match_power_a_range,
            "a_refinement": args.a_refinement,
            "original_A_SPLINE_MINLOG_PK": original_min_a,
            "A_SPLINE_MINLOG_PK": ccl.spline_params.A_SPLINE_MINLOG_PK,
            "A_SPLINE_MIN_PK": ccl.spline_params.A_SPLINE_MIN_PK,
            "A_SPLINE_NA_PK": ccl.spline_params.A_SPLINE_NA_PK,
            "A_SPLINE_NLOG_PK": ccl.spline_params.A_SPLINE_NLOG_PK,
            "input_a_min": float(a[0]),
        },
        "k_sampling": {
            "k_refinement": args.k_refinement,
            "original_N_K": original_nk_density,
            "N_K": ccl.spline_params.N_K,
            "K_MIN_Mpc_inverse": ccl.spline_params.K_MIN,
            "K_MAX_Mpc_inverse": ccl.spline_params.K_MAX,
            "base_count": int(base_response_lk.size),
            "refined_count": int(native_lk.size),
            "nested_nodes_claimed": False,
        },
        "halo_model": {
            "mass_definition": "200m", "mass_function": "Tinker08",
            "bias": "Tinker10", "concentration": "Duffy08",
            "profile": "analytic NFW", "number_count_legs": 0,
            "response": "linear-power slope and amplitude plus normalized I12",
            "hmc_defaults": {
                name: inspect.signature(ccl.halos.HMCalculator.__init__).parameters[
                    name].default
                for name in ("log10M_min", "log10M_max", "nM",
                             "integration_method_M")},
        },
        "window": "CCL circular-disc transverse linear-power integral",
        "disc_variance_units": "Mpc (k dk times a three-dimensional P)",
        "background": "native CCL; only CAMB power is supplied",
        "native_block_seconds_diagnostic_only": elapsed,
        "timing": ({"setup_seconds": setup_seconds,
                    "construction_seconds": elapsed,
                    "scope": "First native SSC block, including halo response, "
                             "survey variance and projection; no disk-cache hit"}
                   if args.timing else None),
        "threads": threads, "timeout_seconds": args.timeout,
        "positive_diagonal": positive_diagonal,
        "symmetric_within_relative_1e-10": symmetric,
        "minimum_ssc_correlation_eigenvalue": smallest,
        "max_asymmetry_over_peak": asymmetry / peak if peak else 0.0,
        "tjpcov": revision(args.tjpcov), "installed_source": source_hashes,
        "versions": {name: importlib.metadata.version(name)
                     for name in ("tjpcov", "pyccl", "sacc", "numpy", "scipy")},
        "python": platform.python_version(), "platform": platform.platform(),
        "ccl_binary_sha256": sha256(ccl._ccllib.__file__),
        "ccl_ssc_source_sha256": sha256(inspect.getfile(
            ccl.halos.halomod_Tk3D_SSC_linear_bias)),
        "script_sha256": sha256(__file__),
        "files": {"ssc.npz": sha256(output / "ssc.npz")},
    }
    write_json(output / "manifest.json", record)
    signal.alarm(0)
    print(f"Saved native SSC at five effective multipoles in {output}")
    print("Cross-code comparison and numerical convergence remain separate checks")


if __name__ == "__main__":
    main()
