"""Isolate galaxy-bias placement with native TJPCov tracers and CCL projection.

This is a scaling diagnostic, not a physical survey covariance. The cNG
input is a positive, high-k-suppressed test trispectrum. The SSC input
uses CCL's actual linear-bias response, including number-count subtraction,
but a constant supplied background variance (1 Mpc). Those common inputs
remove unrelated window/halo differences from the bias-placement test.
"""

import argparse
import signal
from pathlib import Path

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("gaussian_run", type=Path)
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    threads = require_thread_environment()
    if threads > 6 or args.timeout <= 0:
        parser.error("Use at most six threads and a positive timeout")
    if args.output.exists():
        parser.error("Choose a fresh output directory")
    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    gaussian = load_bundle(args.gaussian_run, "tjpcov-gaussian-fields-v1")
    if (gaussian["status"] != "completed" or gaussian["input_manifest_sha256"]
            != sha256(args.inputs / "manifest.json")):
        raise ValueError("Need a completed three-field Gaussian input bundle")
    if inputs["cosmology"]["mnu"] != 0:
        raise ValueError("This diagnostic requires massless neutrinos")
    args.output.mkdir(parents=True)
    write_json(args.output / "run_settings.json", {
        "threads": threads, "timeout_seconds": args.timeout,
        "scope": "Bias placement only; no physical covariance claim",
    })
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)

    import numpy as np
    import pyccl as ccl
    import sacc
    import tjpcov
    from tjpcov.covariance_gaussian_fsky import FourierGaussianFsky

    package = Path(tjpcov.__file__).resolve().parent
    sources = ("covariance_builder.py", "covariance_fourier_cNG.py",
               "covariance_fourier_ssc.py", "covariance_gaussian_fsky.py")
    for name in sources:
        if sha256(package / name) != sha256(args.tjpcov / "tjpcov" / name):
            raise ValueError(f"Installed TJPCov differs from checkout: {name}")
    data = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    cosmos = inputs["cosmology"]
    h = cosmos["H0"]/100
    a_table = 1/(1+data["z"][::-1])
    powers = {
        name: {"a": a_table, "k": data["k_h_mpc"]*h,
               "delta_matter:delta_matter": data[f"p_{name}"][::-1]/h**3}
        for name in ("linear", "nonlinear")
    }
    cosmo = ccl.CosmologyCalculator(
        Omega_c=cosmos["omegam"]-cosmos["omegab"],
        Omega_b=cosmos["omegab"], h=h, n_s=cosmos["ns"],
        A_s=cosmos["As_1e9"]*1e-9, m_nu=0, w0=cosmos["w"],
        wa=cosmos["w0pwa"]-cosmos["w"],
        pk_linear=powers["linear"], pk_nonlin=powers["nonlinear"],
    )
    catalog = sacc.Sacc.load_fits(str(args.gaussian_run / "source_bins.fits"))
    lens = gaussian["field_names"][0]
    source = gaussian["source_name"]
    ell = np.load(args.gaussian_run / "gaussian.npz")["effective_ell"][:2]
    fsky = gaussian["area_sr"]/(4*np.pi)

    # Use the actual TJPCov builder to put b(z)=b into a galaxy tracer.
    # Only lens1 changes between calls. The source tracer should be exactly
    # unaffected; lens2 remains present in SACC but is unused in this test.
    tracers = {}
    for bias in (1, 2):
        options = {
            "cosmo": cosmo, "sacc_file": catalog, "IA": None,
            "use_mpi": False, "fsky": fsky,
            "outdir": str(args.output / f"tracers_b{bias}"),
            f"Ngal_{source}": inputs["source_density_arcmin2"],
            f"sigma_e_{source}": inputs["sigma_e_component"],
        }
        for index, name in enumerate(gaussian["field_names"][:-1]):
            options[f"bias_{name}"] = float(bias if name == lens else 1)
            options[f"Ngal_{name}"] = inputs["lens_density_arcmin2"][index]
        tracers[bias], _ = FourierGaussianFsky(
            {"tjpcov": options}).get_tracer_info()

    # Eight time nodes span the full supplied tracer catalogs. A declining
    # high-k tail keeps the diagnostic shear projection finite at chi=0;
    # a constant matter trispectrum would not have that property.
    amin = 1/(1+max(float(catalog.tracers[name].z.max())
                    for name in (lens, source)))
    if amin < a_table[0]:
        raise ValueError("CAMB does not cover the complete tracer catalogs")
    scale = np.linspace(amin, 1, 8)
    lk = np.linspace(np.log(1e-5), np.log(1e3), 24)
    envelope = 1/(1+(np.exp(lk)/0.1)**2)
    test_t = np.broadcast_to(envelope[:, None]*envelope[None, :],
                             (len(scale), len(lk), len(lk))).copy()
    unit_t = ccl.Tk3D(a_arr=scale, lk_arr=lk,
                      tkk_arr=np.log(test_t), is_logt=True)
    arrays = {"ell": ell, "a": scale, "lk": lk, "test_trispectrum": test_t}
    for bias in (1, 2):
        galaxy = tracers[bias][lens]
        options = dict(ell=ell, fsky=fsky, integration_method="qag_quad")
        # Four galaxy transfer functions alone supply b^4. The second call
        # also inserts the b^4 multiplier used by TJPCov's higher-halo sum.
        arrays[f"cng_tracers_b{bias}"] = ccl.angular_cl_cov_cNG(
            cosmo, galaxy, galaxy, t_of_kk_a=unit_t, **options)
        weighted_t = ccl.Tk3D(a_arr=scale, lk_arr=lk,
                              tkk_arr=np.log(test_t*bias**4), is_logt=True)
        arrays[f"cng_native_placement_b{bias}"] = ccl.angular_cl_cov_cNG(
            cosmo, galaxy, galaxy, t_of_kk_a=weighted_t, **options)
        shear = tracers[bias][source]
        arrays[f"cng_shear_control_b{bias}"] = ccl.angular_cl_cov_cNG(
            cosmo, shear, shear, t_of_kk_a=unit_t, **options)

    # Keep the native halo fits and SSC number-count counterterms. For
    # each fixed response, change only the projection's galaxy tracers.
    # This isolates an extra transfer-function bias independently of the
    # response's physical bias dependence and local-density subtraction.
    mass_def = ccl.halos.MassDef200m
    profile = ccl.halos.HaloProfileNFW(
        mass_def=mass_def,
        concentration=ccl.halos.ConcentrationDuffy08(mass_def=mass_def),
        fourier_analytic=True)
    calculator = ccl.halos.HMCalculator(
        mass_function=ccl.halos.MassFuncTinker08(mass_def=mass_def),
        halo_bias=ccl.halos.HaloBiasTinker10(mass_def=mass_def),
        mass_def=mass_def)
    for bias in (1, 2):
        response = ccl.halos.halomod_Tk3D_SSC_linear_bias(
            cosmo, calculator, prof=profile, a_arr=scale, lk_arr=lk,
            bias1=float(bias), bias2=float(bias),
            bias3=float(bias), bias4=float(bias),
            is_number_counts1=True, is_number_counts2=True,
            is_number_counts3=True, is_number_counts4=True)
        for label, tracer_bias in (("unit_tracers", 1),
                                   ("native_placement", bias)):
            galaxy = tracers[tracer_bias][lens]
            arrays[f"ssc_{label}_b{bias}"] = ccl.angular_cl_cov_SSC(
                cosmo, galaxy, galaxy, ell=ell, t_of_kk_a=response,
                sigma2_B=(scale, np.ones(len(scale))),
                integration_method="qag_quad")

    # Ratios between b=1 and2 alone are not a simple power law for SSC:
    # its number-count counterterm also changes. Compare the SAME response
    # with unit versus biased tracers to isolate the second bias factor.
    metrics = {}
    checks = (
        ("cng_tracer_ratio", "cng_tracers_b2", "cng_tracers_b1", 16),
        ("cng_native_ratio", "cng_native_placement_b2",
         "cng_native_placement_b1", 256),
        ("shear_control_ratio", "cng_shear_control_b2",
         "cng_shear_control_b1", 1),
        ("ssc_extra_bias_ratio", "ssc_native_placement_b2",
         "ssc_unit_tracers_b2", 16),
    )
    for name, numerator, denominator, expected in checks:
        first, second = arrays[numerator], arrays[denominator]
        if (not np.all(np.isfinite(first))
                or not np.all(np.isfinite(second)) or np.any(second == 0)):
            raise ValueError(f"Invalid matrix for {name}")
        ratio = first/second
        metrics[name] = {"expected": expected,
                         "range": [float(ratio.min()), float(ratio.max())],
                         "max_relative_error": float(np.max(
                             np.abs(ratio/expected-1)))}
    np.savez_compressed(args.output / "diagnostic.npz", **arrays)
    record = {
        "schema": "tjpcov-galaxy-bias-diagnostic-v1",
        "scope": "Public projection scaling; not a physical covariance",
        "passed": all(value["max_relative_error"] < 1e-6
                      for value in metrics.values()),
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "gaussian_manifest_sha256": sha256(args.gaussian_run / "manifest.json"),
        "metrics": metrics, "threads": threads,
        "inputs": {"cng_trispectrum": "1Mpc^9/[(1+(k1/0.1)^2)(1+(k2/0.1)^2)]; k in Mpc^-1",
                   "ssc_background_variance": "constant1Mpc",
                   "ssc_response": "native CCL linear bias + count subtraction",
                   "galaxy_bias": [1, 2], "lens": lens, "source": source},
        "tjpcov": revision(args.tjpcov),
        "source_sha256": {name: sha256(package / name) for name in sources},
        "ccl_binary_sha256": sha256(ccl._ccllib.__file__),
        "script_sha256": sha256(__file__),
        "files": {"diagnostic.npz": sha256(args.output / "diagnostic.npz")},
    }
    write_json(args.output / "report.json", record)
    signal.alarm(0)
    if not record["passed"]:
        raise SystemExit("Bias scaling check failed; inspect preserved report")
    print(metrics)


if __name__ == "__main__":
    main()
