"""Bounded native TJPCov shear covariance, with explicit cNG table controls.

Gaussian bands use TJPCov's ell-weighted discrete average. SSC and cNG
use its native effective-ell prescription. This is not a uniformly
band-averaged covariance. All calculators are unmodified public classes.
"""
import argparse
import inspect
import platform
import resource
import signal
import sys
from pathlib import Path

from common import load_bundle, require_thread_environment, revision, sha256, write_json


def cosmology(inputs, a_refinement=1, k_density=None):
    """Install shared CAMB power; retain native CCL background and readers."""
    import numpy as np
    import pyccl as ccl
    m = load_bundle(inputs, "lsst-y1-tjpcov-inputs-v1")
    d = np.load(inputs / "inputs.npz", allow_pickle=False)
    p = m["cosmology"]
    h = p["H0"] / 100
    a = 1 / (1 + d["z"][::-1])
    ccl.spline_params.A_SPLINE_MINLOG_PK = float(a[0])
    ccl.spline_params.A_SPLINE_NLOG_PK = 10 * a_refinement + 1
    ccl.spline_params.A_SPLINE_NA_PK = 39 * a_refinement + 1
    if k_density is not None:
        ccl.spline_params.N_K = k_density
    pk = {name: {"a": a, "k": d["k_h_mpc"] * h,
                 "delta_matter:delta_matter": d[f"p_{name}"][::-1] / h**3}
          for name in ("linear", "nonlinear")}
    c = ccl.CosmologyCalculator(
        Omega_c=p["omegam"]-p["omegab"], Omega_b=p["omegab"], h=h,
        n_s=p["ns"], A_s=p["As_1e9"]*1e-9, m_nu=0,
        w0=p["w"], wa=p["w0pwa"]-p["w"],
        pk_linear=pk["linear"], pk_nonlin=pk["nonlinear"])
    return c, m, d


def catalog_and_config(cosmo, m, d, output):
    import numpy as np
    import sacc
    ell = np.arange(30., 3001., 30.)
    grid = np.arange(15., 3016.)
    edges = np.arange(15., 3016., 30.)
    windows = np.array([(grid >= lo) & (grid <= hi)
                        for lo, hi in zip(edges[:-1], edges[1:])], float).T
    source = f"source{m['source_bin_1based']}"
    cat = sacc.Sacc()
    cat.add_tracer("NZ", source, z=d["source_nz"][:, 0],
                   nz=d["source_nz"][:, 1], quantity="galaxy_shear", spin=2)
    cat.add_ell_cl("cl_ee", source, source, ell, np.zeros(len(ell)),
                   window=sacc.BandpowerWindow(grid, windows))
    cfg = {"tjpcov": {"cosmo": cosmo, "sacc_file": cat, "IA": None,
           "use_mpi": False, "outdir": str(output),
           "fsky": m["area_deg2"]*(np.pi/180)**2/(4*np.pi),
           f"Ngal_{source}": m["source_density_arcmin2"],
           f"sigma_e_{source}": m["sigma_e_component"]}}
    return cat, source, cfg, ell, grid, edges


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("--tjpcov", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--components", nargs="+", choices=("gaussian", "ssc", "cng"),
                   default=["gaussian", "ssc", "cng"])
    p.add_argument("--k-density", type=int, default=16,
                   help="Explicit CCL N_K for cNG only; full default k range retained")
    p.add_argument("--a-refinement", type=int, choices=(1, 2, 4, 8, 16), default=2)
    p.add_argument("--integration-method", choices=("qag_quad", "spline"), default="qag_quad")
    p.add_argument("--timeout", type=int, default=600)
    args = p.parse_args()
    threads = require_thread_environment()
    if not 1 <= threads <= 8 or args.output.exists() or args.k_density < 4:
        p.error("Use at most eight workers, positive k density >=4 and fresh output")
    args.output.mkdir(parents=True)
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)
    import numpy as np
    import pyccl as ccl
    import tjpcov
    import yaml
    from tjpcov.covariance_gaussian_fsky import FourierGaussianFsky
    from tjpcov.covariance_fourier_ssc_fsky import FourierSSCHaloModelFsky
    from tjpcov.covariance_fourier_cNG_fsky import FouriercNGHaloModelFsky
    installed = Path(tjpcov.__file__).resolve().parent
    sources = {}
    for name in ("covariance_fourier_cNG.py", "covariance_fourier_cNG_fsky.py",
                 "covariance_fourier_ssc.py", "covariance_fourier_ssc_fsky.py",
                 "covariance_gaussian_fsky.py", "covariance_builder.py"):
        sources[name] = sha256(installed/name)
        if sources[name] != sha256(args.tjpcov/"tjpcov"/name):
            raise ValueError(f"Installed source differs: {name}")
    default_density = int(ccl.spline_params.N_K)
    arrays, records = {}, {}
    for component in args.components:
        # SSC already has an independent 393->785-node convergence study.
        refinement = 16 if component == "ssc" else args.a_refinement
        density = args.k_density if component == "cng" else default_density
        cosmo, m, d = cosmology(args.inputs, refinement, density)
        cat, source, cfg, ell, grid, edges = catalog_and_config(
            cosmo, m, d, args.output/component)
        if m["cosmology"]["mnu"] or m["gaussian"]["ia"] != "none":
            raise ValueError("This pilot requires massless neutrinos and no IA")
        arrays.update(ell=ell, edges=edges)
        cfg["SSC"] = cfg["cNG"] = {"integration_method": args.integration_method}
        cfg["HOD"] = yaml.safe_load((args.tjpcov/
            "examples/full_3x2pt_cov_example.yml").read_text())["HOD"]
        cls = {"gaussian": FourierGaussianFsky, "ssc": FourierSSCHaloModelFsky,
               "cng": FouriercNGHaloModelFsky}[component]
        calc = cls(cfg)
        a, lk = cosmo.get_pk_spline_a(), cosmo.get_pk_spline_lk()
        print(f"{component}: {len(a)} full time nodes, {len(lk)} k nodes; "
              f"k=[{np.exp(lk[0]):.6g},{np.exp(lk[-1]):.6g}] Mpc^-1", flush=True)
        record = {"a_refinement": refinement, "N_K": density,
                  "default_N_K": default_density, "a_count": len(a),
                  "k_count": len(lk), "k_range_Mpc_inverse": np.exp(lk[[0,-1]]).tolist(),
                  "integration_method": args.integration_method}
        write_json(args.output/"progress.json", {"active": component, **record})
        captured = {}
        code = inspect.unwrap(calc.get_covariance_block).__code__
        def capture(frame, event, arg):
            # Read locals only on return; no library function or data is replaced.
            if frame.f_code is code and event == "return" and component == "cng":
                for key in ("a_arr", "lk_arr", "tkk"):
                    if key in frame.f_locals:
                        captured[key] = frame.f_locals[key]
        if component == "cng":
            sys.setprofile(capture)
        try:
            matrix = calc.get_covariance_block((source,source), (source,source),
                                                include_b_modes=False)
        finally:
            sys.setprofile(None)
        if matrix.shape != (100,100) or not np.isfinite(matrix).all():
            raise ValueError("Invalid native matrix")
        arrays[component] = matrix
        record["max_asymmetry_over_peak"] = float(np.max(abs(matrix-matrix.T))/np.max(abs(matrix)))
        if component == "cng":
            if set(captured) != {"a_arr", "lk_arr", "tkk"}:
                raise ValueError("Read-only table observation failed")
            np.savez_compressed(args.output/"native_trispectrum.npz", **captured)
            record["retained_a_count"] = len(captured["a_arr"])
            record["observer"] = "Read-only return-frame observation; not a timing run"
        if component == "gaussian":
            actual, centers, actual_edges = calc.get_binning_info()
            np.testing.assert_array_equal(actual, grid)
            np.testing.assert_array_equal(centers, ell)
            np.testing.assert_array_equal(actual_edges, edges)
            tracers, noise = calc.get_tracer_info()
            arrays["spectra_ell"] = grid
            arrays["spectra"] = ccl.angular_cl(cosmo, tracers[source], tracers[source], grid)
            arrays["noise"] = np.asarray(noise[source])
            arrays["signal"] = ccl.angular_cl(cosmo, tracers[source], tracers[source], ell)
        records[component] = record
        np.savez_compressed(args.output/"covariance.npz", **arrays)
        print(f"{component}: saved all {matrix.size} entries", flush=True)
    if set(args.components) == {"gaussian","ssc","cng"}:
        arrays["total"] = sum(arrays[x] for x in args.components)
        np.savez_compressed(args.output/"covariance.npz", **arrays)
    write_json(args.output/"manifest.json", {
        "schema":"tjpcov-complete-shear-v1", "status":"completed",
        "scope":"Native Gaussian band averages plus center-sampled SSC/cNG",
        "components":records, "shape":[100,100], "source_name":source,
        "input_manifest_sha256":sha256(args.inputs/"manifest.json"),
        "input_manifest":m, "tjpcov":revision(args.tjpcov),
        "installed_source":sources, "ccl_version":ccl.__version__,
        "ccl_binary_sha256":sha256(ccl._ccllib.__file__),
        "ccl_pk4pt_sha256":sha256(Path(ccl.__file__).parent/"halos/pk_4pt.py"),
        "script_sha256":sha256(__file__), "threads":threads,
        "platform":platform.platform(), "timing":"accuracy only",
        "peak_rss_bytes":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "files":{f.name:sha256(f) for f in args.output.glob("*.npz")}})
    signal.alarm(0)


if __name__ == "__main__":
    main()
