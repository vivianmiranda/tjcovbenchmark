"""Compute a five-band native TJPCov shear covariance and export its spectra.

Run in the separate TJPCov environment. This is the first assembly test,
not a complete LSST covariance or a native-spectrum CoCoA comparison.
No CCL/TJPCov functions are replaced. The SACC means are placeholders;
TJPCov computes the actual spectra with its usual CCL angular_cl calls.
"""

import argparse
import importlib.metadata
import inspect
import platform
from pathlib import Path
from time import perf_counter

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell-range", type=int, nargs=2, default=[30, 150],
                        metavar=("FIRST", "STOP"), help="exclusive upper edge")
    args = parser.parse_args()
    first, stop = args.ell_range
    if first < 2 or stop - first < 20 or stop - first > 200:
        parser.error("This bounded pilot needs 2 <= FIRST and 20..200 ell steps")
    output = args.output.resolve()
    if output.exists():
        parser.error("Choose a fresh --output; cached calculations stay intact")
    threads = require_thread_environment()
    inputs = args.inputs.resolve()
    manifest = load_bundle(inputs, "lsst-y1-tjpcov-inputs-v1")
    if (manifest["cosmology"]["mnu"] != 0 or manifest["gaussian"]["nonlimber"]
            or manifest["gaussian"]["ia"] != "none"
            or manifest["rsd"] or manifest["magnification"]):
        raise ValueError("Need massless neutrinos, Limber and zero IA/RSD/mag")
    if manifest["units"]["power_axes"] != ["redshift", "wavenumber"]:
        raise ValueError("Recheck the supplied matter-power axes")
    if (manifest["units"]["k"] != "h/Mpc"
            or manifest["units"]["power"] != "(Mpc/h)^3"):
        raise ValueError("Recheck the supplied matter-power units")

    import numpy as np
    import pyccl as ccl
    import sacc
    import tjpcov
    from tjpcov.covariance_gaussian_fsky import FourierGaussianFsky
    from tjpcov.wigner_transform import bin_cov

    # Match the installed routines to the checkout named in the report.
    # A clean Git hash alone would not identify a stale installed package.
    module_path = Path(tjpcov.__file__).resolve().parent
    source_hashes = {}
    for name in ("covariance_gaussian_fsky.py", "covariance_builder.py",
                 "wigner_transform.py", "covariance_io.py"):
        installed = sha256(module_path / name)
        if installed != sha256(args.tjpcov / "tjpcov" / name):
            raise ValueError(f"Installed TJPCov differs from checkout: {name}")
        source_hashes[name] = installed
    limber_default = inspect.signature(ccl.angular_cl).parameters[
        "l_limber"].default
    if limber_default != -1:
        raise ValueError("CCL angular_cl default changed; audit before running")

    started = perf_counter()
    data = np.load(inputs / "inputs.npz", allow_pickle=False)
    cosmology = manifest["cosmology"]
    h = cosmology["H0"] / 100
    z = data["z"]
    k = data["k_h_mpc"]
    if (z.ndim != 1 or k.ndim != 1 or z[0] != 0 or k[0] <= 0
            or np.any(np.diff(z) <= 0) or np.any(np.diff(k) <= 0)):
        raise ValueError("Malformed z,k axes in the CAMB export")

    # CCL uses Mpc rather than Mpc/h and increasing scale factor.
    # Reverse the redshift axis, multiply k by h and divide P by h^3.
    # Its background is native CCL; only P(k,z) is supplied from CAMB.
    powers = {}
    for name in ("linear", "nonlinear"):
        values = data[f"p_{name}"]
        if (values.shape != (z.size, k.size)
                or not np.all(np.isfinite(values)) or np.any(values <= 0)):
            raise ValueError(f"Malformed {name} matter-power table")
        powers[name] = {"a": 1 / (1 + z[::-1]), "k": k * h,
                        "delta_matter:delta_matter": values[::-1] / h**3}
    cosmo = ccl.CosmologyCalculator(
        Omega_c=cosmology["omegam"] - cosmology["omegab"],
        Omega_b=cosmology["omegab"], h=h, n_s=cosmology["ns"],
        A_s=cosmology["As_1e9"] * 1e-9, m_nu=0,
        w0=cosmology["w"], wa=cosmology["w0pwa"] - cosmology["w"],
        pk_linear=powers["linear"], pk_nonlin=powers["nonlinear"],
    )

    # TJPCov finds each lower edge from the first positive window value
    # and the last upper edge from the final positive value. Include the
    # endpoint STOP here; bin_cov then excludes it from the final band.
    # Keep this otherwise-unused node in the export, with operator weight 0.
    ell = np.arange(first, stop + 1, dtype=float)
    edges = np.geomspace(first, stop, 6).astype(int)
    if np.any(np.diff(edges) <= 0):
        raise ValueError("Requested bands are not distinct")
    windows = np.zeros((ell.size, 5))
    for band, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        windows[(ell >= lo) & (ell <= hi), band] = 1
    # Weighted centers label the SACC data; native covariance binning is
    # checked separately below because it reconstructs its own operator.
    centers = (ell @ windows) / windows.sum(axis=0)
    nz = data["source_nz"]
    if (nz.ndim != 2 or nz.shape[1] != 2
            or not np.all(np.isfinite(nz)) or np.any(nz < 0)
            or np.any(np.diff(nz[:, 0]) <= 0) or not np.any(nz[:, 1] > 0)):
        raise ValueError("Malformed source n(z)")
    source = f"source{manifest['source_bin_1based']}"
    catalog = sacc.Sacc()
    catalog.add_tracer("NZ", source, z=nz[:, 0], nz=nz[:, 1],
                       quantity="galaxy_shear", spin=2)
    catalog.add_ell_cl("cl_ee", source, source, centers, np.zeros(5),
                       window=sacc.BandpowerWindow(ell, windows))
    area_sr = manifest["area_deg2"] * (np.pi / 180)**2
    config = {"cosmo": cosmo, "sacc_file": catalog, "IA": None,
              "use_mpi": False, "fsky": area_sr / (4 * np.pi),
              f"Ngal_{source}": manifest["source_density_arcmin2"]}
    output.mkdir(parents=True)
    catalog.save_fits(str(output / "source_bins.fits"), overwrite=False)

    # A Gaussian covariance is quadratic in the noise power N. Three
    # public calls with N multiplied by 0,1,2 separate CC, CN and NN.
    # Vary sigma_e by sqrt(factor), since N=sigma_e^2/number_density.
    # This diagnostic adds calls; it is not a one-call production timing.
    calls = []
    call_seconds = []
    for factor in (0, 1, 2):
        options = dict(config)
        options[f"sigma_e_{source}"] = (
            manifest["sigma_e_component"] * np.sqrt(factor))
        options["outdir"] = str(output / f"native_noise_{factor}")
        native = FourierGaussianFsky({"tjpcov": options})
        actual_ell, effective_ell, actual_edges = native.get_binning_info()
        if (not np.array_equal(actual_ell, ell)
                or not np.array_equal(actual_edges, edges)):
            raise ValueError("TJPCov reconstructed different ell nodes or edges")
        before = perf_counter()
        calls.append(native.get_covariance_block(
            (source, source), (source, source), include_b_modes=False))
        call_seconds.append(perf_counter() - before)
        if factor == 1:
            tracer, noise = native.get_tracer_info()
            expected_noise = (manifest["sigma_e_component"]**2
                              * (np.pi / 10800)**2
                              / manifest["source_density_arcmin2"])
            if not np.isclose(noise[source], expected_noise, rtol=1e-14, atol=0):
                raise ValueError("TJPCov shape-noise units do not match")
            spectra = ccl.angular_cl(cosmo, tracer[source], tracer[source], ell)

    # With integer ell nodes and edges, native bin_cov weights each node
    # by ell*dell = ell. Its supplied SACC window amplitudes are not used.
    operators = np.zeros((5, ell.size))
    for band, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        inside = (ell >= lo) & (ell < hi)
        operators[band, inside] = ell[inside] / ell[inside].sum()
    # Check that transcription against TJPCov's own operator routine.
    # This is a binning contract test, not an independent covariance code.
    probe = np.outer(ell, 1 / ell)
    _, native_binned = bin_cov(r=ell, cov=probe, r_bins=actual_edges)
    np.testing.assert_allclose(operators @ probe @ operators.T, native_binned,
                               rtol=5e-14, atol=0)
    components = {"sample_variance": calls[0], "total": calls[1],
                  "noise": (calls[2] - 2 * calls[1] + calls[0]) / 2}
    components["mixed"] = calls[1] - calls[0] - components["noise"]
    for matrix in components.values():
        if matrix.shape != (5, 5) or not np.all(np.isfinite(matrix)):
            raise ValueError("Invalid TJPCov Gaussian component")
    np.savez_compressed(output / "gaussian.npz", ell=ell, edges=edges,
                        effective_ell=effective_ell, operators=operators,
                        spectra=spectra[:, None, None],
                        noise_power=np.array([expected_noise]),
                        pairs=np.array([[0, 0]], dtype=np.int32), **components)
    record = {
        "schema": "tjpcov-gaussian-shear-v1", "status": "completed",
        "scope": "Gaussian EE assembly only; supplied CAMB P; native CCL C_ell",
        "input_manifest": manifest,
        "input_manifest_sha256": sha256(inputs / "manifest.json"),
        "tjpcov": revision(args.tjpcov), "installed_source": source_hashes,
        "versions": {name: importlib.metadata.version(name)
                     for name in ("tjpcov", "pyccl", "sacc", "numpy", "scipy")},
        "python": platform.python_version(), "platform": platform.platform(),
        "ccl_binary_sha256": sha256(ccl._ccllib.__file__),
        "ccl_angular_cl_sha256": sha256(inspect.getfile(ccl.angular_cl)),
        "ccl_l_limber_default": limber_default,
        "background": "native CCL; not imported from CAMB",
        "area_sr": area_sr, "source_name": source,
        "ordering": "one source EE spectrum, then five increasing ell bands",
        "estimator": "native TJPCov ell-weighted average; integer ell nodes",
        "noise_extraction": "N factors 0,1,2; NN=(C2-2*C1+C0)/2; CN=C1-C0-NN",
        "threads": threads,
        "diagnostic_call_seconds": call_seconds,
        "diagnostic_elapsed_seconds": perf_counter() - started,
        "timing_scope": "diagnostic only; includes three noise configurations",
        "script_sha256": sha256(__file__),
        "files": {name: sha256(output / name)
                  for name in ("gaussian.npz", "source_bins.fits")},
    }
    write_json(output / "manifest.json", record)
    print(f"Saved native TJPCov Gaussian components and exact inputs to {output}")


if __name__ == "__main__":
    main()
