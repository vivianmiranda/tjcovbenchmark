"""Export public CCL halo ingredients using the choices fixed by TJPCov.

These are native model ingredients, not a supplied-moment equality test.
Mass is Msun/h, k is h/Mpc, abundance is (h/Mpc)^3 per lnM, profiles
are dimensionless, and power-like moments are (Mpc/h)^3.
"""

import argparse
import importlib.metadata
import platform
import signal
from pathlib import Path

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mass-refinement", type=int, choices=(1, 2, 4),
                        default=1, help="refine native mass intervals")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    threads = require_thread_environment()
    if threads > 6 or args.timeout <= 0:
        parser.error("Use at most six threads and a positive timeout")
    if args.output.exists():
        parser.error("Choose a fresh output directory")
    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    if inputs["cosmology"]["mnu"] != 0:
        raise ValueError("This pilot requires massless neutrinos")
    args.output.mkdir(parents=True)
    write_json(args.output / "run_settings.json", {
        "threads": threads, "timeout_seconds": args.timeout,
        "mass_refinement": args.mass_refinement,
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
    })
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)

    import numpy as np
    import pyccl as ccl
    import tjpcov

    # Check the installed TJPCov source that selects these CCL models.
    # We call the public CCL ingredients, not a copied halo formula.
    package = Path(tjpcov.__file__).resolve().parent
    source = "covariance_fourier_ssc.py"
    if sha256(package / source) != sha256(args.tjpcov / "tjpcov" / source):
        raise ValueError("Installed TJPCov differs from the named checkout")
    data = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    cosmos = inputs["cosmology"]
    h = cosmos["H0"] / 100
    a_table = 1 / (1 + data["z"][::-1])
    wave = data["k_h_mpc"] * h
    if (inputs["units"]["k"] != "h/Mpc"
            or inputs["units"]["power"] != "(Mpc/h)^3"
            or inputs["units"]["power_axes"] != ["redshift", "wavenumber"]
            or np.any(np.diff(a_table) <= 0) or a_table[-1] != 1
            or wave[0] <= 0 or np.any(np.diff(wave) <= 0)):
        raise ValueError("Unexpected power units or grid")
    powers = {}
    for name in ("linear", "nonlinear"):
        values = data[f"p_{name}"]
        if (values.shape != (a_table.size, wave.size)
                or not np.all(np.isfinite(values)) or np.any(values <= 0)):
            raise ValueError(f"Invalid {name} CAMB power")
        powers[name] = {"a": a_table, "k": wave,
                        "delta_matter:delta_matter": values[::-1] / h**3}
    cosmo = ccl.CosmologyCalculator(
        Omega_c=cosmos["omegam"]-cosmos["omegab"],
        Omega_b=cosmos["omegab"], h=h, n_s=cosmos["ns"],
        A_s=cosmos["As_1e9"]*1e-9, m_nu=0, w0=cosmos["w"],
        wa=cosmos["w0pwa"]-cosmos["w"],
        pk_linear=powers["linear"], pk_nonlin=powers["nonlinear"],
    )

    # TJPCov's native SSC fixes these fits, M200m and analytic NFW.
    # Its default HMCalculator has 128 nodes from 1e8 to 1e16 Msun.
    # Refinement adds intervals without moving the original endpoints.
    mass_def = ccl.halos.MassDef200m
    mass_function = ccl.halos.MassFuncTinker08(mass_def=mass_def)
    halo_bias = ccl.halos.HaloBiasTinker10(mass_def=mass_def)
    concentration = ccl.halos.ConcentrationDuffy08(mass_def=mass_def)
    profile = ccl.halos.HaloProfileNFW(
        mass_def=mass_def, concentration=concentration, fourier_analytic=True)
    calculator = ccl.halos.HMCalculator(
        mass_function=mass_function, halo_bias=halo_bias, mass_def=mass_def,
        nM=args.mass_refinement*(128-1)+1)
    pair_profile = ccl.halos.Profile2pt()
    redshift = np.array([0.1, 0.5, 1.0])
    mass = np.geomspace(1e10, 1e15, 41)
    k = np.geomspace(0.001, 10, 33)
    if a_table[0] > 1 / (1 + redshift[-1]):
        raise ValueError("CAMB does not cover the requested redshifts")
    rows = []
    for z in redshift:
        a = float(1 / (1+z))
        mass_msun = mass / h
        k_mpc = k * h
        norm = profile.get_normalization(cosmo, a, hmc=calculator)
        # CCL's Fourier NFW carries mass; divide by M to obtain u(k|M).
        # Its halo integrals carry one factor of mean density per profile.
        # Remove those factors before converting Mpc^3 to (Mpc/h)^3.
        row = {
            "sigma": ccl.sigmaM(cosmo, mass_msun, a),
            "dndlnm": mass_function(cosmo, mass_msun, a)/np.log(10)/h**3,
            "bias": halo_bias(cosmo, mass_msun, a),
            "concentration": concentration(cosmo, mass_msun, a),
            "profile": (profile.fourier(cosmo, k_mpc, mass_msun, a)
                        / mass_msun[:, None]).T,
            "linear": ccl.linear_matter_power(cosmo, k_mpc, a)*h**3,
            "rho": norm/h**2,
            "i01": calculator.I_0_1(cosmo, k_mpc, a, profile)/norm,
            "i11": calculator.I_1_1(cosmo, k_mpc, a, profile)/norm,
            "i02": calculator.I_0_2(
                cosmo, k_mpc, a, profile, prof_2pt=pair_profile)/norm**2*h**3,
            "i12": calculator.I_1_2(
                cosmo, k_mpc, a, profile, prof_2pt=pair_profile)/norm**2*h**3,
        }
        if any(not np.all(np.isfinite(value)) for value in row.values()):
            raise ValueError(f"Non-finite halo ingredient at z={z}")
        rows.append(row)
        print(f"Exported native CCL halo ingredients at z={z}", flush=True)
    arrays = {key: np.asarray([row[key] for row in rows]) for key in rows[0]}
    np.savez_compressed(args.output / "halo.npz", redshift=redshift,
                        mass=mass, k=k, **arrays)
    ccl_package = Path(ccl.__file__).resolve().parent
    sources = ("power.py", "halos/halo_model.py", "halos/halo_model_base.py",
               "halos/hmfunc/tinker08.py", "halos/hbias/tinker10.py",
               "halos/concentration/duffy08.py", "halos/profiles/nfw.py",
               "halos/profiles/profile_base.py", "halos/profiles_2pt.py")
    record = {
        "schema": "tjpcov-halo-ingredients-v1", "status": "completed",
        "scope": "Public CCL ingredients selected by native TJPCov SSC",
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "native_model": {
            "mass_definition": "M200m", "mass_function": "Tinker08",
            "bias": "Tinker10", "concentration": "Duffy08",
            "profile": "analytic truncated NFW",
            "delta_c": float(ccl.halos.get_delta_c(None, None, kind="EdS")),
            "completion": "CCL additive minimum-mass terms in all moments",
        },
        "controls": dict(calculator.precision),
        "mass_refinement": args.mass_refinement,
        "units": {"mass": "Msun/h", "k": "h/Mpc",
                  "dndlnm": "(h/Mpc)^3", "rho": "(Msun/h)/(Mpc/h)^3",
                  "linear_i02_i12": "(Mpc/h)^3", "profile_i01_i11": "1"},
        "axes": {"mass_statistics": ["redshift", "mass"],
                 "power_moments": ["redshift", "k"],
                 "profile": ["redshift", "k", "mass"]},
        "tjpcov": revision(args.tjpcov),
        "tjpcov_source_sha256": sha256(package / source),
        "ccl_source_sha256": {name: sha256(ccl_package / name)
                              for name in sources},
        "ccl_binary_sha256": sha256(ccl._ccllib.__file__),
        "versions": {name: importlib.metadata.version(name)
                     for name in ("tjpcov", "pyccl", "numpy", "scipy")},
        "platform": platform.platform(), "threads": threads,
        "script_sha256": sha256(__file__),
        "timing_scope": "Accuracy only; no benchmark timing",
        "files": {"halo.npz": sha256(args.output / "halo.npz")},
    }
    write_json(args.output / "manifest.json", record)
    signal.alarm(0)


if __name__ == "__main__":
    main()
