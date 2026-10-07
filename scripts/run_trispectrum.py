"""Export the public CCL trispectra selected by TJPCov, without projection.

The native case keeps Tinker08, Duffy08 and separable growth. Optional
public-input changes are labeled diagnostics, never native TJPCov runs.
Outputs share CoCoA's nine k nodes, three redshifts and 45 unordered pairs.
"""

import argparse
import importlib.metadata
import signal
from pathlib import Path
from time import perf_counter

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("cocoa_run", type=Path)
    parser.add_argument("--tjpcov", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mass-refinement", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--direct-growth", action="store_true")
    parser.add_argument("--concentration", choices=("duffy08", "bhattacharya13"),
                        default="duffy08")
    parser.add_argument("--mass-function", choices=("tinker08", "tinker10"),
                        default="tinker08")
    parser.add_argument("--timing", action="store_true",
                        help="record first-use compute time on a quiet machine")
    args = parser.parse_args()
    threads = require_thread_environment()
    if threads > (8 if args.timing else 6) or args.output.exists():
        parser.error("Use at most six accuracy workers or eight timing workers, "
                     "and a new output directory")
    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    cocoa = load_bundle(args.cocoa_run, "cocoa-trispectrum-v1")
    if cocoa["input_manifest_sha256"] != sha256(args.inputs / "manifest.json"):
        raise ValueError("CoCoA used different inputs")
    args.output.mkdir(parents=True)
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(600)

    import numpy as np
    import pyccl as ccl
    import tjpcov

    source = Path(tjpcov.__file__).resolve().parent / "covariance_fourier_cNG.py"
    if sha256(source) != sha256(args.tjpcov / "tjpcov" / source.name):
        raise ValueError("Installed TJPCov differs from the checkout")
    setup_started = perf_counter()
    data = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    cosmos = inputs["cosmology"]
    h = cosmos["H0"]/100
    a_table = 1/(1+data["z"][::-1])
    if (cosmos["mnu"] != 0 or inputs["units"]["k"] != "h/Mpc"
            or inputs["units"]["power"] != "(Mpc/h)^3"
            or inputs["units"]["power_axes"] != ["redshift", "wavenumber"]):
        raise ValueError("Need massless neutrinos and the expected power units")
    powers = {
        name: {"a": a_table, "k": data["k_h_mpc"]*h,
               "delta_matter:delta_matter": data[f"p_{name}"][::-1]/h**3}
        for name in ("linear", "nonlinear")}
    cosmo = ccl.CosmologyCalculator(
        Omega_c=cosmos["omegam"]-cosmos["omegab"], Omega_b=cosmos["omegab"],
        h=h, n_s=cosmos["ns"], A_s=cosmos["As_1e9"]*1e-9, m_nu=0,
        w0=cosmos["w"], wa=cosmos["w0pwa"]-cosmos["w"],
        pk_linear=powers["linear"], pk_nonlin=powers["nonlinear"])
    md = ccl.halos.MassDef200m
    mf_class = (ccl.halos.MassFuncTinker08 if args.mass_function == "tinker08"
                else ccl.halos.MassFuncTinker10)
    concentration_class = (ccl.halos.ConcentrationDuffy08
                           if args.concentration == "duffy08"
                           else ccl.halos.ConcentrationBhattacharya13)
    hmc = ccl.halos.HMCalculator(
        mass_function=mf_class(mass_def=md), mass_def=md,
        halo_bias=ccl.halos.HaloBiasTinker10(mass_def=md),
        nM=args.mass_refinement*127+1)
    profile = ccl.halos.HaloProfileNFW(
        mass_def=md, concentration=concentration_class(mass_def=md),
        fourier_analytic=True)
    pair_profile = ccl.halos.Profile2pt()
    grid = np.load(args.cocoa_run / "trispectrum.npz", allow_pickle=False)
    k, redshift = grid["k"], grid["redshift"]
    first, second = grid["first"], grid["second"]
    wave = k*h
    a = 1/(1+redshift)
    functions = [ccl.halos.halomod_trispectrum_1h,
                 ccl.halos.halomod_trispectrum_2h_13,
                 ccl.halos.halomod_trispectrum_2h_22,
                 ccl.halos.halomod_trispectrum_3h,
                 ccl.halos.halomod_trispectrum_4h]
    setup_seconds = perf_counter() - setup_started
    construction_started = perf_counter()
    terms = []
    exchange = []
    for name, function in zip(cocoa["halo_order"], functions):
        options = {} if name in ("1h", "2h13") else {
            "separable_growth": not args.direct_growth}
        # CCL stores [a,k2,k1]. Save the common (K,Q) pair ordering.
        matrix = function(cosmo, hmc, wave, a, prof=profile, **options)*h**9
        terms.append(matrix[:, second, first])
        exchange.append(float(np.max(np.abs(matrix-matrix.swapaxes(1, 2)))
                              /np.max(np.abs(matrix))))
        print(f"CCL native function completed: {name}", flush=True)
    terms = np.moveaxis(np.array(terms), 0, 1)
    construction_seconds = perf_counter() - construction_started

    # Export CCL's actual moments and power evaluations. Supplying them
    # to CoCoA later isolates the angular integration/assembly from halo
    # abundance and profile choices; no CCL function is replaced.
    rows = []
    for index, aa in enumerate(a):
        norm = profile.get_normalization(cosmo, aa, hmc=hmc)
        single = hmc.I_1_1(cosmo, wave, aa, profile)/norm
        i12 = hmc.I_1_2(cosmo, wave, aa, profile,
                       prof_2pt=pair_profile, diag=False)/norm**2*h**3
        i13 = hmc.I_1_3(cosmo, wave, aa, profile,
                       prof_2pt=pair_profile)/norm**3*h**6
        moments = np.array([
            np.zeros(len(first)),  # I02 is unused by trispectrum assembly.
            i12[second, first], i13[second, first], i13[first, second],
            terms[index, 0],
        ])
        rows.append(dict(
            moments=moments, i11=single[[first, second]],
            power=ccl.linear_matter_power(cosmo, wave, aa)*h**3,
            internal_power=ccl.linear_matter_power(
                cosmo, grid["internal_k"].ravel()*h, aa).reshape(
                    grid["internal_k"].shape)*h**3))
    arrays = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    arrays["terms"] = terms
    if not all(np.isfinite(value).all() for value in arrays.values()):
        raise ValueError("Non-finite native trispectrum or ingredient")
    np.savez_compressed(args.output / "trispectrum.npz", k=k,
                        redshift=redshift, first=first, second=second, **arrays)
    pk_source = Path(ccl.__file__).resolve().parent / "halos/pk_4pt.py"
    write_json(args.output / "manifest.json", {
        "schema": "tjpcov-trispectrum-v1", "status": "completed",
        "scope": "Public CCL matter functions used by TJPCov cNG",
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "cocoa_manifest_sha256": sha256(args.cocoa_run / "manifest.json"),
        "halo_order": cocoa["halo_order"], "units": cocoa["units"],
        "axes": cocoa["axes"], "mass_precision": dict(hmc.precision),
        "mass_function": args.mass_function, "concentration": args.concentration,
        "separable_growth": not args.direct_growth,
        "tinker10_normalization": "CCL norm_all_z=False; not CoCoA normalization",
        "maximum_exchange_asymmetry_over_term_peak": exchange,
        "tjpcov": revision(args.tjpcov), "tjpcov_source_sha256": sha256(source),
        "ccl_source_sha256": sha256(pk_source),
        "ccl_binary_sha256": sha256(ccl._ccllib.__file__),
        "versions": {n: importlib.metadata.version(n) for n in ("pyccl", "tjpcov")},
        "threads": threads, "script_sha256": sha256(__file__),
        "timing_scope": "First five-term matter trispectrum" if args.timing
                        else "accuracy only",
        "timing": ({"setup_seconds": setup_seconds,
                    "construction_seconds": construction_seconds,
                    "scope": "All five native CCL halo terms at nine k nodes "
                             "and three redshifts; includes first halo tables; "
                             "excludes extra ingredient exports and projection"}
                   if args.timing else None),
        "files": {"trispectrum.npz": sha256(args.output / "trispectrum.npz")},
    })


if __name__ == "__main__":
    main()
