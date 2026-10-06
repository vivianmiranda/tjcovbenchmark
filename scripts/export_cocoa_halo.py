"""Export current CoCoA halo ingredients on the saved native CCL grid.

Run in CoCoA's environment. The survey initializer must reproduce the
saved CAMB tables bitwise before any halo comparison is accepted.
"""

import argparse
import signal
import sys
from pathlib import Path

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("native_run", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--integration-accuracy", type=int, choices=range(5),
                        default=0)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    threads = require_thread_environment()
    if threads > 6 or args.timeout <= 0:
        parser.error("Use at most six threads and a positive timeout")
    if args.output.exists():
        parser.error("Choose a fresh output directory")
    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    native = load_bundle(args.native_run, "tjpcov-halo-ingredients-v1")
    if (native["status"] != "completed" or native["input_manifest_sha256"]
            != sha256(args.inputs / "manifest.json")):
        raise ValueError("Native halo export has different or incomplete inputs")
    args.output.mkdir(parents=True)
    write_json(args.output / "run_settings.json", {
        "threads": threads, "timeout_seconds": args.timeout,
        "integration_accuracy": args.integration_accuracy,
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
    })
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)

    import numpy as np

    cocoa = args.cocoa.resolve()
    project = cocoa / "projects/lsst_y1"
    core = cocoa / "external_modules/code/cosmolike_core"
    sys.path[:0] = [str(project), str(project / "covariance"), str(core)]
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize

    settings = configuration(
        gaussian={"nonlimber": False, "ia": "none"},
        integration_accuracy=args.integration_accuracy)
    if settings["cosmology"] != inputs["cosmology"]:
        raise ValueError("CoCoA cosmology changed since CAMB export")
    tables = initialize(interface=ci, settings=settings)
    original = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    np.testing.assert_array_equal(tables["z_2D"], original["z"])
    np.testing.assert_array_equal(10.0**tables["log10k_2D"], original["k_h_mpc"])
    for name, key in (("linear", "lnP_linear"), ("nonlinear", "lnP_nonlinear")):
        power = np.exp(tables[key].reshape(original[f"p_{name}"].shape,
                                          order="F"))
        np.testing.assert_array_equal(power, original[f"p_{name}"])

    data = np.load(args.native_run / "halo.npz", allow_pickle=False)
    redshift, mass, k = data["redshift"], data["mass"], data["k"]
    if (redshift.shape != (3,) or mass.shape != (41,) or k.shape != (33,)
            or not np.array_equal(redshift, [0.1, 0.5, 1.0])
            or np.any(np.diff(mass) <= 0) or np.any(np.diff(k) <= 0)
            or data["sigma"].shape != (3, 41)
            or data["concentration"].shape != (3, 41)):
        raise ValueError("Unexpected native halo grid or axes")

    # CoCoA's internal distances use c/H0 rather than Mpc/h. Convert
    # wavenumber and power using the constants in its initialized structs.
    length = 2997.92458
    rho = 7.4775e21*settings["cosmology"]["omegam"]/length**3
    rows = []
    for index, z in enumerate(redshift):
        a = float(1 / (1+z))
        sigma = np.sqrt([ci.sigma2(float(m), a, 1) for m in mass])
        peak = 1.686 / sigma
        slope = np.array([ci.dlognudlogm(float(m), a) for m in mass])
        multiplicity = np.array([ci.fnu(float(nu), a) for nu in peak])
        concentration = np.array([ci.conc(float(m), a) for m in mass])
        # The native profile uses CoCoA's own concentration. A second
        # diagnostic changes only concentration to CCL's Duffy values,
        # so an NFW formula difference is separated from a fit difference.
        profiles = {}
        for name, values in (("profile", concentration),
                             ("profile_matched_c", data["concentration"][index])):
            profiles[name] = np.array([
                [ci.u_nfw_c(float(c), float(w*length), float(m), a)
                 for m, c in zip(mass, values)] for w in k])

        # The two codes round the mean-density constant differently.
        # At fixed M200m, R scales as rho^(-1/3). Rescaling only this
        # diagnostic's k therefore matches kR as well as concentration;
        # it leaves every native spectrum and moment unchanged.
        radius_match = (rho/data["rho"][index])**(1/3)
        profiles["profile_matched_radius"] = np.array([
            [ci.u_nfw_c(float(c), float(w*length*radius_match), float(m), a)
             for m, c in zip(mass, data["concentration"][index])] for w in k])
        row = {
            "sigma": sigma, "dndlnm": rho/mass*multiplicity*peak*slope,
            "bias": np.array([ci.hb1nu(float(nu), a) for nu in peak]),
            "concentration": concentration, "rho": rho,
            "linear": ci.covariance.covariance_power(
                a=a, k=np.ascontiguousarray(k*length), linear=True)*length**3,
            "bias_matched_sigma": np.array([
                ci.hb1nu(float(1.686/s), a) for s in data["sigma"][index]]),
            "bias_matched_nu": np.array([
                ci.hb1nu(float(native["native_model"]["delta_c"]/s), a)
                for s in data["sigma"][index]]),
            **profiles,
        }
        rows.append(row)
        print(f"Exported CoCoA halo statistics at z={z}", flush=True)
    arrays = {key: np.asarray([row[key] for row in rows]) for key in rows[0]}

    # Request one k per batch row. This computes just the diagonal
    # I02/I12 needed here, without allocating unused k-by-k moments.
    a = np.repeat(1/(1+redshift), len(k))
    waves = np.tile(k, len(redshift))[:, None]*length
    single, moments = ci.covariance.covariance_halo_moments(
        a=a, k=np.ascontiguousarray(waves), lnm_edges=settings["lnm_edges"],
        nquad=settings["halo_mass_nquad"])
    shape = (len(redshift), len(k))
    arrays["i11"] = single.reshape(shape)
    arrays["i02"] = moments[0].reshape(shape)*length**3
    arrays["i12"] = moments[1].reshape(shape)*length**3
    if any(not np.all(np.isfinite(value)) for value in arrays.values()):
        raise ValueError("Non-finite CoCoA halo ingredient")
    np.savez_compressed(args.output / "halo.npz", redshift=redshift,
                        mass=mass, k=k, **arrays)
    record = {
        "schema": "cocoa-halo-ingredients-v1", "status": "completed",
        "scope": "Native halo ingredients and explicit matched-fit diagnostics",
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "native_manifest_sha256": sha256(args.native_run / "manifest.json"),
        "power_table_check": "linear and nonlinear CAMB tables bitwise equal",
        "native_model": {
            "mass_definition": "M200m", "mass_function": "Tinker10",
            "bias": "Tinker10", "concentration": "Bhattacharya13",
            "delta_c": 1.686,
            "completion": "Current production Wynn I11; direct higher moments",
        },
        "controls": {"integration_accuracy": args.integration_accuracy,
                     "halo_mass_nquad": settings["halo_mass_nquad"],
                     "lnm_edges": settings["lnm_edges"].tolist()},
        "units": native["units"], "axes": native["axes"],
        "core": revision(core), "lsst_y1": revision(project),
        "interface_sha256": sha256(ci.__file__),
        "script_sha256": sha256(__file__), "threads": threads,
        "timing_scope": "Accuracy only; no benchmark timing",
        "files": {"halo.npz": sha256(args.output / "halo.npz")},
    }
    write_json(args.output / "manifest.json", record)
    signal.alarm(0)


if __name__ == "__main__":
    main()
