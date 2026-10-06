"""Export the actual CoCoA LSST Y1 survey and CAMB tables, without a covariance.

Run after start_cocoa.sh in the Cocoa environment. Power arrays retain
CoCoA's k in h/Mpc, P in (Mpc/h)^3 and [redshift, wavenumber] axes.
The original catalog midpoints are copied; n(z) is never approximated.
"""

import argparse
import sys
from pathlib import Path

from common import require_thread_environment, revision, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-bin", type=int, choices=range(1, 6), default=3)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error("Choose a fresh output directory; existing inputs stay intact")
    threads = require_thread_environment()

    import numpy as np

    cocoa = args.cocoa.resolve()
    project = cocoa / "projects/lsst_y1"
    core = cocoa / "external_modules/code/cosmolike_core"
    sys.path[:0] = [str(project), str(project / "covariance"), str(core)]
    import camb
    import cosmolike_lsst_y1_interface as interface
    from lsst_y1_covariance import configuration, initialize

    # The project's own setup fixes the cosmology and all accuracy choices.
    # Disable non-Limber and IA explicitly for this Gaussian assembly test.
    settings = configuration(gaussian={"nonlimber": False, "ia": "none"})
    if settings["cosmology"]["mnu"] != 0:
        raise ValueError("This first comparison requires massless neutrinos")
    if settings["photoz_zmid"] != 1:
        raise ValueError("Recheck the catalog convention: expected midpoint z")
    tables = initialize(interface=interface, settings=settings)

    source_file = project / settings["source_file"]
    lens_file = project / settings["lens_file"]
    source = np.loadtxt(source_file)
    lenses = np.loadtxt(lens_file)
    for filename, values in ((source_file, source), (lens_file, lenses)):
        if (values.ndim != 2 or values.shape[1] != 6
                or not np.all(np.isfinite(values)) or np.any(values < 0)
                or np.any(np.diff(values[:, 0]) <= 0)):
            raise ValueError(f"Malformed five-bin LSST catalog: {filename}")

    # The C input stores z fastest. Reshaping with order='F' restores the
    # physical [z,k] axes; neither interpolation nor unit conversion occurs.
    z = np.asarray(tables["z_2D"])
    k = 10.0**np.asarray(tables["log10k_2D"])
    linear = np.exp(np.asarray(tables["lnP_linear"]).reshape(
        (z.size, k.size), order="F"))
    nonlinear = np.exp(np.asarray(tables["lnP_nonlinear"]).reshape(
        (z.size, k.size), order="F"))
    for power in (linear, nonlinear):
        if not np.all(np.isfinite(power)) or np.any(power <= 0):
            raise ValueError("CAMB export contains invalid matter power")
    if np.any(np.diff(z) <= 0) or z[0] != 0 or np.any(np.diff(k) <= 0):
        raise ValueError("Need increasing k,z and a z=0 endpoint")

    output.mkdir(parents=True)
    np.savez_compressed(output / "inputs.npz", z=z, k_h_mpc=k,
                        p_linear=linear, p_nonlinear=nonlinear,
                        source_nz=source[:, [0, args.source_bin]],
                        lens_nz=lenses[:, [0, 1, 2]])
    # Archive the original interchange arrays for future background/growth
    # diagnostics. The first CCL call below imports only the power tables.
    np.savez_compressed(output / "camb_tables.npz", **tables)
    manifest = {
        "schema": "lsst-y1-tjpcov-inputs-v1",
        "survey": "LSST Y1 forecast subset",
        "source_bin_1based": args.source_bin, "lens_bins_1based": [1, 2],
        "area_deg2": settings["area_deg2"],
        "source_density_arcmin2": settings["source_density_arcmin2"][
            args.source_bin - 1],
        "sigma_e_component": settings["sigma_e_component"][args.source_bin - 1],
        "lens_density_arcmin2": settings["lens_density_arcmin2"][:2],
        "bias": settings["bias"][:2], "cosmology": settings["cosmology"],
        "gaussian": settings["gaussian"],
        "rsd": False, "magnification": False, "galaxy_bias": "linear",
        "units": {"k": "h/Mpc", "power": "(Mpc/h)^3",
                  "nz": "z_mid then unrenormalized n(z)",
                  "power_axes": ["redshift", "wavenumber"]},
        "cocoa": revision(cocoa), "core": revision(core),
        "lsst_y1": revision(project), "threads": threads,
        "camb_version": camb.__version__,
        "camb_source": revision(cocoa / "external_modules/code/CAMB"),
        "source_nz_sha256": sha256(source_file),
        "lens_nz_sha256": sha256(lens_file),
        "interface_sha256": sha256(interface.__file__),
        "script_sha256": sha256(__file__),
        "files": {name: sha256(output / name)
                  for name in ("inputs.npz", "camb_tables.npz")},
    }
    write_json(output / "manifest.json", manifest)
    print(f"Saved LSST inputs in {output}; no covariance was computed")


if __name__ == "__main__":
    main()
