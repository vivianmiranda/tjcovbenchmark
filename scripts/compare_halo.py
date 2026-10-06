"""Compare native TJPCov/CCL and CoCoA halo ingredients on identical grids.

The fits and finite-mass completions remain native to each code. Ratios
are diagnostics, not an equality test. The matched-concentration profile
and matched-peak-height bias checks isolate two common analytic formulas.
"""

import argparse
from pathlib import Path

from common import load_bundle, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tjpcov", type=Path)
    parser.add_argument("cocoa", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a fresh comparison directory")

    native_record = load_bundle(args.tjpcov, "tjpcov-halo-ingredients-v1")
    cocoa_record = load_bundle(args.cocoa, "cocoa-halo-ingredients-v1")
    if (native_record["input_manifest_sha256"]
            != cocoa_record["input_manifest_sha256"]):
        raise ValueError("Halo exports used different LSST/CAMB bundles")
    # Matched-profile and matched-peak tests use inputs copied from this
    # exact native export, not merely a cosmology with the same CAMB tables.
    if (cocoa_record["native_manifest_sha256"]
            != sha256(args.tjpcov / "manifest.json")):
        raise ValueError("CoCoA diagnostics used a different native halo export")

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    native = np.load(args.tjpcov / "halo.npz", allow_pickle=False)
    cocoa = np.load(args.cocoa / "halo.npz", allow_pickle=False)
    for key in ("redshift", "mass", "k"):
        if not np.array_equal(native[key], cocoa[key]):
            raise ValueError(f"Different {key} grid")
    for data in (native, cocoa):
        if any(not np.isfinite(data[key]).all() for key in data.files):
            raise ValueError("Nonfinite halo ingredient")

    quantities = (
        ("sigma", r"Mass rms fluctuation $\sigma(M)$", "mass"),
        ("dndlnm", r"Abundance $dn/d\ln M$", "mass"),
        ("bias", "Native halo bias", "mass"),
        ("concentration", "Native concentration", "mass"),
        ("linear", r"Linear power $P_{\rm lin}$", "k"),
        ("i11", r"Moment $I_1^1$", "k"),
        ("i02", r"Moment $I_2^0(k,k)$", "k"),
        ("i12", r"Moment $I_2^1(k,k)$", "k"),
    )
    rows = []
    for row, redshift in enumerate(native["redshift"]):
        values = {"redshift": float(redshift)}
        for key, _, grid in quantities:
            expected = (native["redshift"].size, native[grid].size)
            if native[key].shape != expected or cocoa[key].shape != expected:
                raise ValueError(f"Unexpected {key} axes")
            if np.any(native[key][row] <= 0):
                raise ValueError(f"Cannot form a fractional ratio for {key}")
            ratio = cocoa[key][row] / native[key][row] - 1
            values[key] = {
                "fractional_range": [float(ratio.min()), float(ratio.max())],
                "max_absolute_fractional": float(np.max(np.abs(ratio))),
            }

        # NFW profiles can approach zero, so use an absolute error in the
        # dimensionless Fourier profile u, not an unstable fractional one.
        profile_delta = cocoa["profile_matched_c"][row] - native["profile"][row]
        values["matched_concentration_profile_max_absolute"] = float(
            np.max(np.abs(profile_delta)))
        values["matched_radius_profile_max_absolute"] = float(np.max(np.abs(
            cocoa["profile_matched_radius"][row] - native["profile"][row])))
        values["mean_density_fractional_difference"] = float(
            cocoa["rho"][row]/native["rho"][row]-1)
        values["matched_sigma_bias_max_absolute_fractional"] = float(np.max(
            np.abs(cocoa["bias_matched_sigma"][row]/native["bias"][row]-1)))
        bias_delta = cocoa["bias_matched_nu"][row] / native["bias"][row] - 1
        values["matched_peak_bias_max_absolute_fractional"] = float(
            np.max(np.abs(bias_delta)))
        rows.append(values)

    args.output.mkdir(parents=True)
    archive = {"tjpcov__"+key: native[key] for key in native.files}
    archive.update({"cocoa__"+key: cocoa[key] for key in cocoa.files})
    np.savez_compressed(args.output / "halo.npz", **archive)
    write_json(args.output / "comparison.json", {
        "schema": "cocoa-tjpcov-halo-comparison-v1",
        "scope": "native halo ingredients; different fits and completions",
        "difference": "CoCoA / TJPCov-CCL - 1",
        "units": native_record["units"],
        "tjpcov": native_record, "cocoa": cocoa_record,
        "rows": rows, "script_sha256": sha256(__file__),
        "files": {"halo.npz": sha256(args.output / "halo.npz")},
        "timing_policy": "Correctness campaign; no benchmark times",
    })

    fig, axes = plt.subplots(3, 3, figsize=(12, 10), layout="constrained")
    for axis, (key, title, grid) in zip(axes.ravel(), quantities):
        for row, z in enumerate(native["redshift"]):
            ratio = cocoa[key][row] / native[key][row] - 1
            axis.semilogx(native[grid], 100*ratio, label=f"z={z:g}")
        axis.axhline(0, color="0.5", linewidth=0.7)
        axis.set(title=title, ylabel="CoCoA / TJPCov-CCL − 1 [%]",
                 xlabel=r"$M\ [M_\odot/h]$" if grid == "mass"
                 else r"$k\ [h/{\rm Mpc}]$")
    axes[0, 0].legend()
    axis = axes[-1, -1]
    for row, z in enumerate(native["redshift"]):
        delta = np.abs(cocoa["profile_matched_c"][row] - native["profile"][row])
        axis.semilogx(native["k"], np.max(delta, axis=1), label=f"z={z:g}")
    axis.set(title="NFW: same Duffy concentration",
             xlabel=r"$k\ [h/{\rm Mpc}]$", ylabel="Maximum |Δu| across masses")
    fig.suptitle("Native halo ingredients · shared LSST Y1 CAMB inputs")
    args.figures.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(args.figures / f"halo_ingredients.{suffix}", dpi=180,
                    facecolor="white")
    plt.close(fig)
    print(args.output / "comparison.json")


if __name__ == "__main__":
    main()
