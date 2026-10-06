"""Compare all entries of two native five-multipole shear SSC matrices.

The halo responses and footprint prescriptions differ between codes.
This comparison measures that combined difference, not equality under
shared responses. Input spectra, source identity and ell must still match.
"""

import argparse
from pathlib import Path

from common import load_bundle, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cocoa", type=Path)
    parser.add_argument("tjpcov", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a fresh comparison output directory")
    records = {
        "cocoa": load_bundle(args.cocoa, "cocoa-ssc-shear-v1"),
        "tjpcov": load_bundle(args.tjpcov, "tjpcov-ssc-shear-v1"),
    }
    for key in ("input_manifest_sha256", "gaussian_manifest_sha256"):
        if records["cocoa"][key] != records["tjpcov"][key]:
            raise ValueError(f"The two SSC runs have different {key}")

    import numpy as np
    from scipy.linalg import eigvalsh

    left = np.load(args.cocoa / "ssc.npz", allow_pickle=False)
    right = np.load(args.tjpcov / "ssc.npz", allow_pickle=False)
    np.testing.assert_array_equal(left["ell"], right["ell"])
    matrices = {"cocoa": left["ssc"], "tjpcov": right["ssc"]}
    for name, value in matrices.items():
        if (value.shape != (5, 5) or not np.all(np.isfinite(value))
                or np.any(np.diag(value) <= 0)):
            raise ValueError(f"Invalid {name} native SSC matrix")
    rms = np.sqrt(np.diag(matrices["cocoa"]))
    scale = rms[:, None] * rms[None, :]
    difference = (matrices["tjpcov"] - matrices["cocoa"]) / scale
    # Preserve every raw entry. Eigenvalues describe the symmetric part;
    # record the original antisymmetric residual instead of hiding it.
    modes = {}
    symmetric = {}
    for name, values in matrices.items():
        symmetric[name] = (values + values.T) / 2
        modes[name] = {
            "max_asymmetry_over_cocoa_rms": float(np.max(
                np.abs(values-values.T)/scale)),
            "min_eigenvalue_over_cocoa_rms": float(eigvalsh(
                symmetric[name]/scale)[0]),
        }
    ratios = None
    if all(item["min_eigenvalue_over_cocoa_rms"] > 0 for item in modes.values()):
        values = eigvalsh(symmetric["tjpcov"], symmetric["cocoa"])
        ratios = [float(values[0]), float(values[-1])]
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output / "matrices.npz", ell=left["ell"], **matrices)
    record = {
        "scope": "Native shear SSC at matching effective multipoles",
        "difference": "TJPCov minus CoCoA; CoCoA SSC rms normalization",
        "max_variance_scaled_difference": float(np.max(np.abs(difference))),
        "diagonal_fractional_difference": (
            np.diag(matrices["tjpcov"])/np.diag(matrices["cocoa"])-1).tolist(),
        "generalized_variance_ratio_range": ratios, "modes": modes,
        "runs": records, "script_sha256": sha256(__file__),
        "files": {"matrices.npz": sha256(args.output / "matrices.npz")},
        "interpretation": "Different halo responses and cap/disc windows; "
                          "numerical refinement must be tested separately.",
    }
    write_json(args.output / "comparison.json", record)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    args.figures.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
    for ax, name in zip(axes[:2], ("cocoa", "tjpcov")):
        value = matrices[name]
        local = np.sqrt(np.diag(value))
        image = ax.imshow(value/local[:, None]/local[None, :], origin="lower",
                          vmin=-1, vmax=1, cmap="RdBu_r")
        ax.set_title({"cocoa": "CoCoA", "tjpcov": "TJPCov"}[name])
    fig.colorbar(image, ax=list(axes[:2]), shrink=0.8, label="SSC correlation")
    extent = float(np.max(np.abs(difference)))
    image = axes[2].imshow(100*difference, origin="lower", cmap="RdBu_r",
                           vmin=-100*extent, vmax=100*extent)
    axes[2].set_title("TJPCov − CoCoA")
    fig.colorbar(image, ax=axes[2], shrink=0.8,
                 label="Difference / CoCoA SSC rms product [%]")
    for ax in axes:
        ax.set(xlabel="Increasing effective multipole",
               ylabel="Increasing effective multipole")
    fig.suptitle("LSST Y1 source bin 3 · native SSC · distinct halo/window models")
    for suffix in ("png", "pdf"):
        fig.savefig(args.figures / f"ssc_matrices.{suffix}", dpi=180,
                    facecolor="white")
    plt.close(fig)
    print(f"Saved all SSC entries and diagnostics in {args.output}")


if __name__ == "__main__":
    main()
