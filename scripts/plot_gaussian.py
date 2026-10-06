"""Plot measured TJPCov/CoCoA Gaussian matrices in the OneCov figure style.

Only a passing comparison is displayed. This creates scientific figures,
not timing plots; native SSC/cNG are absent from this Gaussian test.
"""

import argparse
import json
from pathlib import Path

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads((args.comparison / "comparison.json").read_text())
    if not report["passed"]:
        raise ValueError("Review the failed comparison before making its figure")
    matrix_file = args.comparison / "matrices.npz"
    if sha256(matrix_file) != report["files"]["matrices.npz"]:
        raise ValueError("Comparison matrices changed since validation")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    data = np.load(matrix_file, allow_pickle=False)
    cocoa, tjpcov = data["cocoa_total"], data["tjpcov_total"]
    rms = np.sqrt(np.diag(cocoa))
    residual = (tjpcov - cocoa) / rms[:, None] / rms[None, :]
    ell = data["effective_ell"]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    figure, axes = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
    for axis, name, matrix in zip(axes[:2], ("CoCoA", "TJPCov"), (cocoa, tjpcov)):
        scale = np.sqrt(np.diag(matrix))
        image = axis.imshow(matrix / scale[:, None] / scale[None, :],
                            origin="lower", vmin=-1, vmax=1, cmap="RdBu_r")
        axis.set_title(name)
    figure.colorbar(image, ax=list(axes[:2]), shrink=0.8, label="Correlation")
    maximum = float(np.max(np.abs(residual)))
    if maximum == 0:
        axes[2].text(0.5, 0.5, "Exactly zero difference", ha="center",
                     va="center", transform=axes[2].transAxes)
    else:
        image = axes[2].imshow(residual, origin="lower", cmap="RdBu_r",
                                vmin=-maximum, vmax=maximum)
        figure.colorbar(image, ax=axes[2], shrink=0.8,
                        label=r"$\Delta C_{ij}/\sqrt{C^{\rm CoCoA}_{ii}C^{\rm CoCoA}_{jj}}$")
    axes[2].set_title("TJPCov − CoCoA")
    for axis in axes:
        axis.set_xlabel("Increasing multipole band")
        axis.set_ylabel("Increasing multipole band")
    figure.suptitle("LSST Y1 source bin · shared spectra · Gaussian covariance")

    fractions, axis = plt.subplots(figsize=(7, 4), layout="constrained")
    for name, label, color in (("sample_variance", "Sample variance", "#2563a6"),
                                ("mixed", "Signal × noise", "#c05621"),
                                ("noise", "Pure noise", "#23836c")):
        axis.plot(ell, np.diag(data[f"cocoa_{name}"]) / np.diag(cocoa),
                   color=color, label=label)
        axis.plot(ell, np.diag(data[f"tjpcov_{name}"]) / np.diag(tjpcov),
                   linestyle="none", marker="o", fillstyle="none", color=color)
    axis.set(xlabel=r"Effective multipole $\ell$",
             ylabel="Fraction of Gaussian variance", ylim=(-0.03, 1.05),
             title="CoCoA lines / TJPCov markers")
    axis.legend()
    args.output.mkdir(parents=True, exist_ok=True)
    for name, item in (("gaussian_matrices", figure),
                        ("gaussian_components", fractions)):
        for suffix in ("png", "pdf"):
            item.savefig(args.output / f"{name}.{suffix}", dpi=180,
                          facecolor="white")
        plt.close(item)
    print(f"Saved Gaussian figures in {args.output}")


if __name__ == "__main__":
    main()
