"""Show why the pure-noise covariance resists a band-limited projection.

Pure numerics on Legendre kernels: no survey data and no Boltzmann code.
The white-noise part of a harmonic covariance is constant in multipole,
so its projection to two angular bins i and j must reproduce the
bin-projected Dirac delta. With kernels averaged over each bin, the
normalized partial sum

    T_ij(L) = sum_{ell<=L} (ell + 1/2) Pbar_ell(i) Pbar_ell(j)
              sqrt(dx_i dx_j)

approaches the identity matrix as the band limit L grows (dx is the bin
width in cos theta). With kernels evaluated at the bin center the same
sum has no limit. The spin-2 kernels of xi_pm obey the same
completeness relation, so the convergence behavior shown here carries
over; only the exact kernels differ.
"""

import argparse
from pathlib import Path

import numpy as np


# The study's real-space scope: LSST Y1 separations, logarithmic bins.
THETA_MIN_ARCMIN = 2.5
THETA_MAX_ARCMIN = 900.0
NBIN = 26
LMAX = 100_000
ROW_BAND_LIMITS = (2_000, 30_000)
SHOWN_BINS = (0, 12, 25)


def legendre_at(x, lmax):
    """Legendre P_0..P_lmax at each x, by the three-term recurrence."""
    table = np.empty((lmax + 1, x.size))
    table[0] = 1.0
    table[1] = x
    for ell in range(1, lmax):
        table[ell + 1] = ((2*ell + 1)*x*table[ell] - ell*table[ell - 1]) \
                         / (ell + 1)
    return table


def bin_averaged_kernels(edges_rad, lmax):
    """Pbar_ell for each bin, plus the bin widths dx in cos theta.

    The antiderivative of P_ell is (P_{ell+1} - P_{ell-1})/(2 ell + 1),
    so one Legendre table at the bin edges gives every averaged kernel.
    """
    x = np.cos(edges_rad)
    table = legendre_at(x, lmax + 1)
    dx = x[:-1] - x[1:]
    kernels = np.empty((lmax + 1, NBIN))
    kernels[0] = 1.0
    ells = np.arange(1, lmax + 1)
    anti = (table[2:] - table[:-2]) / (2*ells + 1)[:, None]
    kernels[1:] = (anti[:, :-1] - anti[:, 1:]) / dx
    return kernels, dx


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    edges = np.deg2rad(np.geomspace(THETA_MIN_ARCMIN, THETA_MAX_ARCMIN,
                                    NBIN + 1) / 60.0)
    kernels, dx = bin_averaged_kernels(edges, LMAX)
    weight = np.arange(LMAX + 1) + 0.5
    widths_rad = np.diff(edges)
    resolve = np.pi / widths_rad

    # Diagonal partial sums, normalized so the exact answer is one.
    diagonal = {}
    for index in SHOWN_BINS:
        diagonal[index] = np.cumsum(weight * kernels[:, index]**2) * dx[index]

    # The same sum with the kernel evaluated at the narrowest bin's
    # center instead of averaged over it.
    center = np.array([0.5*(edges[0] + edges[1])])
    point = legendre_at(np.cos(center), LMAX)[:, 0]
    point_sum = np.cumsum(weight * point**2) * dx[0]

    # One full matrix row for the narrowest bin at two band limits.
    norm = np.sqrt(dx[0] * dx)
    row_terms = weight[:, None] * kernels[:, :1] * kernels * norm
    rows = {limit: row_terms[:limit + 1].sum(axis=0)
            for limit in ROW_BAND_LIMITS}

    for index in SHOWN_BINS:
        print(f"bin {index + 1}: {np.rad2deg(edges[index])*60:.1f} to "
              f"{np.rad2deg(edges[index + 1])*60:.1f} arcmin, "
              f"pi/width = {resolve[index]:.0f}, "
              f"recovered at L={LMAX}: {diagonal[index][-1]:.4f}")
    print(f"point-evaluated sum at L={LMAX}: {point_sum[-1]:.2f}")
    crossing = np.argmax(point_sum > 1.0)
    print(f"point-evaluated sum first exceeds 1 at L={crossing}")
    for limit in ROW_BAND_LIMITS:
        print(f"row diagonal at L={limit}: {rows[limit][0]:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False})
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.2),
                                layout="constrained")
    ells = np.arange(LMAX + 1)

    figure.suptitle("Pure-noise diagonal: band-limited projection (color) "
                    "against the exact pair-count answer (black)")
    axes[0].plot(ells[1:], point_sum[1:], color="C3",
                 label="band-limited projection")
    axes[0].axhline(1.0, color="black", linewidth=0.8, linestyle="dashed",
                    label="exact pair-count answer")
    axes[0].set_xscale("log")
    axes[0].set_title("Kernel at the bin center: never settles")
    axes[0].set_ylabel("partial sum / exact diagonal")
    axes[0].legend(loc="upper left", fontsize=10)

    def arcmin(value):
        scaled = np.rad2deg(value) * 60.0
        return f"{scaled:.1f}" if scaled < 10.0 else f"{scaled:.0f}"

    labels = [f"{arcmin(edges[i])} to {arcmin(edges[i + 1])} arcmin"
              for i in SHOWN_BINS]
    for (index, series), label, color in zip(diagonal.items(), labels,
                                             ("C0", "C1", "C2")):
        axes[1].plot(ells[1:], series[1:], color=color, label=label)
        axes[1].axvline(resolve[index], color=color, linewidth=0.8,
                        linestyle="dotted")
    axes[1].axhline(1.0, color="black", linewidth=0.8, linestyle="dashed")
    axes[1].set_xscale("log")
    axes[1].set_ylim(0.0, 1.1)
    axes[1].set_title("Bin-averaged kernels: converges, slowly")
    axes[1].legend(loc="upper left", fontsize=10, title="bin, with "
                   r"$\pi$/width dotted", title_fontsize=10)

    centers = np.rad2deg(np.sqrt(edges[:-1] * edges[1:])) * 60.0
    exact = np.zeros(NBIN)
    exact[0] = 1.0
    axes[2].plot(centers, exact, color="black", marker="o", linewidth=0.8,
                 label="exact: CoCoA's analytic term")
    for limit, color in zip(ROW_BAND_LIMITS, ("C4", "C9")):
        axes[2].plot(centers, rows[limit], marker="s", markersize=4,
                     color=color,
                     label=rf"projected, $\ell_{{\max}} = {limit:,}$")
    axes[2].set_xscale("log")
    axes[2].set_title("First matrix row: projection against exact")
    axes[2].set_xlabel(r"$\theta_j$ [arcmin]")
    axes[2].set_ylabel("entry / exact diagonal")
    axes[2].legend(loc="upper right", fontsize=10)

    for axis in axes[:2]:
        axis.set_xlabel(r"band limit $\ell_{\max}$")
    for extension in ("png", "pdf"):
        figure.savefig(args.output / f"noise_projection.{extension}",
                       dpi=200)
    print(f"figure written to {args.output}")


if __name__ == "__main__":
    main()
