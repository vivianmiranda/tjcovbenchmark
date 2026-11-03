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

The band-limited route is what TJPCov's real-space calculator runs at
the studied checkout: covariance_gaussian_fsky.py multiplies
(C_ell + N) into both Wick pairings, pure-noise product included, and
covariance_builder.py projects the result up to a caller-chosen lmax
and bins it. CoCoA replaces the pure-noise product by the exact
pair-count term (gaussian_noise_pair_cov in cosmolike_core). The
figure labels name the two routes accordingly.

Three separate figures are written: the unbinned projection at one bin
center (noise_center), the bin-averaged projection for three bins
(noise_binned), and one matrix row at fixed band limits (noise_row).
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
EXACT_LABEL = "exact pair count (CoCoA)"
PROJECTED_LABEL = "TJPCov/CCL route"


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


def math_band_limit(limit):
    """A band limit as mathtext, with a grouped thousands separator."""
    grouped = f"{limit:,}".replace(",", "{,}")
    return rf"$\ell_{{\max}} = {grouped}$"


def save(figure, output, stem):
    for extension in ("png", "pdf"):
        figure.savefig(output / f"{stem}.{extension}", dpi=200)


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

    plt.rcParams.update({"font.size": 13, "axes.spines.top": False,
                         "axes.spines.right": False})
    ells = np.arange(LMAX + 1)

    def arcmin(value):
        scaled = np.rad2deg(value) * 60.0
        return f"{scaled:.1f}" if scaled < 10.0 else f"{scaled:.0f}"

    # Figure 1: the projection with kernels at one bin center.
    figure, axis = plt.subplots(figsize=(7.0, 4.6), layout="constrained")
    axis.plot(ells[1:], point_sum[1:], color="C3",
              label="projected at the bin center")
    axis.axhline(1.0, color="black", linewidth=0.9, linestyle="dashed",
                 label=EXACT_LABEL)
    axis.set_xscale("log")
    axis.set_title("Unbinned projection of the pure noise: never settles")
    axis.set_xlabel(r"band limit $\ell_{\max}$")
    axis.set_ylabel("partial sum / exact diagonal")
    axis.legend(loc="upper left", frameon=False)
    save(figure, args.output, "noise_center")
    plt.close(figure)

    # Figure 2: the projection with bin-averaged kernels, three bins.
    # Curve names sit next to the curves in their own color; no legend.
    figure, axis = plt.subplots(figsize=(7.0, 4.6), layout="constrained")
    for index, color in zip(SHOWN_BINS, ("C0", "C1", "C2")):
        axis.plot(ells[1:], diagonal[index][1:], color=color)
        axis.axvline(resolve[index], color=color, linewidth=0.9,
                     linestyle="dotted")
    axis.axhline(1.0, color="black", linewidth=0.9, linestyle="dashed")
    axis.text(1.3, 1.045, EXACT_LABEL, color="black")
    axis.text(1.3, 0.84,
              f"{arcmin(edges[25])} to {arcmin(edges[26])} arcmin",
              color="C2")
    axis.text(60.0, 0.62,
              f"{arcmin(edges[12])} to {arcmin(edges[13])} arcmin",
              color="C1")
    axis.text(1.1e3, 0.42,
              f"{arcmin(edges[0])} to {arcmin(edges[1])} arcmin",
              color="C0")
    axis.set_xscale("log")
    axis.set_ylim(0.0, 1.12)
    axis.set_title("Bin-averaged projection: the TJPCov/CCL route")
    axis.set_xlabel(r"band limit $\ell_{\max}$")
    axis.set_ylabel("partial sum / exact diagonal")
    save(figure, args.output, "noise_binned")
    plt.close(figure)

    # Figure 3: one matrix row of the pure-noise term at fixed limits.
    figure, axis = plt.subplots(figsize=(7.0, 4.6), layout="constrained")
    centers = np.rad2deg(np.sqrt(edges[:-1] * edges[1:])) * 60.0
    exact = np.zeros(NBIN)
    exact[0] = 1.0
    axis.plot(centers, exact, color="black", marker="o", linewidth=0.9,
              label=EXACT_LABEL)
    for limit, color in zip(ROW_BAND_LIMITS, ("C4", "C9")):
        axis.plot(centers, rows[limit], marker="s", markersize=4,
                  color=color,
                  label=f"{PROJECTED_LABEL}, {math_band_limit(limit)}")
    axis.set_xscale("log")
    axis.set_title("First matrix row of the pure-noise term")
    axis.set_xlabel(r"$\theta_j$ [arcmin]")
    axis.set_ylabel("entry / exact diagonal")
    axis.legend(loc="upper right", frameon=False)
    save(figure, args.output, "noise_row")
    plt.close(figure)
    print(f"figures written to {args.output}")


if __name__ == "__main__":
    main()
