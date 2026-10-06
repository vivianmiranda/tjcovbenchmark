"""Collect verified trispectrum diagnostics and plot native/model differences.

Every pair is retained. Native ratios use CCL as denominator; the smooth
power control holds CoCoA's halo moments fixed. Timings are not plotted.
"""

import argparse
import json
import shutil
from pathlib import Path

from common import load_bundle, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.figures.exists():
        parser.error("Choose new result and figure directories")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    runs = {}
    manifests = {}
    for label in ("native_m128", "native_m255", "direct_m255",
                  "bhattacharya_m255", "tinker10_m255", "both_m255",
                  "cocoa_i0", "cocoa_i1", "cocoa_i2"):
        folder = args.work / f"trispectrum_{label}"
        schema = "cocoa" if label.startswith("cocoa") else "tjpcov"
        manifests[label] = load_bundle(folder, f"{schema}-trispectrum-v1")
        runs[label] = dict(np.load(folder / "trispectrum.npz", allow_pickle=False))
    for label in ("native", "direct"):
        folder = args.work / f"trispectrum_comparison_{label}"
        manifests[f"comparison_{label}"] = load_bundle(
            folder, "trispectrum-comparison-v1")
        runs[f"comparison_{label}"] = dict(np.load(
            folder / "comparison.npz", allow_pickle=False))
    power_folder = args.work / "trispectrum_power_refinement"
    power_record = json.loads((power_folder / "report.json").read_text())
    assert sha256(power_folder / "diagnostic.npz") == power_record["files"]["diagnostic.npz"]
    power = dict(np.load(power_folder / "diagnostic.npz", allow_pickle=False))
    reference = runs["native_m255"]
    co = runs["cocoa_i0"]
    cocoa_hash = sha256(args.work / "trispectrum_cocoa_i0/manifest.json")
    if (power_record["status"] != "completed"
            or power_record["cocoa_manifest_sha256"] != cocoa_hash
            or power_record["input_manifest_sha256"]
            != manifests["cocoa_i0"]["input_manifest_sha256"]):
        raise ValueError("Power diagnostic uses a different native baseline")
    for label, parent in (("native", "native_m255"), ("direct", "direct_m255")):
        record = manifests[f"comparison_{label}"]
        if (record["cocoa_manifest_sha256"] != cocoa_hash
                or record["native_manifest_sha256"] != sha256(
                    args.work / f"trispectrum_{parent}/manifest.json")):
            raise ValueError("Comparison source manifests do not match")
    first, second, k = co["first"], co["second"], co["k"]
    diagonal = first == second
    for label, data in {**runs, "power": power}.items():
        for key in ("k", "redshift", "first", "second"):
            np.testing.assert_array_equal(data[key], co[key])
        if not all(np.isfinite(v).all() for v in data.values()):
            raise ValueError(f"Non-finite result: {label}")
    np.testing.assert_allclose(power["terms"][0], co["terms"], rtol=1e-12, atol=0)
    inputs = {r.get("input_manifest_sha256", r.get("cocoa", {}).get(
        "input_manifest_sha256")) for r in manifests.values()}
    if len(inputs) != 1:
        raise ValueError("Different input bundles across the diagnostics")
    args.output.mkdir(parents=True)
    args.figures.mkdir(parents=True)

    # The archive includes full triangular arrays, not just plotted slices.
    arrays = {f"{label}__{key}": value for label, data in runs.items()
              for key, value in data.items()}
    arrays.update({f"power__{key}": value for key, value in power.items()})
    np.savez_compressed(args.output / "trispectra.npz", **arrays)
    names = ["1h", "2h (1+3)", "2h (2+2)", "3h", "4h"]
    statistics = {
        "native": manifests["comparison_native"]["native"],
        "shared_inputs_direct_growth": manifests["comparison_direct"]["common_power_and_moments"],
        "power_reader_effect": manifests["comparison_direct"]["power_reader_effect"],
        "cocoa_integration_refinement": [],
        "ccl_mass_refinement": np.max(np.abs(runs["native_m128"]["terms"]
            /reference["terms"]-1), axis=(0, 2)).tolist(),
        "native_diagnostic_changes": {},
        "power_refinement": power_record,
    }
    for left, right in (("cocoa_i0", "cocoa_i1"), ("cocoa_i1", "cocoa_i2")):
        statistics["cocoa_integration_refinement"].append(np.max(
            np.abs(runs[left]["terms"]/runs[right]["terms"]-1), axis=(0, 2)).tolist())
    for label in ("direct_m255", "bhattacharya_m255", "tinker10_m255", "both_m255"):
        statistics["native_diagnostic_changes"][label] = np.max(
            np.abs(runs[label]["terms"]/reference["terms"]-1), axis=(0, 2)).tolist()

    def save(figure, name):
        for suffix in ("png", "pdf"):
            figure.savefig(args.figures / f"{name}.{suffix}", dpi=170)
        plt.close(figure)

    def halo_orders(terms):
        return np.stack((terms[:, 0], terms[:, 1]+terms[:, 2],
                         terms[:, 3], terms[:, 4]), axis=1)

    cocoa_orders, ccl_orders = map(halo_orders, (co["terms"], reference["terms"]))
    if np.any(cocoa_orders[..., diagonal] <= 0) or np.any(ccl_orders[..., diagonal] <= 0):
        raise ValueError("Logarithmic amplitude panels require positive terms")
    colors = ["#31688e", "#d95f02", "#238b45"]
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), layout="constrained")
    for term, name in enumerate(("1 halo", "2 halo", "3 halo", "4 halo")):
        for iz, (z, color) in enumerate(zip(co["redshift"], colors)):
            axes[0, term].loglog(k, cocoa_orders[iz, term, diagonal],
                                color=color, label=f"z={z:g}")
            axes[0, term].loglog(k, ccl_orders[iz, term, diagonal],
                                color=color, linestyle="--", marker=".")
            axes[1, term].semilogx(k, 100*(cocoa_orders[iz, term, diagonal]
                /ccl_orders[iz, term, diagonal]-1), color=color)
        axes[0, term].set(title=name, ylabel=r"$T(k,k)$ [(Mpc/h)$^9$]")
        axes[1, term].axhline(0, color="0.6", lw=0.7)
        axes[1, term].set(xlabel=r"$k$ [h/Mpc]", ylabel="CoCoA / CCL − 1 [%]")
    axes[0, 0].legend()
    fig.suptitle("Matter trispectrum, equal k pairs · solid: CoCoA; dashed: TJPCov/CCL\n"
                 "Native halo prescriptions, same CAMB inputs; no survey projection")
    save(fig, "trispectrum_terms")

    # Unequal pairs expose squeezed configurations absent from diagonal plots.
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), layout="constrained")
    iz = 1
    for term, axis in enumerate(axes.ravel()):
        actual = co["terms"][iz, term] if term < 5 else co["terms"][iz].sum(0)
        expected = reference["terms"][iz, term] if term < 5 else reference["terms"][iz].sum(0)
        values = 100*(actual/expected-1)
        matrix = np.empty((len(k), len(k)))
        matrix[first, second] = values
        matrix[second, first] = values
        bound = np.max(np.abs(matrix))
        picture = axis.imshow(matrix, origin="lower", cmap="RdBu_r",
                              vmin=-bound, vmax=bound, interpolation="nearest",
                              extent=(-3.25, 1.25, -3.25, 1.25))
        axis.set(title=names[term] if term < 5 else "Sum of all halo orders",
                 xlabel=r"log$_{10}$ K [h/Mpc]", ylabel=r"log$_{10}$ Q [h/Mpc]")
        fig.colorbar(picture, ax=axis, label="CoCoA / CCL − 1 [%]", shrink=0.83)
    fig.suptitle("Every wavenumber pair at z=0.5 · native matter trispectra\n"
                 "Each panel has its own scale; this is T(K,Q), not a covariance matrix")
    save(fig, "trispectrum_pairs")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    pair = np.flatnonzero((first == 0) & (second == 5))[0]
    for iz, (z, color) in enumerate(zip(co["redshift"], colors)):
        residual = 100*(power["terms"][:, iz, 4, pair]/power["smooth"][iz, 4, pair]-1)
        axes[0].semilogx(power["node_counts"], residual, "o-", color=color, label=f"z={z:g}")
        axes[1].semilogx(k, 100*(co["power"][iz]/reference["power"][iz]-1),
                         ".-", color=color)
        shared = runs["comparison_direct"]["shared"][iz, 4]
        ref = runs["direct_m255"]["terms"][iz, 4]
        axes[2].plot(np.arange(len(first)), 1e6*(shared/ref-1), ".", color=color)
    for axis in axes:
        axis.axhline(0, color="0.6", lw=0.7)
    axes[0].set(title="Fixed moments; denser power table",
                xlabel="Power k-grid nodes", ylabel="4h / smooth-input control − 1 [%]")
    from matplotlib.ticker import NullLocator
    axes[0].set_xticks(power["node_counts"], ["1.5k", "3k", "6k", "12k", "24k"])
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].legend()
    axes[1].set(title="Small native power differences", xlabel=r"$k$ [h/Mpc]",
                ylabel="CoCoA / CCL linear power − 1 [%]")
    axes[2].set(title="Same P samples and halo moments", xlabel="Unordered k-pair index",
                ylabel="CoCoA / CCL 4h − 1 [parts per million]")
    fig.suptitle("Tracing the 4h discrepancy · K=0.001, Q=0.316 h/Mpc in left panel\n"
                 "Cubic table fill + actual linear reader; CCL direct-growth control at right")
    save(fig, "trispectrum_power_diagnostic")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    iz = 1
    styles = (("native_m255", "Native CCL"),
              ("bhattacharya_m255", "CCL: concentration changed"),
              ("tinker10_m255", "CCL: abundance changed"),
              ("both_m255", "CCL: both changed"))
    for label, legend in styles:
        values = runs[label]["terms"][iz]
        for axis, term in zip(axes, (0, 1, 2)):
            axis.semilogx(k, 100*(co["terms"][iz, term, diagonal]
                /values[term, diagonal]-1), ".-", label=legend)
    for axis, name in zip(axes, names[:3]):
        axis.axhline(0, color="0.6", lw=0.7)
        axis.set(title=name, xlabel=r"$k$ [h/Mpc]", ylabel="CoCoA / selected CCL model − 1 [%]")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Public CCL fit switches at z={co['redshift'][iz]:g}, equal k pairs\n"
                 "Diagnostic models: matching fit names does not match normalization conventions")
    save(fig, "trispectrum_model_diagnostic")
    report = dict(schema="trispectrum-study-v1", halo_order=names,
                  manifests=manifests, statistics=statistics,
                  script_sha256=sha256(__file__),
                  files={"trispectra.npz": sha256(args.output / "trispectra.npz")},
                  figures={p.name: sha256(p) for p in args.figures.iterdir()})
    write_json(args.output / "comparison.json", report)
    shutil.copy2(power_folder / "report.json", args.output / "power_refinement.json")


if __name__ == "__main__":
    main()
