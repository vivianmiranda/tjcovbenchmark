"""Archive and plot saved SSC model diagnostics without running either model.

The two supplied CoCoA diagnostic folders and their parent CCL exports
must be complete and hash-consistent. All matrix entries are retained.
The plots distinguish native-model substitutions from response-formula
tests that use common CCL ingredients. Both output directories must be new.
"""

import argparse
import shutil
from pathlib import Path

from common import load_bundle, sha256, write_json


def read_case(path, parent):
    """Check both calculation stages before accepting their saved arrays."""
    import numpy as np

    record = load_bundle(path, "cocoa-tjpcov-ssc-model-diagnostic-v1")
    ccl = load_bundle(parent, "tjpcov-ssc-model-export-v1")
    if record["status"] != "completed" or ccl["status"] != "completed":
        raise ValueError("Both stages must have completed")
    if record["ccl_export_manifest_sha256"] != sha256(parent / "manifest.json"):
        raise ValueError("The CoCoA diagnostic used a different CCL export")
    if record["cocoa_manifest_sha256"] != ccl["cocoa_manifest_sha256"]:
        raise ValueError("The two stages used different native CoCoA inputs")
    if (ccl["baseline_gate"][
            "max_abs_difference_over_reference_rms_product"] > 1e-8
            or record["formula_gate_max_fractional_difference"] > 1e-10
            or record["cases"]["cocoa_response_cap_window"][
                "max_abs_difference_over_reference_rms_product"] > 1e-10):
        raise ValueError("A required reproduction gate did not pass")

    with np.load(path / "models.npz", allow_pickle=False) as saved:
        data = {name: saved[name] for name in saved.files}
    if any(not np.all(np.isfinite(value)) for value in data.values()):
        raise ValueError("Non-finite saved diagnostic values")
    if data["ell"].shape != (5,) or data["cocoa_native_ssc"].shape != (5, 5):
        raise ValueError("Expected the complete five-multipole SSC pilot")
    if np.any(np.diag(data["cocoa_native_ssc"]) <= 0):
        raise ValueError("Reference SSC variances must be positive")
    if np.any(data["a"] <= 0) or np.any(data["a"] >= 1):
        raise ValueError("Line-of-sight plots require positive redshifts")
    return record, ccl, data


def save_figure(fig, directory, name):
    """Keep matching raster and vector versions of each scientific figure."""
    for extension in ("png", "pdf"):
        fig.savefig(directory / f"{name}.{extension}", dpi=180,
                    facecolor="white")


def matrix_figure(cases, output):
    """Show each complete substitution, with one scale across six panels."""
    import numpy as np
    import matplotlib.pyplot as plt

    choices = (
        ("cocoa_response_disc_window", "Replace background window only"),
        ("ccl_response_cap_window", "Replace matter response only"),
        ("ccl_response_disc_window", "Replace both"),
    )
    changes = []
    for _, _, data in cases:
        baseline = data["cocoa_native_ssc"]
        rms = np.sqrt(np.diag(baseline))
        changes.append([100*(data[key]-baseline)/rms[:, None]/rms[None, :]
                        for key, _ in choices])
    limit = max(float(np.max(np.abs(value)))
                for row in changes for value in row)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout="constrained")
    for row, (_, _, data) in enumerate(cases):
        labels = [f"{value:g}" for value in data["ell"]]
        for column, (_, title) in enumerate(choices):
            ax = axes[row, column]
            values = changes[row][column]
            picture = ax.imshow(values, origin="lower", cmap="RdBu_r",
                                vmin=-limit, vmax=limit)
            ax.set(xticks=range(5), yticks=range(5), xticklabels=labels,
                   yticklabels=labels, xlabel=r"Effective multipole $\ell_j$",
                   ylabel=r"Effective multipole $\ell_i$")
            ax.set_title(title if row == 0 else "High multipoles")
            if column == 0:
                ax.set_ylabel(("Low multipoles\n" if row == 0 else "")
                              + r"Effective multipole $\ell_i$")
            for i in range(5):
                for j in range(5):
                    ax.text(j, i, f"{values[i, j]:+.1f}", ha="center",
                            va="center", fontsize=9,
                            color="white" if abs(values[i, j]) > .6*limit
                            else "black")
    bar = fig.colorbar(picture, ax=axes, shrink=.85, pad=.025)
    bar.set_label(r"$100\,\Delta C_{ij}/\sqrt{C^{\rm CoCoA}_{ii}"
                  r"C^{\rm CoCoA}_{jj}}$ [%]")
    fig.suptitle("LSST Y1 source bin 3 · SSC only · all 25 entries per panel\n"
                 "CCL substitutions with CoCoA geometry and tracer weights "
                 "held fixed", fontsize=14)
    save_figure(fig, output, "ssc_model_swaps")
    plt.close(fig)


def ratio_percent(numerator, denominator, name):
    """Mask only numerically vanishing denominators, never response signs."""
    import numpy as np

    floor = 1e-12*float(np.max(np.abs(denominator)))
    valid = np.abs(denominator) > floor
    ratio = np.full(denominator.shape, np.nan)
    ratio[valid] = 100*(numerator[valid]/denominator[valid]-1)
    note = {"quantity": name, "denominator_floor": floor,
            "relative_floor": 1e-12,
            "masked_samples": int(np.count_nonzero(~valid)),
            "total_samples": int(valid.size)}
    return ratio, note


def sightline_figure(cases, output):
    """Show the response/window trends and where SSC accumulates radially."""
    import numpy as np
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    notes = []
    for column, (_, _, data) in enumerate(cases):
        redshift = 1/data["a"]-1
        order = np.argsort(redshift)
        shown = order[(redshift[order] >= .01) & (redshift[order] <= 2)]
        ax = axes[0, column]
        delta, note = ratio_percent(
            data["ccl_response_at_cocoa_k_Mpc3"],
            data["cocoa_response_Mpc3"], "matter response [Mpc^3]")
        kmin, kmax = np.exp(data["native_response_lk"][[0, -1]])
        wave = data["k_cocoa_Mpc_inverse"]
        inside = (wave >= kmin) & (wave <= kmax)
        for i, ell in enumerate(data["ell"]):
            color = f"C{i}"
            ax.plot(redshift[shown],
                    np.where(inside[:, i], delta[:, i], np.nan)[shown],
                    color=color, label=rf"$\ell={ell:g}$")
            # Dotted segments identify any native k extrapolation in the
            # displayed redshift interval without discarding its signs.
            ax.plot(redshift[shown],
                    np.where(~inside[:, i], delta[:, i], np.nan)[shown],
                    color=color, linestyle=":")
        ax.axhline(0, color=".35", linewidth=.7)
        ax.set(xlabel="Redshift z",
               ylabel="CCL / CoCoA response − 1 [%]",
               title=("Low" if column == 0 else "High")
                     + " multipoles: same k at each z")
        ax.legend(fontsize=8, ncol=2)
        note.update({"redshift_range": [float(redshift.min()),
                                         float(redshift.max())],
                     "ccl_tabulated_k_range_Mpc_inverse": [kmin, kmax],
                     "displayed_redshift_range": [.01, 2],
                     "extrapolated_samples": int(np.count_nonzero(~inside)),
                     "displayed_extrapolated_samples": int(
                         np.count_nonzero(~inside[shown])),
                     "negative_ccl_responses": int(np.count_nonzero(
                         data["ccl_response_at_cocoa_k_Mpc3"] < 0))})
        notes.append(note)

    # The survey-window predictions do not depend on the chosen multipoles.
    # Require their equality before showing one shared window comparison.
    first, second = cases[0][2], cases[1][2]
    for key in ("a", "cocoa_cap_variance_Mpc", "ccl_disc_at_cocoa_a_Mpc"):
        np.testing.assert_array_equal(first[key], second[key])
    redshift = 1/first["a"]-1
    order = np.argsort(redshift)
    shown = order[(redshift[order] >= .01) & (redshift[order] <= 2)]
    delta, note = ratio_percent(first["ccl_disc_at_cocoa_a_Mpc"],
                                first["cocoa_cap_variance_Mpc"],
                                "background variance [Mpc]")
    axes[1, 0].plot(redshift[shown], delta[shown], color="C1")
    axes[1, 0].axhline(0, color=".35", linewidth=.7)
    axes[1, 0].set(xlabel="Redshift z",
                   ylabel="CCL / CoCoA background variance − 1 [%]",
                   title="Window difference before radial weighting")
    notes.append(note)

    # Sum the recorded quadrature contributions for one diagonal in each
    # multipole set. Their physical length conversion cancels in the
    # cumulative fraction; no new halo or projection calculation is made.
    for name, (_, _, data) in zip(("Low", "High"), cases):
        index = len(data["ell"])//2
        contribution = (data["cocoa_shell_native"][index]**2
                        * data["cocoa_dchi_weights_Mpc"]
                        * data["cocoa_cap_variance_Mpc"])
        if np.any(contribution < 0) or contribution.sum() <= 0:
            raise ValueError("Expected nonnegative native SSC diagonal weights")
        fractions = contribution/contribution.sum()
        cumulative = np.cumsum(fractions[order])
        visible = (redshift[order] >= .01) & (redshift[order] <= 2)
        axes[1, 1].plot(redshift[order][visible], 100*cumulative[visible],
                        label=rf"{name}: $\ell={data['ell'][index]:g}$")
        notes.append({"quantity": f"{name} native CoCoA SSC diagonal",
                      "ell": float(data["ell"][index]),
                      "fraction_below_displayed_z": float(
                          fractions[redshift < .01].sum()),
                      "fraction_above_displayed_z": float(
                          fractions[redshift > 2].sum())})
    axes[1, 1].set(xlabel="Redshift z", ylabel="Cumulative SSC variance [%]",
                   title="Where the native CoCoA SSC signal accumulates")
    axes[1, 1].legend(fontsize=9)
    for ax in axes.flat:
        ax.set_xlim(.01, 2)
        ax.grid(alpha=.2)
    masks = sum(note.get("masked_samples", 0) for note in notes)
    fig.suptitle("SSC ingredients along the line of sight · equal physical "
                 r"$k=(\ell+1/2)/\chi_{\rm CoCoA}(z)$" "\n"
                 "Displayed range: z = 0.01–2; all original nodes/signs "
                 "are retained in the archive.\n"
                 f"Near-zero denominator masks: {masks}. "
                 "Dotted: CCL k extrapolation. Ingredient ratios are not "
                 "covariance errors.",
                 fontsize=11)
    save_figure(fig, output, "ssc_model_inputs")
    plt.close(fig)
    return notes


def formula_figure(cases, output):
    """Change the response prescription while retaining actual CCL inputs."""
    import numpy as np
    import matplotlib.pyplot as plt

    data = cases[0][2]
    for key in ("formula_a", "formula_k_Mpc_inverse",
                "formula_ccl_response_Mpc3", "formula_ccl_convention_Mpc3",
                "formula_two_halo_convention_Mpc3",
                "formula_fractional_nonlinear_convention_Mpc3"):
        np.testing.assert_array_equal(data[key], cases[1][2][key])
    redshifts = 1/data["formula_a"]-1
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
    notes = []
    for ax, target in zip(axes, (.1, .5, 1.0)):
        matches = np.flatnonzero(np.isclose(redshifts, target,
                                            rtol=0, atol=1e-12))
        if len(matches) != 1:
            raise ValueError(f"Missing exact saved formula redshift {target}")
        index = matches[0]
        for key, label, style in (
                ("ccl_convention", "Linear-power convention (CCL)", "--"),
                ("two_halo_convention", "Two-halo response", "-"),
                ("fractional_nonlinear_convention",
                 "Fractional halo response × nonlinear power", "-")):
            delta, note = ratio_percent(
                data[f"formula_{key}_Mpc3"][index],
                data["formula_ccl_response_Mpc3"][index], key)
            ax.semilogx(data["formula_k_Mpc_inverse"], 1+delta/100,
                        linestyle=style, label=label)
            note["redshift"] = target
            notes.append(note)
        ax.set(xlabel=r"Physical wavenumber k [Mpc$^{-1}$]",
               ylabel="Response / CCL linear-power response",
               title=f"z = {target:g}")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=7, loc="best")
    fig.suptitle("Common CCL halo moments and power spectra · "
                 "response prescriptions only\n"
                 "Actual CCL and CoCoA APIs; these are not native CoCoA "
                 "halo-model predictions.", fontsize=12)
    save_figure(fig, output, "ssc_response_conventions")
    plt.close(fig)
    return notes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("low", "high", "low-ccl", "high-ccl", "results", "figures"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    if args.results.exists() or args.figures.exists():
        parser.error("Choose fresh results and figures directories")

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    cases = [read_case(args.low, args.low_ccl),
             read_case(args.high, args.high_ccl)]
    if (cases[0][1]["input_manifest_sha256"]
            != cases[1][1]["input_manifest_sha256"]):
        raise ValueError("Low/high cases do not share the same survey inputs")
    args.results.mkdir(parents=True)
    args.figures.mkdir(parents=True)

    # Each final diagnostic already contains every parent CCL array.
    # Verify this and archive that complete union once per multipole set.
    # Preserve all four original manifests byte-for-byte alongside it.
    files = {}
    arrays = {}
    provenance = args.results / "provenance"
    provenance.mkdir()
    for name, source in (("low", args.low), ("high", args.high),
                         ("low_ccl", args.low_ccl),
                         ("high_ccl", args.high_ccl)):
        destination = provenance / f"{name}.json"
        shutil.copy2(source / "manifest.json", destination)
        files[f"provenance/{name}.json"] = sha256(destination)
    for name, case, parent in zip(("low", "high"), cases,
                                  (args.low_ccl, args.high_ccl)):
        with np.load(parent / "models.npz", allow_pickle=False) as saved:
            for key in saved.files:
                np.testing.assert_array_equal(saved[key], case[2][key])
        arrays.update({f"{name}__{key}": value
                       for key, value in case[2].items()})
    np.savez_compressed(args.results / "models.npz", **arrays)
    files["models.npz"] = sha256(args.results / "models.npz")
    matrix_figure(cases, args.figures)
    sightline = sightline_figure(cases, args.figures)
    formulas = formula_figure(cases, args.figures)
    write_json(args.results / "comparison.json", {
        "schema": "tjpcov-ssc-model-figures-v1", "status": "completed",
        "script_sha256": sha256(Path(__file__)), "files": files,
        "figures": {path.name: sha256(path)
                    for path in sorted(args.figures.iterdir())},
        "matrix_normalization": "100*(substitution-CoCoA)/sqrt(Cii*Cjj); "
                                "C is native CoCoA SSC, not total covariance",
        "estimator": "point multipoles; no band averaging",
        "cases": {name: record for name, (record, _, _) in
                  zip(("low", "high"), cases)},
        "sightline_masks_and_ranges": sightline,
        "formula_masks": formulas,
        "scope": "Saved correctness diagnostics; no performance measurement",
    })
    print(f"Archived complete arrays in {args.results}")
    print(f"Saved three PNG/PDF figure pairs in {args.figures}")


if __name__ == "__main__":
    main()
