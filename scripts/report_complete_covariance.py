"""Verify and plot complete CoCoA/TJPCov source-bin Fourier covariances.

The TJPCov Gaussian/SSC base and chosen cNG run can come from different
native runs, provided every input fingerprint and coordinate agrees.
All sums use actual saved components. Matrices retain their signed entries
and native antisymmetry; symmetric parts are used only for variance modes.

This is the same four-panel scientific comparison as the OneCov report:
100 (CoCoA - TJPCov) / sqrt(TJPCov_total[ii] TJPCov_total[jj]). Gaussian
uses discrete ell-weighted bands; SSC and cNG use nominal band centres.
The plotted total therefore retains TJPCov's mixed estimator convention.
"""

import argparse
import os
from pathlib import Path

from common import load_bundle, sha256, write_json


COMPONENTS = ("gaussian", "ssc", "cng", "total")


def read_run(folder, code):
    """Verify an input manifest and return only the compact covariance archive."""
    import numpy as np

    record = load_bundle(folder, f"{code}-complete-shear-v1")
    if record["status"] not in ("completed", "nonpositive-total"):
        raise ValueError(f"The run did not finish its numerical calculation: {folder}")
    with np.load(folder / "covariance.npz", allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    for name in ("ell", "edges"):
        if name not in arrays or not np.isfinite(arrays[name]).all():
            raise ValueError(f"Missing or invalid {name}: {folder}")
    np.testing.assert_array_equal(arrays["ell"], np.arange(30., 3001., 30.))
    np.testing.assert_array_equal(arrays["edges"], np.arange(15., 3016., 30.))
    for name in COMPONENTS:
        if name in arrays and (arrays[name].shape != (100, 100)
                               or not np.isfinite(arrays[name]).all()):
            raise ValueError(f"Invalid {name} matrix: {folder}")
    return arrays, record


def full_components(arrays):
    """Sum the three actual native components without matrix repair."""
    result = {name: arrays[name] for name in COMPONENTS[:-1]}
    result["total"] = result["gaussian"]+result["ssc"]+result["cng"]
    return result


def rms_product(matrix):
    """Return a diagonal RMS product, or None for nonpositive variances."""
    import numpy as np

    diagonal = np.diag(matrix)
    if np.any(diagonal <= 0) or not np.isfinite(diagonal).all():
        return None
    rms = np.sqrt(diagonal)
    return rms[:, None]*rms[None, :]


def matrix_difference(candidate, reference, total_scale):
    """Retain whole-matrix, diagonal and signed normalized differences."""
    import numpy as np

    delta = candidate-reference
    scaled = delta/total_scale
    relative = np.diag(delta)[np.diag(reference) != 0]
    relative /= np.diag(reference)[np.diag(reference) != 0]
    norm = float(np.linalg.norm(reference))
    diagonal_reference = np.diag(reference)
    return {
        "max_total_rms_scaled_difference_percent": float(100*np.max(np.abs(scaled))),
        "rms_total_rms_scaled_difference_percent": float(100*np.sqrt(np.mean(scaled**2))),
        "signed_total_rms_scaled_range_percent": [float(100*scaled.min()),
                                                   float(100*scaled.max())],
        "max_relative_diagonal_difference_percent": float(100*np.max(abs(relative)))
            if len(relative) else None,
        "rms_relative_diagonal_difference_percent": float(100*np.sqrt(np.mean(relative**2)))
            if len(relative) else None,
        "signed_relative_diagonal_range_percent": [float(100*relative.min()),
                                                    float(100*relative.max())]
            if len(relative) else None,
        "zero_reference_diagonals": int(np.count_nonzero(diagonal_reference == 0)),
        "fractional_frobenius_difference_percent": float(100*np.linalg.norm(delta)/norm)
            if norm else None,
        "absolute_frobenius_difference": float(np.linalg.norm(delta)),
        "max_asymmetry_total_rms_scaled": {
            "candidate": float(np.max(np.abs(candidate-candidate.T)/total_scale)),
            "reference": float(np.max(np.abs(reference-reference.T)/total_scale)),
        },
    }


def variance_modes(candidate, reference, scale):
    """Check native total symmetric parts and generalized variance ratios."""
    import numpy as np
    from scipy.linalg import cholesky, eigvalsh

    symmetric = {"candidate": (candidate+candidate.T)/(2*scale),
                 "reference": (reference+reference.T)/(2*scale)}
    result = {"convention": "Symmetric parts only for eigenvalues; saved matrices unchanged"}
    positive = {}
    for name, matrix in symmetric.items():
        values = eigvalsh(matrix)
        result[f"{name}_scaled_eigenvalue_range"] = [float(values[0]), float(values[-1])]
        try:
            cholesky(matrix, lower=True)
            positive[name] = bool(values[0] > 0)
        except np.linalg.LinAlgError:
            positive[name] = False
    result["positive_definite"] = positive
    result["both_positive_definite"] = all(positive.values())
    result["generalized_variance_ratio_range"] = None
    result["max_generalized_variance_change_percent"] = None
    if positive["reference"]:
        ratios = eigvalsh(symmetric["candidate"], symmetric["reference"])
        result["generalized_variance_ratio_range"] = [float(ratios[0]), float(ratios[-1])]
        result["max_generalized_variance_change_percent"] = float(100*np.max(abs(ratios-1)))
    return result


def refinement_difference(coarse, fine, common_gaussian, common_ssc):
    """Compare every cNG entry and its resulting full total to the finer run."""
    import numpy as np

    coarse_total = common_gaussian+common_ssc+coarse
    fine_total = common_gaussian+common_ssc+fine
    scale = rms_product(fine_total)
    if scale is None:
        raise ValueError("A refinement total has nonpositive diagonal variances")
    result = {"difference": "coarse minus fine",
              "normalization": "finer total diagonal RMS product",
              "cng": matrix_difference(coarse, fine, scale),
              "total": matrix_difference(coarse_total, fine_total, scale),
              "variance_modes": variance_modes(coarse_total, fine_total, scale)}
    component_scale = rms_product(fine)
    result["cng_diagonal_rms_normalization_available"] = component_scale is not None
    if component_scale is not None:
        normalized = (coarse-fine)/component_scale
        result["cng_rms_scaled_difference_percent"] = {
            "maximum": float(100*np.max(np.abs(normalized))),
            "rms": float(100*np.sqrt(np.mean(normalized**2))),
        }
    else:
        result["cng_rms_scaled_difference_percent"] = None
        result["cng_normalization_note"] = "Nonpositive cNG diagonal; no absolute-value or clipped RMS used"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tjpcov-base", type=Path, required=True)
    parser.add_argument("--tjpcov-cng", type=Path, required=True)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--cocoa-refined", type=Path)
    parser.add_argument("--cng-refinement", type=Path, action="append", default=[],
                        help="Repeat in coarse-to-fine order; chosen cNG is appended last if absent")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.figures.exists():
        parser.error("Choose new result and figure directories")

    import numpy as np
    from scipy.linalg import eigvalsh
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        if os.environ.get(name) != "1":
            raise ValueError(f"Set {name}=1 before importing numerical libraries")

    folders = {"tjpcov_base": args.tjpcov_base, "tjpcov_cng": args.tjpcov_cng,
               "cocoa": args.cocoa}
    if args.cocoa_refined is not None:
        folders["cocoa_refined"] = args.cocoa_refined
    refinement_folders = []
    for path in args.cng_refinement:
        if path.resolve() not in [p.resolve() for p in refinement_folders]:
            refinement_folders.append(path)
    if refinement_folders and args.tjpcov_cng.resolve() not in [
            p.resolve() for p in refinement_folders]:
        refinement_folders.append(args.tjpcov_cng)
    for index, path in enumerate(refinement_folders):
        folders[f"cng_refinement_{index}"] = path

    arrays, manifests = {}, {}
    for label, path in folders.items():
        code = "cocoa" if label.startswith("cocoa") else "tjpcov"
        arrays[label], manifests[label] = read_run(path, code)
    input_hash = manifests["cocoa"]["input_manifest_sha256"]
    for label, record in manifests.items():
        if record["input_manifest_sha256"] != input_hash:
            raise ValueError(f"Different cosmology/survey/power bundle: {label}")
        for key in ("ell", "edges"):
            np.testing.assert_array_equal(arrays[label][key], arrays["cocoa"][key],
                                          err_msg=f"Coordinates differ: {label}/{key}")
    if manifests["cocoa"]["source_bin_1based"] != 3:
        raise ValueError("Expected LSST source bin 3")
    for label, record in manifests.items():
        if not label.startswith("cocoa"):
            if (record["source_name"] != "source3"
                    or record["input_manifest"]["source_bin_1based"] != 3):
                raise ValueError(f"A native run has a different source selection: {label}")

    cocoa = full_components(arrays["cocoa"])
    tjpcov = full_components({"gaussian": arrays["tjpcov_base"]["gaussian"],
                             "ssc": arrays["tjpcov_base"]["ssc"],
                             "cng": arrays["tjpcov_cng"]["cng"]})
    # CoCoA's saved mean is band-averaged. The separately saved centre
    # signal, not that mean, matches the TJPCov runner's signal samples.
    co_signal = arrays["cocoa"]["signal_centres"]
    tj_signal = arrays["tjpcov_base"]["signal"]
    if (co_signal.shape != (100,) or tj_signal.shape != (100,)
            or not np.isfinite(co_signal).all() or not np.isfinite(tj_signal).all()):
        raise ValueError("Invalid native point-sampled shear signals")
    scale = rms_product(tjpcov["total"])
    if scale is None:
        raise ValueError("TJPCov total has nonpositive variances; cannot define the requested residual")

    report = {
        "schema": "cocoa-tjpcov-complete-shear-comparison-v1",
        "scope": "LSST Y1 source bin 3; complete 100x100 Fourier pilot",
        "status": "checked", "ndata": 100, "masked_entries": 0,
        "difference": "CoCoA minus TJPCov",
        "normalization": "TJPCov total diagonal RMS product",
        "estimator": "Gaussian: discrete ell-weighted bands, upper edge excluded; SSC/cNG: nominal centres",
        "input_manifest_sha256": input_hash,
        "native_component_sources": {"gaussian": "tjpcov_base", "ssc": "tjpcov_base",
                                     "cng": "tjpcov_cng"},
        "provenance": {label: {"manifest_sha256": sha256(path / "manifest.json"),
                               "manifest": manifests[label]}
                       for label, path in folders.items()},
        "script_sha256": sha256(__file__), "components": {},
        "cng_refinements": [], "cocoa_refinement": None,
    }
    nonzero = tj_signal != 0
    relative = (co_signal[nonzero]-tj_signal[nonzero])/tj_signal[nonzero]
    report["signal_centres"] = {
        "difference": "CoCoA minus TJPCov",
        "maximum_fractional_difference_percent": float(100*np.max(abs(relative)))
            if len(relative) else None,
        "rms_fractional_difference_percent": float(100*np.sqrt(np.mean(relative**2)))
            if len(relative) else None,
        "zero_reference_entries": int(np.count_nonzero(~nonzero)),
    }
    with np.errstate(invalid="raise", divide="raise"):
        report["variance_modes"] = variance_modes(cocoa["total"], tjpcov["total"], scale)
        reference = (tjpcov["total"]+tjpcov["total"].T)/(2*scale)
        for name in COMPONENTS:
            item = matrix_difference(cocoa[name], tjpcov[name], scale)
            item["total_variance_mode_change_range_percent"] = None
            if report["variance_modes"]["positive_definite"]["reference"]:
                delta = cocoa[name]-tjpcov[name]
                modes = eigvalsh((delta+delta.T)/(2*scale), reference)
                item["total_variance_mode_change_range_percent"] = [float(100*modes[0]),
                                                                    float(100*modes[-1])]
            report["components"][name] = item
        for index in range(len(refinement_folders)-1):
            left, right = f"cng_refinement_{index}", f"cng_refinement_{index+1}"
            entry = refinement_difference(arrays[left]["cng"], arrays[right]["cng"],
                                           tjpcov["gaussian"], tjpcov["ssc"])
            entry.update(coarse=left, fine=right,
                         coarse_settings=manifests[left]["components"]["cng"],
                         fine_settings=manifests[right]["components"]["cng"])
            report["cng_refinements"].append(entry)
        if "cocoa_refined" in arrays:
            fine = full_components(arrays["cocoa_refined"])
            fine_scale = rms_product(fine["total"])
            if fine_scale is None:
                raise ValueError("Refined CoCoA total has nonpositive diagonal variances")
            report["cocoa_refinement"] = {
                "difference": "baseline minus refined", "normalization": "refined CoCoA total diagonal RMS product",
                "components": {name: matrix_difference(cocoa[name], fine[name], fine_scale)
                               for name in COMPONENTS},
                "variance_modes": variance_modes(cocoa["total"], fine["total"], fine_scale),
                "baseline_settings": manifests["cocoa"]["settings"],
                "refined_settings": manifests["cocoa_refined"]["settings"],
            }

    report["both_positive_definite"] = report["variance_modes"]["both_positive_definite"]
    report["status"] = "completed" if report["both_positive_definite"] else "nonpositive-total"
    args.output.mkdir(parents=True)
    args.figures.mkdir(parents=True)
    compact = {"ell": arrays["cocoa"]["ell"], "edges": arrays["cocoa"]["edges"],
               "cocoa_signal_centres": co_signal, "tjpcov_signal_centres": tj_signal,
               "cocoa_signal_band_average": arrays["cocoa"]["signal"]}
    for label, data in (("cocoa", cocoa), ("tjpcov", tjpcov)):
        compact.update({f"{label}_{name}": data[name] for name in COMPONENTS})
    for index in range(len(refinement_folders)):
        label = f"cng_refinement_{index}"
        compact[f"{label}_cng"] = arrays[label]["cng"]
    if "cocoa_refined" in arrays:
        compact.update({f"cocoa_refined_{name}": value for name, value in
                        full_components(arrays["cocoa_refined"]).items()})
    np.savez_compressed(args.output / "matrices.npz", **compact)
    report["files"] = {"matrices.npz": sha256(args.output / "matrices.npz")}
    # Save the diagnostic record before plotting, including any failed
    # positivity check. There is no tolerance reset or hidden repair.
    write_json(args.output / "report.json", report)
    make_figures(cocoa, tjpcov, scale, compact["ell"], args.figures, report)
    report["figures"] = {path.name: sha256(path) for path in sorted(args.figures.iterdir())}
    write_json(args.output / "report.json", report)
    print(f"Saved complete matrix diagnostics and figures in {args.output} and {args.figures}")
    if not report["both_positive_definite"]:
        raise SystemExit("A total symmetric part is not positive definite; raw matrices and report preserved")


def make_figures(cocoa, tjpcov, scale, ell, folder, report):
    """Use the OneCov four-panel layout plus native structure and variances."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
                         "axes.labelsize": 22, "axes.titlesize": 23,
                         "xtick.labelsize": 17, "ytick.labelsize": 17,
                         "axes.linewidth": 1.0})
    ticks = np.array([0, 24, 49, 74, 99])
    labels = ell[ticks].astype(int)

    def matrix_axis(axis):
        axis.set(xlabel=r"$\ell$", ylabel=r"$\ell'$", xticks=ticks, yticks=ticks,
                 xticklabels=labels, yticklabels=labels)

    def save(figure, name):
        for suffix in ("png", "pdf"):
            figure.savefig(folder / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
        plt.close(figure)

    def residual_panel(figure, axis, values, title):
        limit = float(np.max(abs(values)))
        if limit == 0:
            axis.imshow(np.zeros_like(values), cmap="RdBu_r", vmin=-1, vmax=1,
                        origin="upper", interpolation="nearest", rasterized=True)
            axis.text(.5, .5, "Exactly zero difference", transform=axis.transAxes,
                      ha="center", va="center", fontsize=18)
        else:
            image = axis.imshow(values, origin="upper", cmap="RdBu_r", vmin=-limit,
                                vmax=limit, interpolation="nearest", rasterized=True)
            bar = figure.colorbar(image, ax=axis, shrink=.88, pad=.03)
            bar.set_label("Difference / TJPCov total rms product [%]", size=16)
            bar.ax.tick_params(labelsize=14)
        axis.set_title(title)
        matrix_axis(axis)

    figure, axes = plt.subplots(2, 2, figsize=(14, 12), layout="constrained")
    for axis, name, title in zip(axes.ravel(), COMPONENTS,
                                 ("Gaussian", "SSC", "Connected non-Gaussian", "Total")):
        residual_panel(figure, axis, 100*(cocoa[name]-tjpcov[name])/scale, title)
    warning = "" if report["both_positive_definite"] else "\nWARNING: a total is not positive definite"
    figure.suptitle("LSST Y1 source bin 3 · CoCoA − TJPCov\n"
                   "100 × 100 Fourier covariance · native halo prescriptions\n"
                   "Every entry retained; each panel has its own colour scale"+warning,
                   fontsize=23)
    save(figure, "complete_shear_difference")

    figure, axes = plt.subplots(1, 3, figsize=(18, 5.8), layout="constrained")
    component_scales = [rms_product(data["cng"]) for data in (cocoa, tjpcov)]
    correlations = all(item is not None for item in component_scales)
    if not correlations:
        component_scales = [rms_product(data["total"]) for data in (cocoa, tjpcov)]
        if any(item is None for item in component_scales):
            component_scales = [scale, scale]
    normalized = [data["cng"]/denominator for data, denominator in
                  zip((cocoa, tjpcov), component_scales)]
    limit = max(1.0, *(float(np.max(abs(values))) for values in normalized))
    for axis, name, values in zip(axes[:2], ("CoCoA", "TJPCov"), normalized):
        image = axis.imshow(values, origin="upper", cmap="RdBu_r", vmin=-limit,
                            vmax=limit, interpolation="nearest", rasterized=True)
        axis.set_title(name)
        matrix_axis(axis)
    label = "cNG / own cNG diagonal rms product" if correlations else "cNG / total diagonal rms product"
    bar = figure.colorbar(image, ax=list(axes[:2]), shrink=.86)
    bar.set_label("cNG / own cNG rms product" if correlations
                  else "cNG / total rms product", size=16)
    residual_panel(figure, axes[2], 100*(cocoa["cng"]-tjpcov["cng"])/scale,
                   "CoCoA − TJPCov")
    figure.suptitle("LSST Y1 source bin 3 · connected covariance structure\n"
                   "Native signed matrices; no eigenvalue clipping", fontsize=22)
    save(figure, "connected_covariance_structure")
    report["cng_structure_figure_normalization"] = label

    figure, axes = plt.subplots(2, 3, figsize=(17, 10), layout="constrained")
    colors = {"cocoa": "#2563a6", "tjpcov": "#c05621"}
    for column, (name, title) in enumerate(zip(COMPONENTS[:-1],
                                               ("Gaussian", "SSC", "Connected non-Gaussian"))):
        for code, data in (("cocoa", cocoa), ("tjpcov", tjpcov)):
            diagonal = np.diag(data[name])
            axes[0, column].plot(ell, diagonal, color=colors[code],
                                 ls="-" if code == "cocoa" else "--",
                                 label="CoCoA" if code == "cocoa" else "TJPCov")
            axes[1, column].plot(ell, 100*diagonal/np.diag(data["total"]),
                                 color=colors[code], ls="-" if code == "cocoa" else "--")
        both = np.concatenate([np.diag(data[name]) for data in (cocoa, tjpcov)])
        if np.all(both > 0):
            axes[0, column].set_yscale("log")
        elif np.any(both != 0):
            axes[0, column].set_yscale("symlog", linthresh=np.max(abs(both))*1.e-4)
        axes[0, column].set(title=title, ylabel="Diagonal covariance")
        axes[0, column].legend(fontsize=17)
        axes[1, column].set(ylabel="Fraction of own total variance [%]")
        for row in (0, 1):
            axes[row, column].set(xlabel=r"Multipole $\ell$", xscale="log")
            axes[row, column].grid(alpha=.2)
    figure.suptitle("LSST Y1 source bin 3 · native covariance components\n"
                   "Gaussian: discrete band averages · SSC/cNG: band centres", fontsize=23)
    save(figure, "complete_shear_components")


if __name__ == "__main__":
    main()
