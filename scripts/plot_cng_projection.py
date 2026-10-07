"""Publish the eight-multipole shared-input CoCoA/CCL cNG projection.

Both matrix panels use the same CCL cNG diagonal RMS normalization. The
third retains every signed CoCoA-minus-CCL residual in percent. This is a
shared-Tk3D projection diagnostic, not a native halo-model comparison.
The finite radial endpoints differ, so its residual is not attributed
solely to the quadrature algorithm.
"""

import argparse
import shutil
from pathlib import Path

from common import load_bundle, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--export", type=Path,
                        help="Verified shared-input export to package beside the comparison")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.figures.exists():
        parser.error("Choose fresh result and figure directories")
    source = load_bundle(args.comparison, "cng-shared-projection-v1")
    if source["status"] != "completed":
        raise ValueError("The shared-input projection did not complete")
    exported = None
    if args.export is not None:
        exported = load_bundle(args.export, "cng-shared-projection-inputs-v1")
        if (exported["status"] != "completed"
                or sha256(args.export / "manifest.json")
                != source["export_manifest_sha256"]):
            raise ValueError("The supplied export does not belong to this completed comparison")

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with np.load(args.comparison / "comparison.npz", allow_pickle=False) as archive:
        cocoa = archive["cocoa"]
        ccl = archive["tjpcov"]
        ell = archive["ell"]
        saved_difference = archive["scaled_difference"]
    if (cocoa.shape != (8, 8) or ccl.shape != (8, 8) or ell.shape != (8,)
            or not all(np.isfinite(values).all() for values in (cocoa, ccl, ell))
            or np.any(np.diff(ell) <= 0)):
        raise ValueError("Malformed eight-multipole projection archive")
    diagonal = np.diag(ccl)
    if np.any(diagonal <= 0):
        raise ValueError("CCL cNG diagonals must be positive for the requested RMS normalization")
    scale = np.sqrt(diagonal[:, None]*diagonal[None, :])
    difference = (cocoa-ccl)/scale
    np.testing.assert_allclose(saved_difference, difference, rtol=5.e-14, atol=1.e-17)
    np.testing.assert_allclose(np.max(abs(difference)), source["max_abs_rms_residual"],
                               rtol=5.e-14, atol=1.e-17)
    if exported is not None:
        with np.load(args.export / "projection_inputs.npz", allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive["ell"], ell)
            np.testing.assert_array_equal(archive["reference"], ccl)
    args.output.mkdir(parents=True)
    args.figures.mkdir(parents=True)
    # Copy the verified compact scientific archive exactly, preserving
    # all native entries and any small antisymmetric component.
    shutil.copy2(args.comparison / "comparison.npz", args.output / "comparison.npz")
    shutil.copy2(args.comparison / "manifest.json", args.output / "source_manifest.json")
    copied = ["comparison.npz", "source_manifest.json"]
    if exported is not None:
        # This is the compact eight-mode evaluation, not the much larger
        # native Tk3D grid. Its parent native-manifest hash stays recorded.
        shutil.copy2(args.export / "manifest.json", args.output / "inputs_manifest.json")
        shutil.copy2(args.export / "projection_inputs.npz", args.output / "projection_inputs.npz")
        if sha256(args.output / "projection_inputs.npz") != exported["files"]["projection_inputs.npz"]:
            raise ValueError("Copied projection inputs do not match their original fingerprint")
        copied.extend(("inputs_manifest.json", "projection_inputs.npz"))
    report = {
        "schema": "cng-shared-projection-publication-v1", "status": "completed",
        "scope": "Eight multipoles; shared CCL Tk3D, background, shear kernel and spin factors",
        "comparison": "CoCoA fixed line-of-sight quadrature versus CCL native projection",
        "normalization": "CCL cNG diagonal RMS product, common to both matrix panels",
        "residual": "100 times (CoCoA minus CCL) divided by CCL cNG diagonal RMS product",
        "limits": "Common physical inputs with finite-endpoint and resolution differences; not a formal identical-domain exactness test",
        "shape": [8, 8], "masked_entries": 0, "ell": ell.tolist(),
        "maximum_rms_scaled_difference_percent": float(100*np.max(abs(difference))),
        "rms_scaled_difference_percent": float(100*np.sqrt(np.mean(difference**2))),
        "maximum_rms_scaled_asymmetry": {
            "cocoa": float(np.max(abs(cocoa-cocoa.T)/scale)),
            "ccl": float(np.max(abs(ccl-ccl.T)/scale)),
        },
        "source_manifest_sha256": sha256(args.comparison / "manifest.json"),
        "export_manifest_sha256": source["export_manifest_sha256"],
        "source": source, "export": exported,
        "export_packaging": None if exported is None else {
            "manifest_copy": "inputs_manifest.json",
            "original_manifest_name": "manifest.json",
            "scientific_file_names_and_hashes": {
                "projection_inputs.npz": exported["files"]["projection_inputs.npz"]},
            "native_manifest_sha256": exported["native_manifest_sha256"],
            "native_trispectrum_copied": False,
            "note": "The compact sampled inputs are included; the full native Tk3D is identified through its parent manifest",
        },
        "script_sha256": sha256(__file__),
        "files": {name: sha256(args.output / name)
                  for name in copied},
    }
    plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
                         "axes.labelsize": 19, "axes.titlesize": 20,
                         "xtick.labelsize": 14, "ytick.labelsize": 14})
    figure, axes = plt.subplots(1, 3, figsize=(18, 6.4), layout="constrained")
    normalized = (cocoa/scale, ccl/scale)
    bound = max(float(np.max(abs(values))) for values in normalized)
    for axis, title, values in zip(axes[:2],
                                  ("CoCoA fixed quadrature", "CCL native projection"),
                                  normalized):
        image = axis.imshow(values, origin="upper", cmap="RdBu_r", vmin=-bound,
                            vmax=bound, interpolation="nearest", rasterized=True)
        axis.set_title(title)
    figure.colorbar(image, ax=list(axes[:2]), shrink=.82, pad=.03,
                    label="Covariance / CCL cNG rms product")
    residual = 100*difference
    maximum = float(np.max(abs(residual)))
    image = axes[2].imshow(residual, origin="upper", cmap="RdBu_r",
                           vmin=-maximum if maximum else -1,
                           vmax=maximum if maximum else 1,
                           interpolation="nearest", rasterized=True)
    axes[2].set_title("CoCoA − CCL")
    if maximum:
        figure.colorbar(image, ax=axes[2], shrink=.82, pad=.03,
                        label="Difference / CCL cNG rms product [%]")
    else:
        axes[2].text(.5, .5, "Exactly zero difference", ha="center", va="center",
                     transform=axes[2].transAxes)
    for axis in axes:
        axis.set(xlabel=r"$\ell_1$", ylabel=r"$\ell_2$",
                 xticks=np.arange(8), yticks=np.arange(8),
                 xticklabels=[f"{value:g}" for value in ell],
                 yticklabels=[f"{value:g}" for value in ell])
        axis.tick_params(axis="x", labelrotation=45)
    figure.suptitle("LSST Y1 source bin 3 · shared CCL trispectrum, background and shear kernel\n"
                   "Eight multipoles · every matrix entry retained", fontsize=22)
    figure.supxlabel("Radial endpoint differences are retained; the residual is not purely quadrature error",
                       fontsize=16)
    figure_files = []
    for suffix in ("png", "pdf"):
        path = args.figures / f"connected_projection.{suffix}"
        figure.savefig(path, dpi=180, bbox_inches="tight")
        figure_files.append(path)
    plt.close(figure)
    report["figures"] = {path.name: sha256(path) for path in figure_files}
    write_json(args.output / "manifest.json", report)
    print(f"Published shared cNG projection in {args.output} and {args.figures}")


if __name__ == "__main__":
    main()
