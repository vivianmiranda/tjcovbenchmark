"""Stage verified refresh results under the README's existing public paths.

This reads completed calculations; it never runs CoCoA, CCL or TJPCov.
All 63 refresh stages must have finished successfully. Reports and arrays
are checked against their recorded hashes before copying. Existing tracked
results, figures and README text are never overwritten by this script.
The output is a fresh staging tree for visual review and later publication.
"""

import argparse
import json
import shutil
from pathlib import Path

from common import sha256, write_json
from refresh_comparison import SCRIPTS, plan, verify_record


def completed_stages(campaign):
    """Check the last recorded attempt of every stage, including plot stages."""
    latest = {}
    journals = sorted(campaign.glob("refresh_*.json"))
    for filename in journals:
        record = json.loads(filename.read_text())
        if record["schema"] != "tjpcov-production-power-refresh-v1":
            raise ValueError(f"Unexpected refresh journal: {filename}")
        for stage in record["stages"]:
            if "status" in stage:
                latest[stage["name"]] = stage
    expected = plan(campaign, Path("unused-cocoa"), Path("unused-tjpcov"))
    for item in expected:
        name = item["name"]
        stage = latest.get(name)
        if stage is None or stage["status"] != "completed":
            raise ValueError(f"No completed final attempt for stage {name}")
        if (stage["script"] != item["script"]
                or stage["script_sha256"] != sha256(SCRIPTS / item["script"])):
            raise ValueError(f"Stage script changed since execution: {name}")
        if item["record"] is not None:
            report = Path(item["record"])
            if (stage["record"] != str(report)
                    or verify_record(report) != stage["record_sha256"]):
                raise ValueError(f"Stage result changed since execution: {name}")
    return journals, {item["name"]: latest[item["name"]] for item in expected}


def halo_refinements(campaign):
    """Summarize saved ingredient ratios, without evaluating a halo model."""
    import numpy as np

    checks = []
    for coarse, fine in (("halo_native_m128", "halo_native_m255"),
                         ("halo_cocoa_i0", "halo_cocoa_i1"),
                         ("halo_cocoa_i1", "halo_cocoa_i2")):
        with np.load(campaign / coarse / "halo.npz") as left, np.load(
                campaign / fine / "halo.npz") as right:
            for axis in ("redshift", "mass", "k"):
                np.testing.assert_array_equal(left[axis], right[axis])
            delta = {}
            for key in ("sigma", "dndlnm", "bias", "concentration",
                        "linear", "i11", "i02", "i12"):
                if (not np.all(np.isfinite(left[key]))
                        or not np.all(np.isfinite(right[key]))
                        or np.any(right[key] <= 0)):
                    raise ValueError(f"Invalid positive halo ingredient: {key}")
                delta[key] = float(np.max(np.abs(left[key]/right[key]-1)))
        checks.append({
            "coarse": coarse, "fine": fine,
            "normalization": "largest absolute coarse/fine - 1",
            "max_fractional_change": delta,
        })
    return {"scope": "mass quadrature; not sigma or full-matrix convergence",
            "checks": checks}


def trispectrum_locations(campaign):
    """Report the actual largest 4h discrepancy and its weight in the sum."""
    import numpy as np

    with np.load(campaign / "results/trispectrum/trispectra.npz") as data:
        cocoa = data["cocoa_i0__terms"]
        ccl = data["native_m255__terms"]
        first, second = data["cocoa_i0__first"], data["cocoa_i0__second"]
        k, redshift = data["cocoa_i0__k"], data["cocoa_i0__redshift"]
        rows = []
        for iz, z in enumerate(redshift):
            if np.any(ccl[iz, 4] == 0):
                raise ValueError("A zero 4h denominator needs a zero-safe metric")
            relative = cocoa[iz, 4]/ccl[iz, 4]-1
            pair = int(np.argmax(np.abs(relative)))
            total = ccl[iz, :, pair].sum()
            if total == 0:
                raise ValueError("Cannot report a 4h fraction of a zero sum")
            rows.append({
                "redshift": float(z), "pair_index": pair,
                "K_h_Mpc": float(k[first[pair]]),
                "Q_h_Mpc": float(k[second[pair]]),
                "signed_cocoa_over_ccl_minus_one": float(relative[pair]),
                "ccl_4h_fraction_of_sum_at_pair": float(ccl[iz, 4, pair]/total),
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output", type=Path, required=True,
                        help="new staging root; never the repository root")
    args = parser.parse_args()
    campaign, output = args.campaign.resolve(), args.output.resolve()
    if output.exists():
        parser.error("Use a fresh staging tree; current public results stay intact")
    journals, stages = completed_stages(campaign)
    inputs = json.loads((campaign / "lsst_y1/manifest.json").read_text())
    preparation = inputs["power_preparation"]
    if (preparation["accuracy_boost"] != 1
            or preparation["refinement"] != 8
            or preparation["native_k_nodes"] != 1500
            or preparation["installed_k_nodes"] != 11993
            or not preparation["original_nodes_retained_exactly"]
            or not preparation["installed_tables_reproduced_exactly"]):
        raise ValueError("This publication is for the verified 11993-node baseline")
    input_hash = sha256(campaign / "lsst_y1/manifest.json")
    copied = {}

    def copy_file(source, target):
        """Copy bytes, check their identity and record the publication path."""
        destination = output / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        digest = sha256(source)
        if sha256(destination) != digest:
            raise RuntimeError(f"Publication copy differs: {target}")
        copied[str(target)] = digest

    def copy_report(source, target, name):
        report_path = source / name
        verify_record(report_path)
        report = json.loads(report_path.read_text())
        if ("input_manifest_sha256" in report
                and report["input_manifest_sha256"] != input_hash):
            raise ValueError(f"Different power bundle in {report_path}")
        if ("interface_sha256" in report
                and report["interface_sha256"] != inputs["interface_sha256"]):
            raise ValueError(f"CoCoA binary changed in {report_path}")
        copy_file(report_path, target / name)
        for filename in report.get("files", {}):
            copy_file(source / filename, target / filename)
        return report

    def figures(source, target, names, expected=None):
        for name in names:
            for suffix in ("png", "pdf"):
                filename = f"{name}.{suffix}"
                if not (source / filename).is_file():
                    raise FileNotFoundError(source / filename)
                if (expected is not None
                        and sha256(source / filename) != expected[filename]):
                    raise ValueError(f"Figure changed after plotting: {filename}")
                copy_file(source / filename, target / filename)

    results, images = Path("results"), Path("figures")
    summary = {"schema": "tjpcov-publication-summary-v1",
               "input_manifest_sha256": input_hash, "inputs": inputs,
               "completed_stages": len(stages), "gaussian": {}, "ssc": {},
               "timing_policy": "Correctness refresh; no timing comparison",
               "scope": "Gaussian and SSC small matrices; halo ingredients "
                        "and separated trispectra. Projected cNG/full totals "
                        "and real-space comparisons remain pending."}
    for kind in ("shear", "3x2"):
        for band in ("low", "high"):
            native = campaign / f"gaussian_{kind}_{band}"
            comparison = campaign / f"comparison_{kind}_{band}"
            public = f"gaussian_{band}" if kind == "shear" else native.name
            copy_report(native, results / public, "manifest.json")
            summary["gaussian"][public] = copy_report(
                comparison, results / public, "comparison.json")
            figures(campaign / "figures" / native.name, images / public,
                    ("gaussian_matrices", "gaussian_components"))

    # Keep the existing flat SSC links. Its archive retains all 16 native
    # matrices. Save the separate four-way sampling diagnostic alongside it.
    ssc = campaign / "results/ssc"
    copy_file(ssc / "comparison.json", results / "ssc_native.json")
    copy_file(ssc / "ssc_native.npz", results / "ssc_native.npz")
    summary["ssc"]["refinements"] = json.loads(
        (ssc / "comparison.json").read_text())["refinements"]
    summary["ssc"]["sampling"] = copy_report(
        campaign / "ssc_sampling_a2", results / "ssc_sampling", "report.json")
    figures(campaign / "figures/ssc", images / "ssc_sampling", ("ssc_sampling",))
    for band in ("low", "high"):
        for label in ("base", "fine"):
            public = f"ssc_{band}_refined" if label == "fine" else f"ssc_{band}_base"
            source = campaign / f"ssc_comparison_{band}_{label}"
            summary["ssc"][public] = copy_report(
                source, results / public, "comparison.json")
            figures(campaign / "figures" / source.name, images / public,
                    ("ssc_matrices",))
    summary["ssc"]["models"] = copy_report(
        campaign / "results/ssc_models", results / "ssc_models", "comparison.json")
    figures(campaign / "figures/ssc_models", images / "ssc_models",
            ("ssc_model_swaps", "ssc_model_inputs", "ssc_response_conventions"),
            expected=summary["ssc"]["models"]["figures"])

    summary["halo"] = copy_report(campaign / "halo_comparison_i2",
                                  results / "halo", "comparison.json")
    for name in ("halo_native_m128", "halo_native_m255",
                 "halo_cocoa_i0", "halo_cocoa_i1", "halo_cocoa_i2"):
        copy_report(campaign / name, results / "halo" / name, "manifest.json")
    figures(campaign / "figures/halo_comparison_i2", images / "halo",
            ("halo_ingredients",))
    refinements = halo_refinements(campaign)
    write_json(output / "results/halo/refinements.json", refinements)
    copied["results/halo/refinements.json"] = sha256(
        output / "results/halo/refinements.json")
    summary["halo_refinements"] = refinements

    summary["trispectrum"] = copy_report(
        campaign / "results/trispectrum", results / "trispectrum", "comparison.json")
    # Preserve the old report URL; its own file map stays valid by copying
    # the diagnostic array beside it.
    copy_file(campaign / "trispectrum_power_refinement/report.json",
              results / "trispectrum/power_refinement.json")
    copy_file(campaign / "trispectrum_power_refinement/diagnostic.npz",
              results / "trispectrum/diagnostic.npz")
    summary["trispectrum_4h_maximum_locations"] = trispectrum_locations(campaign)
    figures(campaign / "figures/trispectrum", images / "trispectrum",
            ("trispectrum_terms", "trispectrum_pairs",
             "trispectrum_power_diagnostic", "trispectrum_model_diagnostic"),
            expected=summary["trispectrum"]["figures"])
    summary["galaxy_bias"] = copy_report(campaign / "galaxy_bias_placement",
                                         results / "galaxy_bias", "report.json")

    # Native CAMB and installed arrays remain in the fresh campaign bundle;
    # their verified hashes and full configuration are in this public record.
    # This avoids duplicating dense input arrays in every result directory.
    write_json(output / "results/refresh_summary.json", summary)
    copied["results/refresh_summary.json"] = sha256(
        output / "results/refresh_summary.json")
    for path in journals:
        copy_file(path, results / "refresh_provenance" / path.name)
    write_json(output / "publication.json", {
        "schema": "tjpcov-publication-staging-v1", "status": "prepared",
        "script_sha256": sha256(__file__), "files": copied,
        "pending": "Visual figure review, measured README rewrite and explicit "
                   "copy into tracked folders; this script publishes nothing.",
    })
    print(f"Verified staging tree: {output}")
    print("README numbers: results/refresh_summary.json; review before copying.")


if __name__ == "__main__":
    main()
