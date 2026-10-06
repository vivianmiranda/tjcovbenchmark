"""Compare native halo orders and isolate assembly with supplied CCL inputs.

Run in the CoCoA environment. The shared-input result still uses CoCoA's
actual angular kernels and assembly. CCL supplies every halo moment and
every linear power sample; the 1h value is copied, so its equality is a
consistency check rather than an independent mass-integration test.
"""

import argparse
from pathlib import Path

from common import load_bundle, require_thread_environment, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cocoa_run", type=Path)
    parser.add_argument("native_run", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_thread_environment()
    if args.output.exists():
        parser.error("Choose a fresh output directory")
    cocoa = load_bundle(args.cocoa_run, "cocoa-trispectrum-v1")
    native = load_bundle(args.native_run, "tjpcov-trispectrum-v1")
    if cocoa["input_manifest_sha256"] != native["input_manifest_sha256"]:
        raise ValueError("Different cosmological inputs")
    if native["cocoa_manifest_sha256"] != sha256(args.cocoa_run / "manifest.json"):
        raise ValueError("CCL internal power uses another CoCoA quadrature")

    import numpy as np
    import cosmolike_lsst_y1_interface as ci

    if sha256(ci.__file__) != cocoa["interface_sha256"]:
        raise ValueError("CoCoA library differs from the native export")
    co = np.load(args.cocoa_run / "trispectrum.npz", allow_pickle=False)
    cc = np.load(args.native_run / "trispectrum.npz", allow_pickle=False)
    for name in ("k", "redshift", "first", "second"):
        np.testing.assert_array_equal(co[name], cc[name])
    first, second = co["first"], co["second"]
    pairs = np.array([co["k"][first], co["k"][second]])
    shared = []
    shared_power = []
    for index, z in enumerate(co["redshift"]):
        pk = cc["power"][index][[first, second]]
        tree = ci.covariance.covariance_tree_averages(
            k=pairs, pk=pk, corner=co["corner"], weight=co["weight"],
            ps=np.ascontiguousarray(cc["internal_power"][index]))
        for table, source in ((shared, cc), (shared_power, co)):
            table.append(ci.covariance.covariance_halo_trispectrum(
                pk=pk, tree=tree, i11=np.ascontiguousarray(source["i11"][index]),
                moments=np.ascontiguousarray(source["moments"][index])))
    shared, shared_power = np.array(shared), np.array(shared_power)
    if not np.isfinite(shared).all() or not np.isfinite(shared_power).all():
        raise ValueError("Non-finite shared-input result")

    def metrics(actual, reference):
        rows = []
        for index, z in enumerate(co["redshift"]):
            delta = actual[index]-reference[index]
            if np.any(reference[index] == 0):
                raise ValueError("Use a zero-safe metric for a zero term")
            relative = delta/reference[index]
            rows.append({
                "redshift": float(z),
                "fractional_range": np.stack(
                    [relative.min(axis=1), relative.max(axis=1)], axis=1).tolist(),
                "maximum_absolute_fractional": np.max(
                    np.abs(relative), axis=1).tolist(),
                "maximum_diagonal_fractional": np.max(
                    np.abs(relative[:, first == second]), axis=1).tolist(),
                "maximum_offdiagonal_fractional": np.max(
                    np.abs(relative[:, first != second]), axis=1).tolist(),
                "sum_fractional_range": [float(x) for x in (
                    np.min(actual[index].sum(axis=0)/reference[index].sum(axis=0)-1),
                    np.max(actual[index].sum(axis=0)/reference[index].sum(axis=0)-1))],
            })
        return rows

    args.output.mkdir(parents=True)
    np.savez_compressed(args.output / "comparison.npz", k=co["k"],
                        redshift=co["redshift"], first=first, second=second,
                        cocoa=co["terms"], ccl=cc["terms"],
                        shared=shared, shared_power=shared_power)
    report = {
        "schema": "trispectrum-comparison-v1", "status": "completed",
        "cocoa": cocoa, "ccl": native, "halo_order": cocoa["halo_order"],
        "native": metrics(co["terms"], cc["terms"]),
        "common_power_and_moments": metrics(shared, cc["terms"]),
        "power_reader_effect": metrics(shared_power, co["terms"]),
        "cocoa_manifest_sha256": sha256(args.cocoa_run / "manifest.json"),
        "native_manifest_sha256": sha256(args.native_run / "manifest.json"),
        "script_sha256": sha256(__file__),
        "files": {"comparison.npz": sha256(args.output / "comparison.npz")},
        "scope": "Matter trispectra, not projected covariance; copied 1h consistency check",
    }
    write_json(args.output / "manifest.json", report)
    print({"native": report["native"],
           "common_power_and_moments": report["common_power_and_moments"]})


if __name__ == "__main__":
    main()
