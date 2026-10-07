"""Save CoCoA's five native matter trispectra and their actual ingredients.

The nine wavenumbers cover 0.001--10 h/Mpc at three redshifts. All 45
unordered pairs are retained, including equal pairs. Power and trispectra
are converted from core units to (Mpc/h)^3 and (Mpc/h)^9 respectively.
The angular nodes are also exported for a shared-power CCL diagnostic.
"""

import argparse
import signal
import sys
from pathlib import Path

from common import (load_bundle, require_thread_environment, revision,
                    sha256, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--integration-accuracy", type=int, choices=range(3),
                        default=0)
    args = parser.parse_args()
    threads = require_thread_environment()
    if threads > 6 or args.output.exists():
        parser.error("Use at most six threads and a new output directory")
    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(600)

    import numpy as np

    cocoa = args.cocoa.resolve()
    project = cocoa / "projects/lsst_y1"
    core = cocoa / "external_modules/code/cosmolike_core"
    sys.path[:0] = [str(project), str(project / "covariance"), str(core)]
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.geometry import angular_rule

    settings = configuration(gaussian={"nonlimber": False, "ia": "none"},
                             integration_accuracy=args.integration_accuracy)
    if settings["cosmology"] != inputs["cosmology"]:
        raise ValueError("Cosmology differs from the input bundle")
    tables = initialize(interface=ci, settings=settings)
    original = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    np.testing.assert_array_equal(tables["z_2D"], original["z"])
    np.testing.assert_array_equal(10**tables["log10k_2D"], original["k_h_mpc"])
    for name in ("linear", "nonlinear", "linear_cb"):
        actual = np.exp(tables[f"lnP_{name}"].reshape(
            original[f"p_{name}"].shape, order="F"))
        np.testing.assert_array_equal(actual, original[f"p_{name}"])

    k = np.geomspace(0.001, 10, 9)
    redshift = np.array([0.1, 0.5, 1.0])
    first, second = np.triu_indices(len(k))
    pairs = np.array([k[first], k[second]])
    _, weight, corner = angular_rule(
        nquad=settings["tree_nquad"], npanel=settings["tree_npanel"],
        interface=ci.covariance)

    # The stable |K+Q| expression avoids cancellation near opposite,
    # equal vectors. CCL will read its power at these exact same nodes.
    internal_k = np.sqrt((pairs[0]-pairs[1])[:, None]**2
                         +2*pairs[0, :, None]*pairs[1, :, None]*corner)
    length = 2997.92458
    rows = []
    for z in redshift:
        a = float(1/(1+z))
        i11, moments = ci.covariance.covariance_halo_moments(
            a=np.array([a]), k=np.ascontiguousarray(k[None, :]*length),
            lnm_edges=settings["lnm_edges"], nquad=settings["halo_mass_nquad"])
        power = ci.covariance.covariance_power(
            a=a, k=k*length, linear=True)*length**3
        internal = ci.covariance.covariance_power(
            a=a, k=internal_k*length, linear=True)*length**3
        # Each extra mass profile contributes one power of volume.
        # Moment roles are I02,I12,I13(K,Q,Q),I13(Q,K,K),I04.
        moments = moments[:, 0, :]*np.array(
            [length**3, length**3, length**6, length**6, length**9])[:, None]
        pk = power[[first, second]]
        single = i11[0, [first, second]]
        tree = ci.covariance.covariance_tree_averages(
            k=pairs, pk=pk, corner=corner, weight=weight, ps=internal)
        terms = ci.covariance.covariance_halo_trispectrum(
            pk=pk, i11=single, moments=np.ascontiguousarray(moments), tree=tree)
        rows.append(dict(terms=terms, moments=moments, i11=single,
                         power=power, internal_power=internal, tree=tree))
        print(f"CoCoA five halo terms: z={z}", flush=True)
    arrays = {name: np.array([row[name] for row in rows]) for name in rows[0]}
    if not all(np.isfinite(value).all() for value in arrays.values()):
        raise ValueError("Non-finite CoCoA trispectrum ingredient")
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output / "trispectrum.npz", k=k,
                        redshift=redshift, first=first, second=second,
                        corner=corner, weight=weight, internal_k=internal_k,
                        **arrays)
    write_json(args.output / "manifest.json", {
        "schema": "cocoa-trispectrum-v1", "status": "completed",
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "halo_order": ["1h", "2h13", "2h22", "3h", "4h"],
        "units": {"k": "h/Mpc", "power": "(Mpc/h)^3", "terms": "(Mpc/h)^9"},
        "axes": {"terms": ["redshift", "halo_order", "unordered_pair"]},
        "integration_accuracy": args.integration_accuracy,
        "halo_mass_nquad": settings["halo_mass_nquad"],
        "tree_nquad": settings["tree_nquad"], "tree_npanel": settings["tree_npanel"],
        "lnm_edges": settings["lnm_edges"].tolist(),
        "core": revision(core), "lsst_y1": revision(project),
        "interface_sha256": sha256(ci.__file__), "threads": threads,
        "script_sha256": sha256(__file__),
        "power_table_check": "linear and nonlinear inputs bitwise equal",
        "timing_scope": "accuracy only",
        "files": {"trispectrum.npz": sha256(args.output / "trispectrum.npz")},
    })


if __name__ == "__main__":
    main()
