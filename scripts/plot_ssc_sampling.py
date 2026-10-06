"""Plot saved native SSC sampling results; never reinterpret them as timings."""

import argparse
import json
from pathlib import Path

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    archive = args.report.parent / "ssc_native.npz"
    if sha256(archive) != report["files"][archive.name]:
        raise ValueError("Saved SSC arrays changed after validation")

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = np.load(archive, allow_pickle=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for name, label, style in (("ssc_low_a16", "785 response-a nodes", "-"),
                                ("ssc_low_a2", "99 response-a nodes", "s"),
                                ("ssc_low_qag", "50 response-a nodes", "o")):
        redshift = 1/data[name+"__window_a"]-1
        use = redshift <= 0.15
        # Plot values actually evaluated by CCL. The line only joins the
        # fine samples; it is not a reconstruction of CCL's interpolant.
        axes[0].semilogy(redshift[use], data[name+"__disc_variance"][use],
                          style, label=label, fillstyle="none")
    axes[0].set(xlabel="Redshift z", ylabel="CCL disc variance [Mpc]",
                title="Rapid change near the observer")
    axes[0].legend(fontsize=8)
    reference = np.diag(data["ssc_low_a16__ssc"])
    ell = data["ssc_low_a16__ell"]
    for name, label in (("ssc_low_qag", "50 nodes"),
                        ("ssc_low_a2", "99 nodes"),
                        ("ssc_low_a4", "197 nodes"),
                        ("ssc_low_a8", "393 nodes")):
        residual = 100*(np.diag(data[name+"__ssc"])/reference-1)
        axes[1].plot(ell, residual, "o-", label=label)
    axes[1].axhline(0, color="0.4", linewidth=0.6)
    axes[1].set(xlabel="Effective multipole", ylabel="SSC variance / fine − 1 [%]",
                title="Native TJPCov: time-grid refinement")
    axes[1].legend(fontsize=8)
    fig.suptitle("LSST Y1 source bin 3 · CCL background-variance sampling")
    args.output.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(args.output / f"ssc_sampling.{suffix}", dpi=180,
                    facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
