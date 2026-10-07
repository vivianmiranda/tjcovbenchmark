"""Verify complete-pilot adapters with the actual production C kernels.

Shared TJPCov spectra isolate Gaussian assembly. Replaying eight selected
CoCoA cNG modes through its dedicated connected projector checks the
memory-saving one-source contraction used by run_cocoa_complete.py.
"""
import argparse
from pathlib import Path
from common import load_bundle, require_thread_environment, sha256, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cocoa",type=Path)
    p.add_argument("native",type=Path)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():
        p.error("Choose a fresh output file")
    require_thread_environment()
    import numpy as np
    import cosmolike_lsst_y1_interface as ci
    cm=load_bundle(args.cocoa,"cocoa-complete-shear-v1")
    tm=load_bundle(args.native,"tjpcov-complete-shear-v1")
    if cm["input_manifest_sha256"]!=tm["input_manifest_sha256"]:
        raise ValueError("Input bundles differ")
    co=np.load(args.cocoa/"covariance.npz",allow_pickle=False)
    tj=np.load(args.native/"covariance.npz",allow_pickle=False)
    inp=np.load(args.cocoa/"projection_inputs.npz",allow_pickle=False)
    np.testing.assert_array_equal(co["gaussian_ell"],tj["spectra_ell"])
    gauss=ci.covariance.covariance_gaussian_fourier(
        spectra=np.ascontiguousarray(tj["spectra"][:,None,None]),
        noise=np.atleast_1d(tj["noise"]),pairs=np.array([[0,0]],dtype=np.int32),
        operators=np.ascontiguousarray(co["operators"]),ell_min=15,area_sr=cm["area_sr"])
    scale=np.sqrt(np.outer(tj["gaussian"].diagonal(),tj["gaussian"].diagonal()))
    gerr=float(np.max(abs(gauss-tj["gaussian"])/scale))
    ell=inp["ell"]
    spin=np.sqrt((ell-1)*ell*(ell+1)*(ell+2))/(ell+.5)**2
    n=len(ell)
    t=np.zeros((4*n,4*n,inp["geometry"].shape[1]))
    t[:n,:n]=inp["matter_trispectrum"]*spin[:,None,None]**2*spin[None,:,None]**2
    val=ci.covariance.covariance_project_connected(
        probes=np.array([0],dtype=np.int32),projected=t,
        pair_window=np.ascontiguousarray(inp["window"][None,:]**2),
        measure=np.ascontiguousarray(inp["geometry"][3]/(cm["area_sr"]*inp["geometry"][2]**6)))
    index=cm["projection_archive"]["indices"]
    expected=co["cng"][np.ix_(index,index)]
    scale=np.sqrt(abs(np.outer(expected.diagonal(),expected.diagonal())))
    cerr=float(np.max(abs(val-expected)/scale))
    passed=gerr<1e-12 and cerr<1e-12
    write_json(args.output,{"passed":passed,"gaussian_shared_spectra_rms_residual":gerr,
        "cng_connected_projector_replay_rms_residual":cerr,
        "cocoa_manifest_sha256":sha256(args.cocoa/"manifest.json"),
        "native_manifest_sha256":sha256(args.native/"manifest.json"),
        "interface_sha256":sha256(ci.__file__),"script_sha256":sha256(__file__)})
    print(f"Shared Gaussian residual {gerr:.5g}; dedicated cNG projector replay {cerr:.5g}")
    if not passed:
        raise ValueError("Production adapter verification failed")


if __name__=="__main__":
    main()
