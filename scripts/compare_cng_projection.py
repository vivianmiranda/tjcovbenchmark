"""Separate cNG projection from halo physics, using the actual native CCL Tk3D.

The export stage samples that table and CCL's native tracer/background on
CoCoA's saved radial nodes. The project stage uses CoCoA's C contraction
and precomputed GSL rule. This is an eight-multipole shared-input diagnostic,
not a native CoCoA halo-model prediction or an independent reference code.
"""
import argparse
import signal
import sys
from pathlib import Path

from common import load_bundle, require_thread_environment, sha256, write_json


def export(args):
    import numpy as np
    import pyccl as ccl
    from run_complete_native import cosmology, catalog_and_config
    native = load_bundle(args.native, "tjpcov-complete-shear-v1")
    cocoa = load_bundle(args.cocoa_run, "cocoa-complete-shear-v1")
    for saved in (native, cocoa):
        if saved["input_manifest_sha256"] != sha256(args.inputs/"manifest.json"):
            raise ValueError("Different input bundles")
    cfg = native["components"]["cng"]
    if args.epsrel is not None:
        ccl.gsl_params.INTEGRATION_LIMBER_EPSREL = args.epsrel
    cosmo, m, d = cosmology(args.inputs, cfg["a_refinement"], cfg["N_K"])
    _, source, setup, _, _, _ = catalog_and_config(cosmo,m,d,args.output/"tracer")
    # The Gaussian builder gives the same shear tracer without requiring
    # irrelevant HOD values; it does not compute a Gaussian covariance here.
    from tjpcov.covariance_gaussian_fsky import FourierGaussianFsky
    tracer = FourierGaussianFsky(setup).get_tracer_info()[0][source]
    tab = np.load(args.native/"native_trispectrum.npz", allow_pickle=False)
    tk = ccl.Tk3D(a_arr=tab["a_arr"], lk_arr=tab["lk_arr"],
                   tkk_arr=tab["tkk"], is_logt=False)
    raw = np.load(args.native/"covariance.npz", allow_pickle=False)
    co = np.load(args.cocoa_run/"projection_inputs.npz", allow_pickle=False)
    grid_settings = cocoa["settings"]
    grid_manifest = args.cocoa_run/"manifest.json"
    if args.radial_ssc is not None:
        radial = load_bundle(args.radial_ssc, "cocoa-ssc-shear-v1")
        if radial["input_manifest_sha256"] != native["input_manifest_sha256"]:
            raise ValueError("The saved refined radial grid used different inputs")
        if radial["settings"]["a_edges"] != grid_settings["a_edges"]:
            raise ValueError("Refined radial panel support differs")
        co = np.load(args.radial_ssc/"ssc.npz", allow_pickle=False)
        grid_settings = radial["settings"]
        grid_manifest = args.radial_ssc/"manifest.json"
    take = np.linspace(0,len(raw["ell"])-1,8,dtype=int)
    ell = raw["ell"][take]
    ref = ccl.angular_cl_cov_cNG(cosmo, tracer, tracer, ell=ell, t_of_kk_a=tk,
        tracer3=tracer, tracer4=tracer, fsky=setup["tjpcov"]["fsky"],
        integration_method=cfg["integration_method"])
    expected = raw["cng"][np.ix_(take,take)]
    scale = np.sqrt(abs(np.outer(np.diag(expected),np.diag(expected))))
    gate = float(np.max(abs(ref-expected)/scale))
    if args.epsrel is None and gate > 1e-8:
        raise ValueError(f"Saved native Tk3D replay failed: {gate}")
    a = co["geometry"][0]
    chi = cosmo.comoving_radial_distance(a)
    kernel = tracer.get_kernel(chi[::-1])[0,::-1]
    spin = tracer.get_f_ell(ell)[0]/(ell+.5)**2
    np.testing.assert_array_equal(tracer.get_bessel_derivative(),[-1])
    np.testing.assert_array_equal(tracer.get_angles_derivative(),[2])
    power = np.asarray([tk((ell+.5)/distance,float(aa)) for aa,distance in zip(a,chi)])
    if not np.isfinite(power).all():
        raise ValueError("Non-finite supplied trispectrum")
    np.savez_compressed(args.output/"projection_inputs.npz", a=a, chi=chi,
        hubble=cosmo.h_over_h0(a), kernel=kernel, spin=spin, ell=ell,
        trispectrum=power, reference=ref,
        area_sr=np.asarray(m["area_deg2"]*(np.pi/180)**2))
    write_json(args.output/"manifest.json", {
        "schema":"cng-shared-projection-inputs-v1", "status":"completed",
        "native_manifest_sha256":sha256(args.native/"manifest.json"),
        "cocoa_manifest_sha256":sha256(args.cocoa_run/"manifest.json"),
        "input_manifest_sha256":sha256(args.inputs/"manifest.json"),
        "replay_max_rms_residual":gate, "projection_epsrel_control":args.epsrel,
        "h":m["cosmology"]["H0"]/100,
        "cocoa_settings":grid_settings,
        "radial_grid_manifest_sha256":sha256(grid_manifest),
        "units":{"distance":"Mpc","kernel":"Mpc^-1","trispectrum":"Mpc^9"},
        "axes":{"trispectrum":["a","ell2","ell1"]},
        "ccl_binary_sha256":sha256(ccl._ccllib.__file__),
        "radial_endpoints_a":list(grid_settings["a_edges"])[::len(grid_settings["a_edges"])-1],
        "native_tracer_chi_bounds_Mpc":[float(tracer._trc[0].chi_min),float(tracer._trc[0].chi_max)],
        "gauss_panel_chi_bounds_Mpc":cosmo.comoving_radial_distance(np.asarray(grid_settings["a_edges"])[[0,-1]]).tolist(),
        "script_sha256":sha256(__file__),
        "files":{"projection_inputs.npz":sha256(args.output/"projection_inputs.npz")}})


def project(args):
    import numpy as np
    import cosmolike_lsst_y1_interface as ci
    saved = load_bundle(args.export,"cng-shared-projection-inputs-v1")
    d = np.load(args.export/"projection_inputs.npz",allow_pickle=False)
    rule = ci.covariance.covariance_integration_rule(nquad=saved["cocoa_settings"]["radial_nquad"])
    edges = saved["cocoa_settings"]["a_edges"]
    a = np.concatenate([(lo+hi)/2+(hi-lo)*rule[0]/2 for lo,hi in zip(edges[:-1],edges[1:])])
    wa = np.concatenate([(hi-lo)*rule[1]/2 for lo,hi in zip(edges[:-1],edges[1:])])
    np.testing.assert_allclose(a,d["a"],rtol=5e-14,atol=0)
    dchi = wa * (2997.92458/saved["h"]) / (a*a*d["hubble"])
    measure = dchi/(d["area_sr"]*d["chi"]**6)
    # CCL axis order is [a,ell2,ell1]; keep every native signed entry.
    t = d["trispectrum"].transpose(1,2,0)
    t = t*d["spin"][:,None,None]**2*d["spin"][None,:,None]**2
    def contract(mask):
        return ci.covariance.covariance_project(
            left=np.ascontiguousarray(t.reshape(64,-1)),
            right=np.ascontiguousarray(d["kernel"][None,:]**4),
            weight=np.ascontiguousarray(measure*mask)).reshape(8,8)
    matrix = contract(np.ones(len(a)))
    ref = d["reference"]
    rms = np.sqrt(abs(np.outer(np.diag(ref),np.diag(ref))))
    scaled = (matrix-ref)/rms
    np.savez_compressed(args.output/"comparison.npz",cocoa=matrix,tjpcov=ref,
                        ell=d["ell"],scaled_difference=scaled)
    write_json(args.output/"manifest.json",{
        "schema":"cng-shared-projection-v1","status":"completed",
        "scope":"Same CCL trispectrum, tracer, spin and background; distinct radial integrators",
        "export_manifest_sha256":sha256(args.export/"manifest.json"),
        "max_abs_rms_residual":float(np.max(abs(scaled))),
        "max_fractional_residual":float(np.max(abs(matrix/ref-1))),
        "interface_sha256":sha256(ci.__file__), "script_sha256":sha256(__file__),
        "files":{"comparison.npz":sha256(args.output/"comparison.npz")}})
    print(f"Shared cNG projection maximum RMS-normalized residual: {np.max(abs(scaled)):.8g}",flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="mode",required=True)
    ex=sub.add_parser("export")
    ex.add_argument("--epsrel",type=float,default=None,
                    help="Explicit CCL integration tolerance control, holding Tk3D fixed")
    ex.add_argument("--radial-ssc",type=Path,default=None,
                    help="Reuse the actual radial nodes of a completed refined CoCoA SSC run")
    for name in ("inputs","native","cocoa_run"):
        ex.add_argument(name,type=Path)
    pr=sub.add_parser("project")
    pr.add_argument("export",type=Path)
    for q in (ex,pr):
        q.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    require_thread_environment()
    if args.output.exists():
        p.error("Choose a fresh output directory")
    args.output.mkdir(parents=True)
    signal.signal(signal.SIGALRM,signal.SIG_DFL)
    signal.alarm(600)
    (export if args.mode=="export" else project)(args)
    signal.alarm(0)


if __name__=="__main__":
    main()
