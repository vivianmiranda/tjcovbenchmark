"""Trace native SSC differences with public CCL and supplied-input CoCoA APIs.

Run the ``ccl`` stage in the TJPCov environment, then ``cocoa`` in the
CoCoA environment. No library is patched and no covariance is independently
implemented. Both stages require fresh outputs and have a 600 s deadline.
The saved native runs select low or high effective multipoles.

The main 2x2 test holds CoCoA geometry/tracer weights fixed while swapping
native response and window predictions. Those predictions include their
own halo, power-reader and background choices. A separate projection uses
CCL ingredients throughout. Small formula-only response tests retain CCL
halo ingredients; they do not constitute native CoCoA predictions.
"""

import argparse
import inspect
import signal
import sys
from pathlib import Path

from common import load_bundle, require_thread_environment, sha256, write_json


def matrix_difference(matrix, reference):
    """Report every-entry differences in units of reference diagonal rms."""
    import numpy as np

    diagonal = np.diag(reference)
    if np.any(diagonal <= 0):
        raise ValueError("Need positive reference variances")
    rms = np.sqrt(diagonal)
    return {
        "max_abs_difference_over_reference_rms_product": float(np.max(
            np.abs((matrix-reference)/rms[:, None]/rms[None, :]))),
        "diagonal_fractional_difference": (np.diag(matrix)/diagonal-1).tolist(),
    }


def save_stage(output, record, arrays):
    """Keep partial scientific arrays and provenance if a later gate fails."""
    import numpy as np

    if any(not np.all(np.isfinite(value)) for value in arrays.values()):
        raise ValueError("Non-finite exported diagnostic array")
    np.savez_compressed(output / "models.npz", **arrays)
    record["files"] = {"models.npz": sha256(output / "models.npz")}
    write_json(output / "manifest.json", record)


def export_ccl(args, record):
    """Recreate native CCL SSC and export its ingredients on CoCoA nodes."""
    import numpy as np
    import pyccl as ccl
    import sacc
    import tjpcov
    from tjpcov.covariance_fourier_ssc_fsky import FourierSSCHaloModelFsky

    inputs = load_bundle(args.inputs, "lsst-y1-tjpcov-inputs-v1")
    gaussian = load_bundle(args.gaussian_run, "tjpcov-gaussian-shear-v1")
    native = load_bundle(args.native_ssc, "tjpcov-ssc-shear-v1")
    cocoa = load_bundle(args.cocoa_ssc, "cocoa-ssc-shear-v1")
    for run in (native, cocoa):
        if (run["status"] != "completed"
                or run["input_manifest_sha256"]
                != sha256(args.inputs / "manifest.json")
                or run["gaussian_manifest_sha256"]
                != sha256(args.gaussian_run / "manifest.json")):
            raise ValueError("Native SSC runs must share these completed inputs")
    if (inputs["cosmology"]["mnu"] != 0
            or inputs["gaussian"]["ia"] != "none"
            or inputs["gaussian"]["nonlimber"]
            or inputs["rsd"] or inputs["magnification"]):
        raise ValueError("Need massless, Limber, zero-IA/RSD/magnification shear")

    package = Path(tjpcov.__file__).resolve().parent
    for name, expected in native["installed_source"].items():
        if (sha256(package / name) != expected
                or sha256(args.tjpcov / "tjpcov" / name) != expected):
            raise ValueError(f"Changed TJPCov source: {name}")
    if (sha256(ccl._ccllib.__file__) != native["ccl_binary_sha256"]
            or sha256(inspect.getfile(ccl.halos.halomod_Tk3D_SSC_linear_bias))
            != native["ccl_ssc_source_sha256"]):
        raise ValueError("CCL changed since the saved native SSC run")

    # Match the refined native run, including its CAMB-support adaptation.
    for name in ("A_SPLINE_MINLOG_PK", "A_SPLINE_MIN_PK",
                 "A_SPLINE_NA_PK", "A_SPLINE_NLOG_PK"):
        setattr(ccl.spline_params, name, native["a_domain_control"][name])
    for key, name in (("N_K", "N_K"), ("K_MIN_Mpc_inverse", "K_MIN"),
                      ("K_MAX_Mpc_inverse", "K_MAX")):
        setattr(ccl.spline_params, name, native["k_sampling"][key])
    raw = np.load(args.inputs / "inputs.npz", allow_pickle=False)
    cosmos = inputs["cosmology"]
    h = cosmos["H0"]/100
    length = 2997.92458/h  # one CoCoA length unit, in physical Mpc
    if (inputs["units"]["k"] != "h/Mpc"
            or inputs["units"]["power"] != "(Mpc/h)^3"
            or inputs["units"]["power_axes"] != ["redshift", "wavenumber"]):
        raise ValueError("Unexpected CAMB units or axes")
    powers = {
        name: {"a": 1/(1+raw["z"][::-1]), "k": raw["k_h_mpc"]*h,
               "delta_matter:delta_matter": raw[f"p_{name}"][::-1]/h**3}
        for name in ("linear", "nonlinear")}
    cosmo = ccl.CosmologyCalculator(
        Omega_c=cosmos["omegam"]-cosmos["omegab"],
        Omega_b=cosmos["omegab"], h=h, n_s=cosmos["ns"],
        A_s=cosmos["As_1e9"]*1e-9, m_nu=0,
        w0=cosmos["w"], wa=cosmos["w0pwa"]-cosmos["w"],
        pk_linear=powers["linear"], pk_nonlin=powers["nonlinear"])
    original = np.load(args.native_ssc / "ssc.npz", allow_pickle=False)
    saved_cocoa = np.load(args.cocoa_ssc / "ssc.npz", allow_pickle=False)
    np.testing.assert_array_equal(cosmo.get_pk_spline_a(),
                                  original["native_response_a"])
    np.testing.assert_array_equal(cosmo.get_pk_spline_lk(),
                                  original["native_response_lk"])

    catalog = sacc.Sacc.load_fits(str(args.gaussian_run / "source_bins.fits"))
    source = gaussian["source_name"]
    if (len(catalog.tracers) != 1
            or catalog.tracers[source].quantity != "galaxy_shear"
            or "lens" in source):
        raise ValueError("This diagnostic requires one shear tracer")
    fsky = gaussian["area_sr"]/(4*np.pi)
    calculator = FourierSSCHaloModelFsky({"tjpcov": {
        "cosmo": cosmo, "sacc_file": catalog, "IA": None, "use_mpi": False,
        "outdir": str(args.output / "tracer_setup"), "fsky": fsky,
        f"Ngal_{source}": inputs["source_density_arcmin2"],
        f"sigma_e_{source}": inputs["sigma_e_component"]}})
    tracers, _ = calculator.get_tracer_info()
    tracer = tracers[source]
    ell = calculator.get_ell_eff()
    for data in (original, saved_cocoa):
        np.testing.assert_array_equal(ell, data["ell"])
    if ell.shape != (5,):
        raise ValueError("Expected the existing five-multipole pilot")

    mass_def = ccl.halos.MassDef200m
    hmc = ccl.halos.HMCalculator(
        mass_function=ccl.halos.MassFuncTinker08(mass_def=mass_def),
        halo_bias=ccl.halos.HaloBiasTinker10(mass_def=mass_def),
        mass_def=mass_def, **native["halo_model"]["hmc_defaults"])
    profile = ccl.halos.HaloProfileNFW(
        mass_def=mass_def, fourier_analytic=True,
        concentration=ccl.halos.ConcentrationDuffy08(mass_def=mass_def))
    response = ccl.halos.halomod_Tk3D_SSC_linear_bias(
        cosmo=cosmo, hmc=hmc, prof=profile,
        a_arr=original["native_response_a"],
        lk_arr=original["native_response_lk"])

    def project(window_a, variance):
        return ccl.covariances.angular_cl_cov_SSC(
            cosmo, tracer1=tracer, tracer2=tracer, tracer3=tracer,
            tracer4=tracer, ell=ell, t_of_kk_a=response,
            sigma2_B=(window_a, variance),
            integration_method=native["integration_method"])

    variance = ccl.sigma2_B_disc(cosmo, a_arr=original["window_a"], fsky=fsky)
    baseline = project(original["window_a"], variance)
    record.update({
        "input_manifest_sha256": sha256(args.inputs / "manifest.json"),
        "gaussian_manifest_sha256": sha256(args.gaussian_run / "manifest.json"),
        "native_manifest_sha256": sha256(args.native_ssc / "manifest.json"),
        "cocoa_manifest_sha256": sha256(args.cocoa_ssc / "manifest.json"),
        "native_model": native, "cocoa_model": cocoa,
        "length_unit_Mpc": length,
        "baseline_gate": matrix_difference(baseline, original["ssc"]),
        "units": {"response": "Mpc^3", "window_variance": "Mpc",
                  "distance": "Mpc", "source_kernel": "Mpc^-1",
                  "covariance": "dimensionless"},
    })
    arrays = {"ell": ell, "native_ccl_ssc": original["ssc"],
              "ccl_baseline": baseline, "cocoa_native_ssc": saved_cocoa["ssc"]}
    save_stage(args.output, record, arrays)
    if record["baseline_gate"][
            "max_abs_difference_over_reference_rms_product"] > 1e-8:
        raise ValueError("Native CCL reproduction gate failed; inspect export")
    print("Native CCL reproduction passed", flush=True)

    # Extract the actual response factors instead of reconstructing their
    # physics or taking a square root of an outer product. Pk2D supplies
    # CCL's interpolation. Every requested row is checked against Tk3D.
    response_a, lk1, lk2, factors = response.get_spline_arrays()
    if len(factors) != 2:
        raise ValueError("Expected the native separable SSC response")
    np.testing.assert_array_equal(lk1, lk2)
    np.testing.assert_array_equal(factors[0], factors[1])
    factor = ccl.Pk2D(a_arr=response_a, lk_arr=lk1, pk_arr=factors[0],
                     is_logp=False, extrap_order_lok=response.extrap_order_lok,
                     extrap_order_hik=response.extrap_order_hik)
    a = saved_cocoa["geometry"][0]
    cocoa_chi = saved_cocoa["geometry"][2]*length
    ccl_chi = cosmo.comoving_radial_distance(a)
    if (np.any(np.diff(a) <= 0) or np.any(cocoa_chi <= 0)
            or np.any(ccl_chi <= 0) or a[0] < response_a[0]
            or a[-1] > response_a[-1]):
        raise ValueError("Unsupported saved radial grid")

    def sample_response(distance):
        wave = (ell[None, :]+0.5)/distance[:, None]
        values = []
        worst = 0.0
        for aa, kk in zip(a, wave):
            value = factor(kk, float(aa))
            expected = response(kk, float(aa))
            error = np.max(np.abs(np.outer(value, value)-expected))
            worst = max(worst, float(error/np.max(np.abs(expected))))
            values.append(value)
        if worst > 1e-10:
            raise ValueError("Pk2D factor does not reproduce native Tk3D")
        return wave, np.asarray(values), worst

    k_cocoa, at_cocoa_chi, gate1 = sample_response(cocoa_chi)
    k_ccl, at_ccl_chi, gate2 = sample_response(ccl_chi)
    np.testing.assert_array_equal(tracer.get_bessel_derivative(), [-1])
    np.testing.assert_array_equal(tracer.get_angles_derivative(), [2])
    # The public CCL shear tracer has j_l(k chi)/(k chi)^2. The Limber
    # replacement k chi=l+1/2 gives this spin factor for one shear leg.
    # CCL's ccl_cl_tracer_t_get_f_ell uses the exact factorial expression
    # only through ell=10, an asymptotic expression through ell=1000,
    # and (ell+1/2)^2 above it. Keep the actual public-API factor in the
    # shared projection; its difference from the exact expression is a
    # recorded model convention, not a failed equality check.
    spin = tracer.get_f_ell(ell)[0]/(ell+0.5)**2
    exact_spin = np.sqrt((ell-1)*ell*(ell+1)*(ell+2))/(ell+0.5)**2
    record["native_spin_versus_exact"] = {
        "max_absolute_difference": float(np.max(np.abs(spin-exact_spin))),
        "max_fractional_difference": float(np.max(
            np.abs(spin/exact_spin-1))),
        "source": "CCL src/ccl_tracers.c:ccl_cl_tracer_t_get_f_ell",
        "projection_uses": "Public tracer.get_f_ell, not exact polynomial",
    }
    np.testing.assert_array_equal(
        tracer.get_transfer(np.log(k_ccl[0]), float(a[0])), np.ones((1, 5)))
    # get_kernel expects increasing physical distance; the saved a grid
    # has decreasing distance, so reverse it only for this API call.
    kernel = tracer.get_kernel(ccl_chi[::-1])[0, ::-1]

    # Refine only the variance input to the CCL projection at the same
    # Gauss nodes used by CoCoA. Preserve the native endpoints, including
    # a=1. This exposes residual window interpolation in the projection
    # comparison rather than mistaking it for a halo-model difference.
    dense_a = np.unique(np.concatenate((original["window_a"], a)))
    dense_variance = ccl.sigma2_B_disc(cosmo, a_arr=dense_a, fsky=fsky)
    disc = dense_variance[np.searchsorted(dense_a, a)]
    dense_matrix = project(dense_a, dense_variance)
    arrays.update({
        "a": a, "cocoa_chi_Mpc": cocoa_chi, "ccl_chi_Mpc": ccl_chi,
        "ccl_E": cosmo.h_over_h0(a), "ccl_kernel_Mpc_inverse": kernel,
        "ccl_spin": spin, "exact_spin": exact_spin,
        "k_cocoa_Mpc_inverse": k_cocoa,
        "k_ccl_Mpc_inverse": k_ccl,
        "ccl_response_at_cocoa_k_Mpc3": at_cocoa_chi,
        "ccl_response_at_ccl_k_Mpc3": at_ccl_chi,
        "ccl_disc_at_cocoa_a_Mpc": disc,
        "ccl_dense_window_a": dense_a,
        "ccl_dense_window_variance_Mpc": dense_variance,
        "ccl_dense_window_ssc": dense_matrix,
        "native_response_a": response_a, "native_response_lk": lk1,
        "native_response_Mpc3": factors[0],
        "cocoa_response_Mpc3": saved_cocoa["response"]*length**3,
        "cocoa_cap_variance_Mpc": saved_cocoa["variance"]*length,
        "cocoa_kernel_Mpc_inverse": saved_cocoa["window"]/length,
    })
    record["response_factor_gates"] = [gate1, gate2]
    record["dense_window_versus_native"] = matrix_difference(
        dense_matrix, baseline)
    save_stage(args.output, record, arrays)
    print("Exported response, window and geometry along the line of sight",
          flush=True)

    # A separate small grid isolates the response convention with common
    # CCL halo fits. These are inputs to CoCoA's public algebra kernel,
    # not replacements of either code's native response above.
    fa = np.array([0.5, 2/3, 1/1.1, 1.0])
    fk = np.geomspace(0.001, 10.0, 33)*h
    ft = ccl.halos.halomod_Tk3D_SSC_linear_bias(
        cosmo=cosmo, hmc=hmc, prof=profile, a_arr=fa, lk_arr=np.log(fk))
    _, _, _, direct = ft.get_spline_arrays()
    pk = cosmo.get_linear_power()
    pnl = cosmo.get_nonlin_power()
    step = cocoa["settings"]["response_step"]
    fields = []
    pairs = ccl.halos.Profile2pt()
    for aa in fa:
        norm = profile.get_normalization(cosmo, float(aa), hmc=hmc)
        single = hmc.I_1_1(cosmo, fk, float(aa), profile)/norm
        i02 = hmc.I_0_2(cosmo, fk, float(aa), profile, prof_2pt=pairs)/norm**2
        i12 = hmc.I_1_2(cosmo, fk, float(aa), profile, prof_2pt=pairs)/norm**2
        low = hmc.I_1_1(cosmo, fk*np.exp(-step), float(aa), profile)/norm
        high = hmc.I_1_1(cosmo, fk*np.exp(step), float(aa), profile)/norm
        # Follow CoCoA's centered derivative of I11^2 P_linear while
        # retaining the same public CCL moments and power on both sides.
        slope = np.log(high**2*pk(fk*np.exp(step), float(aa))
                       /(low**2*pk(fk*np.exp(-step), float(aa))))/(2*step)
        fields.append([pk(fk, float(aa)), pnl(fk, float(aa)), single,
                       i02, i12, slope, pk(fk, float(aa), derivative=True)])
    arrays.update({"formula_a": fa, "formula_k_Mpc_inverse": fk,
                   "formula_ccl_inputs": np.asarray(fields).transpose(1, 0, 2),
                   "formula_ccl_response_Mpc3": direct[0]})
    record["formula_diagnostic"] = {
        "scope": "common CCL ingredients; response conventions only",
        "units": "Mpc^3 for power, I02, I12 and response; I11/slope unitless",
        "input_axis_order": "[quantity,a,k]",
        "input_roles": ["Plinear", "Pnonlinear", "I11", "I02", "I12",
                        "dln(I11^2 Plinear)/dlnk", "dlnPlinear/dlnk"],
        "centered_logk_step": step,
    }
    record["status"] = "completed"
    save_stage(args.output, record, arrays)


def project_cocoa(args, record):
    """Use CoCoA C kernels for controlled swaps and shared-input projection."""
    import numpy as np

    exported = load_bundle(args.ccl_export, "tjpcov-ssc-model-export-v1")
    cocoa = load_bundle(args.cocoa_ssc, "cocoa-ssc-shear-v1")
    if (exported["status"] != "completed"
            or exported["cocoa_manifest_sha256"]
            != sha256(args.cocoa_ssc / "manifest.json")):
        raise ValueError("CCL export must use this completed CoCoA SSC run")
    project = args.cocoa.resolve() / "projects/lsst_y1"
    sys.path.insert(0, str(project))
    import cosmolike_lsst_y1_interface as ci
    if sha256(ci.__file__) != cocoa["interface_sha256"]:
        raise ValueError("CoCoA binary changed since the saved native SSC run")
    backend = ci.covariance
    original = np.load(args.cocoa_ssc / "ssc.npz", allow_pickle=False)
    supplied = np.load(args.ccl_export / "models.npz", allow_pickle=False)
    ell = original["ell"]
    a, _, distance, dchi = original["geometry"]
    np.testing.assert_array_equal(ell, supplied["ell"])
    np.testing.assert_array_equal(a, supplied["a"])
    length = exported["length_unit_Mpc"]
    spin = np.sqrt((ell-1)*ell*(ell+1)*(ell+2))/(ell+0.5)**2

    def shell_response(distances, window, response, shear_spin):
        # Each C shell response is W^2 D / chi^2, with two shear legs.
        # The C projection then multiplies two shells through dchi sigma_B^2.
        pair = np.ascontiguousarray(np.broadcast_to(window**2, response.T.shape))
        return backend.covariance_ssc_shell_response(
            distance=np.ascontiguousarray(distances), signal=np.zeros(len(ell)),
            pair_window=pair, mean_window=np.zeros_like(pair),
            power_response=np.ascontiguousarray(response.T*shear_spin[:, None]**2))

    def project_shell(shell, weight):
        return backend.covariance_project(
            left=shell, right=shell, weight=np.ascontiguousarray(weight))

    record.update({
        "ccl_export_manifest_sha256": sha256(args.ccl_export / "manifest.json"),
        "cocoa_manifest_sha256": sha256(args.cocoa_ssc / "manifest.json"),
        "interface_sha256": sha256(ci.__file__),
        "cases": {}, "units": {"covariance": "dimensionless"},
        "scope": "SSC-only controlled predictions; no G/cNG or band averaging",
    })
    arrays = {key: supplied[key] for key in supplied.files}
    shells = {}
    # The first case exactly reconstructs the saved CoCoA calculation.
    # Each later case swaps only a supplied table; no halo fit is changed.
    for response_name, values in (
            ("cocoa", original["response"]),
            ("ccl", supplied["ccl_response_at_cocoa_k_Mpc3"]/length**3)):
        shells[response_name] = shell_response(
            distance, original["window"], values, spin)
        for window_name, variance in (
                ("cap", original["variance"]),
                ("disc", supplied["ccl_disc_at_cocoa_a_Mpc"]/length)):
            name = f"{response_name}_response_{window_name}_window"
            matrix = project_shell(shells[response_name], dchi*variance)
            arrays[name] = matrix
            record["cases"][name] = matrix_difference(matrix, original["ssc"])
            save_stage(args.output, record, arrays)
            if name == "cocoa_response_cap_window":
                error = record["cases"][name][
                    "max_abs_difference_over_reference_rms_product"]
                if error > 1e-10:
                    raise ValueError("CoCoA native reassembly gate failed")
                print("Native CoCoA reassembly passed", flush=True)

    # Obtain exactly CoCoA's precomputed GSL rule and map it to the saved
    # panels. Only its quadrature machinery is retained for this check;
    # distances, H/H0, tracer kernel, response and window all come from CCL.
    settings = cocoa["settings"]
    rule = backend.covariance_integration_rule(nquad=settings["radial_nquad"])
    edges = settings["a_edges"]
    nodes = np.concatenate([(left+right)/2+(right-left)*rule[0]/2
                            for left, right in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(right-left)*rule[1]/2
                              for left, right in zip(edges[:-1], edges[1:])])
    np.testing.assert_allclose(nodes, a, rtol=5e-14, atol=0)
    ccl_dchi = weights/(a*a*supplied["ccl_E"])
    shared_shell = shell_response(
        supplied["ccl_chi_Mpc"]/length,
        supplied["ccl_kernel_Mpc_inverse"]*length,
        supplied["ccl_response_at_ccl_k_Mpc3"]/length**3,
        supplied["ccl_spin"])
    common = project_shell(
        shared_shell, ccl_dchi*supplied["ccl_disc_at_cocoa_a_Mpc"]/length)
    arrays.update({"cocoa_projection_of_ccl_inputs": common,
                   "shared_ccl_shell_core_units": shared_shell,
                   "shared_ccl_dchi_weights_Mpc": ccl_dchi*length,
                   "cocoa_shell_native": shells["cocoa"],
                   "cocoa_shell_ccl_response": shells["ccl"],
                   "cocoa_dchi_weights_Mpc": dchi*length})
    record["projection_check"] = matrix_difference(
        common, supplied["ccl_dense_window_ssc"])
    record["remaining_after_response_and_window_swaps"] = matrix_difference(
        arrays["ccl_response_disc_window"], supplied["ccl_dense_window_ssc"])
    # A few near-observer samples extrapolate to extreme k and can have
    # negative CCL responses. Keep them in the full result. Measure their
    # entire shell contribution with the same C projector, to distinguish
    # a large unweighted response ratio from a large covariance effect.
    negative = np.any(supplied["ccl_response_at_ccl_k_Mpc3"] < 0, axis=1)
    contribution = project_shell(
        shared_shell, ccl_dchi*supplied["ccl_disc_at_cocoa_a_Mpc"]
        /length*negative)
    arrays["negative_response_shell_contribution"] = contribution
    rms = np.sqrt(np.diag(common))
    record["negative_response_shells"] = {
        "count": int(negative.sum()),
        "maximum_redshift": float(np.max(1/a[negative]-1)) if negative.any() else None,
        "max_abs_contribution_over_common_rms_product": float(np.max(
            np.abs(contribution/rms[:, None]/rms[None, :]))),
        "scope": "All shells with any negative CCL response; retained in full result",
    }
    record["interpretation"] = {
        "response_swap": "native responses include different halo fits and readers",
        "window_swap": "native window predictions include power/background choices",
        "projection_check": "same CCL physical inputs; GSL versus CCL radial "
                            "integration and interpolation remain distinct",
        "remaining": "CoCoA versus CCL geometry/tracer weights and projection; "
                     "not uniquely a quadrature error",
        "interaction": "response and window changes multiply; they need not add",
    }
    save_stage(args.output, record, arrays)
    print("Saved the four swaps and common-CCL-input projection", flush=True)

    # The C response kernel accepts dimensional inputs in any consistent
    # length unit. Keep Mpc^3 here, avoiding unnecessary conversion rounds.
    common_inputs = supplied["formula_ccl_inputs"]
    shape = common_inputs.shape[1:]
    flattened = common_inputs.reshape((7, -1))
    linear_rule = flattened[[0, 1, 2, 3, 4, 6]].copy()
    linear_rule[2] = 1.0  # CCL's response helper uses Plinear, not I11^2 Plinear.
    raw_rule = np.ascontiguousarray(flattened[:6])
    for name, values, fractional in (
            ("ccl_convention", linear_rule, False),
            ("two_halo_convention", raw_rule, False),
            ("fractional_nonlinear_convention", raw_rule, True)):
        result = backend.covariance_halo_response(
            inputs=np.ascontiguousarray(values), growth_coefficient=47/21,
            dilation_coefficient=1/3, fractional=fractional)
        arrays[f"formula_{name}_Mpc3"] = result[1].reshape(shape)
    expected = supplied["formula_ccl_response_Mpc3"]
    gate = float(np.max(np.abs(arrays["formula_ccl_convention_Mpc3"]/expected-1)))
    record["formula_gate_max_fractional_difference"] = gate
    save_stage(args.output, record, arrays)
    if gate > 1e-10:
        raise ValueError("Shared-ingredient CCL response convention gate failed")
    record["formula_scope"] = exported["formula_diagnostic"]
    record["status"] = "completed"
    save_stage(args.output, record, arrays)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest="stage", required=True)
    ccl_parser = stages.add_parser("ccl", help="export in the TJPCov environment")
    for name in ("inputs", "gaussian_run", "native_ssc", "cocoa_ssc"):
        ccl_parser.add_argument(name, type=Path)
    ccl_parser.add_argument("--tjpcov", type=Path, required=True)
    cocoa_parser = stages.add_parser("cocoa", help="project in the CoCoA environment")
    cocoa_parser.add_argument("ccl_export", type=Path)
    cocoa_parser.add_argument("cocoa_ssc", type=Path)
    cocoa_parser.add_argument("--cocoa", type=Path, required=True)
    for stage in (ccl_parser, cocoa_parser):
        stage.add_argument("--output", type=Path, required=True)
        stage.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    threads = require_thread_environment()
    if threads > 6 or args.timeout <= 0:
        parser.error("Need at most six threads and a positive deadline")
    if args.output.exists():
        parser.error("Choose a fresh output directory")
    args.output.mkdir(parents=True)
    record = {
        "schema": ("tjpcov-ssc-model-export-v1" if args.stage == "ccl"
                   else "cocoa-tjpcov-ssc-model-diagnostic-v1"),
        "status": "started", "threads": threads,
        "timeout_seconds": args.timeout, "script_sha256": sha256(__file__),
        "timing_policy": "Correctness only; no benchmark durations",
    }
    write_json(args.output / "manifest.json", record)
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(args.timeout)
    try:
        (export_ccl if args.stage == "ccl" else project_cocoa)(args, record)
    except Exception as error:
        record["status"] = "failed"
        record["error"] = f"{type(error).__name__}: {error}"
        write_json(args.output / "manifest.json", record)
        raise
    finally:
        signal.alarm(0)
    print(f"Saved {args.stage} SSC model diagnostic to {args.output}")


if __name__ == "__main__":
    main()
