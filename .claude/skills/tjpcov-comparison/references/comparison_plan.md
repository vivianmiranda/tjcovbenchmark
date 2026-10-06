# Comparison plan and parity with OneCovariance

Added as item 5 in the overall work order, 2026-10-06. Source review,
README and environment preparation may proceed while CoCoA validation
runs. Numerical TJPCov comparisons wait for a validated, committed Wynn
baseline and the current OneCov figures/README refresh. Do not interrupt
the single active numerical job or install into its environment.

## Common inputs

Export inputs from CoCoA's LSST Y1 covariance adapter. Start with the
OneCov study's source bin 3 and lens bins 1 and 2 (one-based survey bin
identities, not new labels for a truncated catalog). Preserve full n(z)
support and the original midpoint convention.

| Quantity | Initial case |
| --- | --- |
| Area | 12300 deg^2 |
| Source density | 2 arcmin^-2 in the selected bin |
| Lens density | 3.6 arcmin^-2 per selected bin |
| Shape noise | 0.26 per ellipticity component |
| Lens biases | 1.72716 and 1.65168 |
| Cosmology | Omega_m=0.3, Omega_b=0.05, h=0.7, n_s=0.965, A_s=2.1e-9 |
| Initial modeling | Massless neutrinos, Limber, zero IA/RSD/magnification |

Verify these values against the saved export, not just this plan. Preserve
input hashes, component order, units, code revisions and package versions.
Shared CAMB tables remove a background/power ambiguity; shared tables do
not make halo abundance or SSC response prescriptions identical.

## Ordered tests

1. **Environment and input adapter.** Install only with authorization in
   the separate Conda/.local environment. Check imports, freeze resolved
   versions and create one SACC source-bin case. Validate counts, units,
   windows, ell nodes, ordering and the actual CCL object consumed by TJPCov.
   Confirm a minimal native Gaussian call before extending the runner.

2. **Gaussian assembly.** Reproduce the OneCov five-band shear case and
   30x30 two-lens/one-source case for integer ell 30--149 and 1500--1619.
   TJPCov f_sky recomputes C_ell: export those actual spectra to CoCoA and
   match the actual discrete TJPCov bin operator. Test signal, mixed and
   pure noise separately using supported inputs; document extraction if
   repeated native calls are needed. Do not substitute SACC means for
   spectra or secretly intercept CCL calls. Follow with native spectra.

3. **SSC.** Compare background variance, responses and radial projection
   separately. First use common supplied response/window tables in the
   real CoCoA and CCL projection APIs; label this CCL/TJPCov ingredient
   diagnostic. Then run the native TJPCov SSC class. Probe redshifts 0,
   0.5 and 1, k=0.001--10 h/Mpc, low/high ell and off-diagonal terms.
   Explain circular-disc versus spherical-cap variance before interpreting
   a difference. Check local-density subtraction and where each galaxy
   bias enters. Refine native integration and table controls separately.

4. **Halo ingredients.** Repeat sigma(M), dn/dlnM, bias, concentration,
   NFW, I11, I02 and I12 checks at z=0.1, 0.5 and 1; masses 1e10--1e15
   Msun/h and k=0.001--10 h/Mpc. Inspect CCL's mass-integration completion,
   mass definition and profile normalization. Distinguish native Tinker08
   versus CoCoA Tinker10 abundance from a numerical integration effect.
   Audit sigma power range, extrapolation and derivatives explicitly.

5. **Trispectrum orders.** Use the same nine k nodes over 0.001--10 h/Mpc,
   all 45 unordered pairs and the three redshifts. Retain 1h, 2h13, 2h22,
   3h and 4h separately. Export the exact public CCL calls used by TJPCov.
   Compare native ingredients and a common-input diagnostic only where
   that API supports it. Check unequal K,Q, both 1+3 partitions, pair
   exchange and equal-pair/internal-q limits. Establish angular and
   mass convergence independently; do not inherit OneCov findings.

6. **cNG projection.** Repeat the 8x8 shared-trispectrum projection and
   native five-band comparisons, using CCL Tk3D and angular_cl_cov_cNG.
   Match Mpc/h versus Mpc units and all powers of density and bias.
   Native TJPCov SSC/cNG evaluate ell centers, while its Gaussian path
   averages bands. Compare that native convention explicitly, then a
   separately labeled shared-estimator result. Do not silently mix them.

7. **Complete Fourier matrix.** Repeat the source-bin 100x100 case at
   ell=30,60,...,3000. Keep G, SSC, cNG and total separate, retaining every
   entry. Establish the Gaussian/non-Gaussian estimator convention before
   claiming the result is a uniformly band-averaged covariance. Expand to
   selected 3x2pt only after documenting the galaxy HOD mismatch.

8. **Real space.** Target the same eight annuli over 2.5--250 arcmin and
   the full 16x16 xi+/xi− matrix. Begin with native Gaussian covariance.
   Audit omitted B-mode noise, bin reconstruction/weights, Wigner/Bessel
   switch, worker count and memory first. Use each code's actual controls
   and demonstrate its own cutoff convergence. Native SSC/cNG are absent;
   a later projection of TJPCov Fourier terms would be a named diagnostic,
   not a native TJPCov real-space total. Do not fill absent components
   with zeros or invent an independent transform for the comparison.

9. **Timing.** Once each component's accuracy/modeling is understood,
   measure equivalent work sequentially. Separate initialization, shared
   halo tables, projection and assembly; retain full construction totals.
   No warm cache or precomputed data may disappear from a full-runtime
   claim. Use fresh outdirs to avoid filename-based block reuse. Record
   native worker counts, not just OMP_NUM_THREADS. Repeated means/scatter
   go into tables. Only estimate full-survey cost after measuring shared
   and incremental costs; a one-bin time alone is insufficient.

## Matrix diagnostics and reporting

- Save every matrix component in machine-readable form before plotting.
  Common order and input hashes must be asserted, not assumed.
- Use a four-panel G/SSC/cNG/total difference plot like OneCov. Normalize
  differences by sqrt(Cref_total[ii]*Cref_total[jj]), state the sign and
  reference, and retain zero-valued component differences explicitly.
- Also compare component norms, diagonal variances, correlations, total
  positivity and generalized eigenvalues against a positive reference.
  Report uncut and survey-cut results separately; never hide negative modes.
- Small-block convergence does not certify every LSST bin or a Fisher
  forecast. Fisher checks require specified parameters and derivatives.
- README: current measured results only, plots for science and tables for
  times. Unsupported/pending comparisons are named honestly. Historical
  cutoff investigations remain in skill references.

## README and figure template

Use the current sibling OneCov-benchmark-/README.md and its plotting
scripts as the presentation template. Keep the same progression through
Gaussian, SSC, halo ingredients, sigma/abundance, separated trispectra,
cNG projection, complete Fourier and real-space covariance results.
Use the same Cocoa-style installation and reproduction steps.

| OneCov figure | Corresponding TJPCov figure |
| --- | --- |
| gaussian_matrices.png | Correlation matrices and magnified residuals for the same selected bins. |
| gaussian_components.png | Gaussian signal, mixed-noise and pure-noise variance contributions. |
| ssc_projection.png | Shared-input SSC projection and matrix differences. |
| halo_ingredients.png | Native halo ingredients and their ratios at matching redshifts. |
| sigma_abundance.png | Mass variance and abundance differences, with their causes tested separately. |
| trispectrum_terms.png | Separate 1h, 2h, 3h and 4h terms, including unequal wavenumbers. |
| connected_projection.png | Shared-trispectrum cNG projection and residuals. |
| complete_shear_difference.png | Four-panel G/SSC/cNG/total Fourier covariance differences. |
| real_shear_difference.png | The complete selected xi+/xi− matrix, including their cross-covariance, for implemented components. |

Retain panel order, probe separators, units, legend conventions and
normalization across the two studies. Set residual color ranges from
the actual results and label them; matching appearance must not hide
different discrepancy sizes. Captions must name the compared codes,
inputs, difference sign, reference normalization and scope.

Only show panels supported by actual TJPCov calculations. In particular,
do not mimic a four-component native real-space result where SSC/cNG
are unavailable. Explain that capability difference beside the figure.
Keep timings in tables. Add each plot to the README when its comparison
is validated; do not fill it with prospective or copied numerical results.

## Next concrete action

After dependency installation is authorized and the numerical slot is
free, validate the environment and implement the SACC/export adapter plus
the smallest Gaussian case. No numerical result exists yet.
