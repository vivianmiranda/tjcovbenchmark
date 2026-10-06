# Comparison plan and parity with OneCovariance

Current queue item 4. On 2026-10-06 the user authorized numerical
TJPCov comparisons using up to six OpenMP threads, alongside the existing
CoCoA validation. The separate environment is installed and checked.
Correctness runs remain sequential within this comparison; benchmark
timings wait for a quiet machine. No production code is changed here.

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

## Current numerical checkpoint — 2026-10-06

Gaussian shear 5x5 and galaxy/shear 30x30 matrices pass at both low and
high multipoles. Native spectra are exported to CoCoA and native discrete
TJPCov bin weights are matched. Every component agrees within 7e-16 of
the total Gaussian diagonal rms product. All totals are positive.

Nineteen native SSC cases and four public-CCL sampling diagnostics have
completed. See `references/ssc_sampling_diagnostic.md` for the isolated
background-variance sensitivity and remaining model differences. The
complete records and arrays are in `results/ssc_native.json` and `.npz`.
These are 5x5 point-sampled SSC components, not full G+SSC+cNG totals.

Next: native halo ingredients and separated trispectra, then the bounded
cNG and complete-matrix comparisons below. Real-space worker, memory and
noise audits remain gates before a native real-space run. No overlapping
elapsed times from this campaign may be used as benchmarks.

## Gaussian preparation checkpoint — 2026-10-06

`scripts/export_lsst_y1.py` uses the actual project configuration and
initialize functions, saving the source-3/lens-1,2 n(z) columns, survey
numbers and CAMB tables with input units/axes and hashes. It does not
import the OneCov scripts or depend on OneCov's untracked work folders.

`scripts/run_gaussian.py` constructs a SACC source case with five bands,
passes a real CCL CosmologyCalculator to FourierGaussianFsky, and calls
the untouched native calculator. Linear/nonlinear CAMB P are converted
from h/Mpc and (Mpc/h)^3 to CCL's Mpc units; the background remains native
CCL and this distinction is recorded. It checks installed TJPCov source
bytes against the named checkout and checks CCL's Limber default.

The runner exports the CCL C_ell values evaluated with the same cosmology,
tracer, ell nodes and defaults as TJPCov. A SACC mean is never a substitute.
It separates CC/CN/NN by native calls at shape-noise powers 0,N,2N:
NN=(C(2N)-2C(N)+C(0))/2 and CN=C(N)-C(0)-NN. The calls are diagnostic;
their elapsed times must not be put in a production timing table.

SACC endpoint detail: windows include each nominal upper boundary because
TJPCov gets the final edge from the last positive window node. Its native
bin_cov excludes the upper edge. Therefore the exported ell grid includes
one final zero-weight node. The comparison checks actual reconstructed
edges, uses ell weights, and verifies its operator transcription against
native bin_cov with a supplied test array. General windows/noninteger ell
grids are not supported by this first script.

`scripts/compare_gaussian.py` calls CoCoA's actual production interface
with these arrays. It compares all four Gaussian pieces, checks total
positivity and generalized variance ratios, and preserves every entry.
`scripts/plot_gaussian.py` uses the OneCov correlation/residual/component
layout, with actual residual ranges. It will not plot a failed comparison
as a passing result. No scientific figures were produced during the initial preparation.
The completed numerical campaign is summarized above.

The initial import, source-bin and galaxy/shear extension gates are now
complete. Gaussian assembly agreement with shared spectra does not yet
establish agreement of independently generated native spectra.
