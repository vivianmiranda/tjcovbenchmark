# Trispectrum and SSC attribution, 2026-10-06

This is the historical 1,500-node comparison before adoption of global
natural-cubic power refinement. Follow [the refresh record](global_power_refresh.md)
for the new production campaign; do not reuse these numerical conclusions
as its results. This benchmark used the guarded
Wynn/bias-0.8 binary with SHA256
`22e2e3499ac4de7bdcd52eaee62fb148a30d3773d1c49e684d827e8943d079cc`.
CCL is 3.3.3; TJPCov is commit
`2f59302af33607aec712185632d8274e59e6b33c`.
No code was patched to improve cross-code agreement. Runs used at most
six OpenMP threads and one BLAS thread, sequentially. Timings are pending.

## Separated trispectra

Scripts export actual CoCoA C kernels and actual public CCL functions
called by TJPCov: 1h, 2h13, 2h22, 3h, 4h. The common grid is nine k nodes
from 0.001 to 10 h/Mpc, all 45 unordered pairs, z=0.1,0.5,1. CCL arrays
are [a,k2,k1]; both I13 orientations are checked in shared-input assembly.
Trispectrum units are (Mpc/h)^9. These are not projected covariances.

Native maximum 4h discrepancy is about 61.6%, at K=.001,Q=.316 h/Mpc.
At z=1 this term contributes about 61% of the CCL sum: it is not merely
a negligible-term ratio. Largest summed-term difference is 38.81%.

Using identical CCL moments and power samples in CoCoA kernels gives
agreement within 4e-16 (2h13), 4.7e-9 (2h22), 4.1e-9 (3h), and 5e-6 (4h).
CCL direct-growth mode is used for this algebra comparison. The 1h input
is copied, so its equality is only a consistency check. Changing native
CCL separable growth alone changes 4h by at most .119%.

Native linear P differs by only .0114%, but cancellation in the squeezed
tree terms amplifies its local interpolation structure. Cubic-filling
the same CAMB samples into nested denser tables, with the actual linear
CoCoA reader and fixed moments, reduces max 4h differences from a smooth
input control: 61.63%,19.39%,1.73%,.4875%,.3883% at
1500,2999,5997,11993,23985 nodes. This is not a new CAMB solution, nor
proof that the last grid is converged. No production refinement adopted.

Boundary caveat: the dense end intervals alter extrapolation secants;
the direct cubic control retains original slopes. The archive records
both all-pair and interior-only residuals. All 36 interior pairs stay
inside CAMB support at every angular sample, and every maximum above
occurs among them. Never conflate the +61.6% CoCoA/control ratio with
the approximately -38.1% control/CoCoA ratio.

At z=.5,k=10, native 1h differs by26.15%. Changing CCL concentration
from Duffy08 to Bhattacharya13 reduces this to2.81%; changing abundance
as well gives -7.93%. Fit names do not fix multiplicity normalization,
growth or collapse-threshold conventions. These are sensitivities, not
evidence for preferring one fit. CCL128→255 mass nodes changes terms by
at most .00312%; CoCoA96→128→256 rules change them by .000303%/.000207%.

Files: scripts/export_cocoa_trispectrum.py, run_trispectrum.py,
compare_trispectrum.py, diagnose_trispectrum_power.py, plot_trispectrum.py.
Verified arrays/manifests are in results/trispectrum; four figure pairs
are in figures/trispectrum. Source-manifest hashes link the exported
CCL internal powers to the exact CoCoA quadrature. Plotting checks
coordinates, finite values, native baseline and positivity on log axes.

## SSC attribution

scripts/diagnose_ssc_models.py has CCL and CoCoA stages. It reproduces
native CCL SSC, exports response/window/background/tracer ingredients,
and calls actual CoCoA C projection for table swaps. CoCoA's native
reassembly is bitwise exact. This is point-sampled five-ell shear SSC,
not the full band-averaged survey covariance.

At fixed CoCoA geometry/tracers/projection, replacing only the window
prediction lowers diagonal SSC by8.55–8.83% at low ell,7.27–7.30% high;
response-only lowers1.86–19.66% low,13.78–13.99% high. Both lower
10.33–26.49% low,20.06–20.28% high. Changes multiply, not add. Window
swap includes cap/disc plus P/background readers; response swap includes
halo fits/readers and response prescription. Neither is a single-factor
physical calibration.

Both swaps leave max .187% low/.351% high relative to CCL SSC rms.
Using CCL physical inputs throughout with CoCoA GSL projection leaves
.00204%/.000251%. CCL variance regridded to include GSL nodes changes
its own projection by .00152%/.000309%, recorded separately.

The CCL shear transfer prefactor uses piecewise large-ell approximations;
do not equate it exactly to sqrt[(ell-1)ell(ell+1)(ell+2)]. The diagnostic
uses public get_f_ell and preserves this difference. The original failed
exact-prefactor gate remains in work/ssc_models_low_ccl; the successful
CCL exports are low_ccl_v2 and high_ccl.

A few near-observer extrapolated CCL responses are negative. They are
not dropped: an actual C projection of all affected shells contributes
at most3.38e-18 low/1.22e-12 high relative to common SSC rms. Full arrays
preserve their signs. Main plots explicitly show z=.01–2.

Using common CCL halo moments/powers, CoCoA response kernel reproduces
the CCL linear-response convention to6.7e-16. Switching to the two-halo
slope/amplitude and then transferring to nonlinear power gives the
separate convention curves. At k10, combined change is7.40%,8.28%,12.69%
for z.1,.5,1. These are local responses, not projected SSC percentages.

## Galaxy-bias placement

Source audit finds bias in NumberCountsTracer plus explicit bias product
in TJPCov higher-halo cNG; SSC response is also bias-weighted before
projection through biased tracers. scripts/diagnose_galaxy_bias.py uses
the actual builder and public CCL projectors. Four galaxy legs with
b1→b2 give16; adding TJPCov's extra higher-halo multiplier gives256;
shear control1. SSC holds its actual bias/subtraction-aware response
fixed while changing projection tracers, giving another factor16.
This is a placement test with supplied diagnostic inputs, not a native
full physical matrix. Shear comparisons are unaffected. No issue filed
or library changed; report the inspected revision, not all TJPCov releases.

## Remaining gates

- Carry interpolation changes through actual cNG projections and complete
  small G/SSC/cNG/total matrices; do not extrapolate a 4h percentage to
  total covariance accuracy. Keep native and shared-ingredient results.
- Check consistent Gaussian band averaging versus SSC/cNG ell-center
  evaluation before a total comparison.
- Investigate real-space noise/capability limits separately.
- Quiet timings last. No construction runtime conclusion from these
  supplied-input response and assembly tests.
