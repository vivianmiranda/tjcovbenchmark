# Projected cNG and complete Fourier covariance — 2026-10-07

## Accepted scope and archives

The current source-bin pilot compares actual TJPCov public Gaussian,
SSC and cNG calculators with CoCoA production covariance C kernels.
No CoCoA, CCL or TJPCov scientific source was modified for this comparison.
Use the 11,993-node global power/Wynn baseline, not historical 1500-node
power or retired selective-4h results.

Published results:

- `results/complete_fourier/report.json` embeds all source manifests,
  input/code hashes, component diagnostics and refinement comparisons.
  `matrices.npz` preserves all 100x100 entries for G/SSC/cNG/total.
- `results/cng_projection/manifest.json` and `comparison.npz` preserve
  the shared-input 8x8 projection comparison. `inputs_manifest.json` and
  `projection_inputs.npz` preserve the sampled physical inputs. The large
  native trispectrum archive is identified by its original manifest/hash,
  not duplicated in the compact published record.
- `results/cng_time_grid/manifest.json` and `time_grid.npz` preserve all
  entries of the native 99-to-197 time-grid control at fixed 234 k nodes.
  `assembly_verification.json` records actual C-kernel replays.
- Figures are under `figures/complete_fourier/` and
  `figures/cng_projection/`. The classic four-panel figure uses separate
  G/SSC/cNG/total color scales and retains every matrix entry. Additional
  figures show native cNG structure and diagonal component contributions.

The local campaign is `work/projected_cng_20261007/`. Native Gaussian/SSC
come from `native_k32_a2`; the accepted native cNG component comes from
`native_k96_a2`. The collector verifies their common input fingerprint,
ell, edges, signal, source code and scientific settings before combining
them. CoCoA cases are `cocoa_i0` and `cocoa_i1`. These are accuracy runs,
not quiet timing measurements; native cNG table capture uses read-only
return-frame observation, which is explicitly excluded from timing use.

## Physical inputs and estimator

Common input manifest SHA256:
`50ec38096881eea6de1e27ae06c187c86f4d3f931d9b39ef508f6f8817bdda33`.
The original bundle is `work/global_power_11993/lsst_y1`.
Use LSST Y1 source bin 3, area 12300 deg2, source density 2 arcmin^-2,
shape dispersion 0.26 per component, massless neutrinos, Limber and zero
IA, magnification and RSD. CoCoA initialization must reproduce the saved
CAMB power inputs exactly. Each native calculation retains its own
background, power interpolation, halo model and survey-window prediction.

- **Gaussian:** integer multipoles with edges 15,45,...,3015 and weights
  proportional to ell. TJPCov's native binning excludes the upper edge.
  SACC windows contain that endpoint so TJPCov can reconstruct the edge;
  it does not receive nonzero covariance weight. Native Gaussian uses
  CCL's default N_K=167 (1220 k nodes) and 99 a-grid nodes here.
- **SSC and cNG:** point-evaluated at ell=30,60,...,3000. Do not describe
  these matrices as band averages or silently apply the Gaussian operator
  to them. This mixed estimator is the tested TJPCov convention.
- **Signals:** compare native `signal` to CoCoA `signal_centres`.
  CoCoA's separate `signal` is the Gaussian band mean and is not the same
  estimator. The collector checks center signals separately.
- **SSC:** retain the previously tested native a refinement 16 (785 a
  nodes) and default 1220 k nodes. Its finite-disc background variance
  remains distinct from CoCoA's spherical-cap prescription.

This is one 100x100 source-bin covariance, not a complete 1560x1560 LSST
Y1 matrix, a native real-space total, or a calibrated likelihood covariance.

## Native cNG sampling and controls

The full physical k interval remains 5e-5 to 1000 Mpc^-1 in every native
cNG control. N_K is CCL's sampling-density control, not the resulting
table length. Explicit N_K=16,32,64,96 gives 117,234,468,701 k nodes.
The accepted pilot uses N_K=96/701 nodes; CCL's installed default is
N_K=167/1220 nodes. Do not label the pilot as default sampling or publish
its elapsed time as the default runtime. No cutoff was reduced to make
the pilot fit laptop memory.

The a-refinement-2 full grid has 99 nodes, with 69 retained by the native
tracer-support selection. The a-refinement-4 control has 197 full nodes.
All native calls retain `qag_quad` and public CCL numerical defaults
unless the individual manifest explicitly identifies a diagnostic.

| k nodes, coarse to fine | Maximum cNG change / fine cNG rms product | Maximum total change / fine total rms product | Maximum total mode-variance change |
| --- | ---: | ---: | ---: |
| 117 to 234 | 5.12795% | 0.0301142% | 0.0426291% |
| 234 to 468 | 3.35706% | 0.0167665% | 0.0170107% |
| 468 to 701 | 1.11497% | 0.00435435% | 0.00443272% |

Keep native G and SSC fixed while assessing these cNG refinements. The
last cNG diagonal change alone is 0.4164%; the 1.11497% figure includes
all off-diagonal entries. It is a maximum rms-normalized entry, not the
RMS of the residual matrix. At fixed 234 k nodes, 99-to-197 a nodes changes
cNG by at most 0.0992720% of the fine cNG diagonal rms product. That time
control does not establish convergence at the final 701-node k grid.

CoCoA i0-to-i1 (96-to-128 point line-of-sight/integral rules, with power
tables fixed) changes total rms-normalized entries by at most 0.000240107%
and total generalized mode variances by at most 0.000685918%. Its G, SSC
and cNG maximum total-rms changes are respectively 9.34058e-6%,
0.000229820% and 1.61442e-6%.

The cNG sampling result is useful for the *total covariance* but does not
establish percent-level convergence for every cNG entry or all cosmologies.
Do not weaken tolerances or replace these distinctions with a generic
"converged" label.

## Production assembly checks

`run_cocoa_complete.py` calls the actual Gaussian Fourier and covariance
projection C APIs. Its matter helpers retain a 100x100-by-shell cNG
array (about 54 MB at 672 shells). The final one-probe contraction uses
the C weighted projector over flattened matrix entries, avoiding a
four-probe allocation roughly sixteen times larger. This changes storage,
not the physics or projection weights.

`verify_complete_assembly.py` independently replays the production C
Gaussian kernel on TJPCov's actual spectra and the dedicated connected
projection kernel on saved CoCoA samples:

- Shared-spectrum Gaussian maximum rms residual: 4.37225e-16.
- Dedicated connected-projector replay residual: 3.37669e-16.

These are checks of actual code paths, not a public third-party reference
implementation. They separate assembly agreement from native inputs.

## Complete native comparison

| Component | Maximum absolute difference / TJPCov total rms product | Fractional diagonal range relative to that component |
| --- | ---: | ---: |
| Gaussian | 0.563153% | +0.01609% to +0.57665% |
| SSC | 0.748644% | +9.60162% to +41.27150% |
| cNG | 0.0153905% | -0.869402% to +18.66481% |
| Total | 1.132167% | +0.0194833% to +1.132167% |

Native center signals differ by at most 0.470280%, which explains why the
native Gaussian comparison cannot inherit the earlier shared-spectrum
roundoff conclusion. The native halo/window differences remain physical
and numerical modeling differences, not solely interpolation error.

Both total symmetric parts pass Cholesky. Generalized CoCoA/TJPCov total
variance ratios are 1.0000649195 to 1.0823604380, so the largest coherent
mode changes variance by 8.23604%. This exceeds the largest diagonal
change and is mainly SSC; the isolated SSC difference contributes up to
7.77162% in a total-whitened mode. Do not equate this mode test with a
parameter-level Fisher constraint. No cNG positivity requirement is
imposed. Native antisymmetries around 8e-18 on the total rms scale are
preserved in the archives; only eigenvalue diagnostics use symmetric parts.

## Shared-input cNG projection

Final inputs are `shared_final_i2` and comparison `projection_final_i2`.
They use the accepted native 701-k-node CCL Tk3D, CCL distances/background,
CCL source kernel and the same spin factors at eight multipoles:
30,450,870,1290,1710,2130,2550,3000. Native Tk3D replay is bitwise equal
before projecting. The actual trispectrum samples have axes (a,ell2,ell1)
and must be transposed to (ell2,ell1,a) for the C projector.

The final CoCoA grid reuses the actual GSL 256-point nodes from the saved
SSC integration-level-2 grid, seven panels and 1792 nodes in total. Panel
redshift edges are 3.5,2,1.5,1,0.7,0.4,0.2,0.00001. Do not substitute
NumPy-generated quadrature nodes while claiming to reproduce this grid.
The maximum difference is 0.0196581239% of the CCL cNG diagonal rms
product, with all 64 entries retained; RMS residual is 0.00467131%.

This is a common-physical-input finite-endpoint and resolution check:

- Native CCL tracer support in chi is [0,3393.33681165] Mpc.
- CoCoA panels span chi [6794.49762417,0.0428268375] Mpc in integration
  order. The CCL kernel is zero above its actual support, but CoCoA omits
  the observer interval below z=1e-5. Do not claim identical domains.
- CCL Tk3D is piecewise linear in a, introducing derivative discontinuities
  inside fixed CoCoA panels. Earlier 96/128-node controls were not monotonic
  (about 0.064/0.074%); the 256-node result does not prove arbitrary-grid
  convergence. Keeping the interpolation knots explicit matters.
- A previous general `EPSREL` diagnostic also changed background/kernel
  interpolation and cannot be described as a projection-only refinement.
  Future isolated controls may vary only the public Limber tolerance;
  require exported kernel/background equality before attributing effects.

## Reproduction and remaining work

The README gives the four commands for a fresh native complete pilot.
Optional `--cocoa-refined` and repeated `--cng-refinement` arguments to
`report_complete_covariance.py` collect the corresponding controls.
`report_cng_time_grid.py` publishes the separate time-grid comparison.
`plot_cng_projection.py --comparison ... --export ... --output ...
--figures ...` preserves shared inputs and produces the projection figure.
Read each script's arguments and use fresh directories to avoid native
cache hits. Only the coordinating parent starts numerical jobs.

Next measured scopes are quiet complete-pilot/Gaussian/halo timing,
full-survey expansion, and native real-space Gaussian after noise,
worker-count and memory audits. The inspected TJPCov path still has no
native real-space SSC/cNG. Do not extrapolate full survey speed ratios
from matrix dimensions or reuse OneCov's scaling. Keep the established
current CoCoA Wynn/global-power prescription unchanged.
