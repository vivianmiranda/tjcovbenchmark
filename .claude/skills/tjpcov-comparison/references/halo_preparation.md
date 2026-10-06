# Native halo ingredient comparison — 2026-10-06

This stage uses the same LSST CAMB input bundle as the Gaussian and SSC
pilots. It calls public CCL routines with the choices selected inside
TJPCov's `FourierSSCHaloModel.get_covariance_block`, and the actual CoCoA
halo readers. It does not implement a third halo model or a trispectrum.

## Bounded calculation and provenance

The output grid has three redshifts (0.1, 0.5, 1), 41 logarithmic masses
from 1e10 to 1e15 Msun/h and 33 logarithmic wavenumbers from 0.001 to
10 h/Mpc. Both endpoints are retained. The scripts require an explicit
OpenMP allocation of at most six threads, one BLAS thread, a fresh output
directory and a positive Unix wall-time limit (600 seconds by default).

`run_halo.py` saves CCL quantities to `halo.npz` and a hash-checked
`tjpcov-halo-ingredients-v1` manifest. It checks the installed TJPCov
selector source against the checkout and records all relevant CCL Python
files, its compiled library, package versions and mass integration controls.
The output contains no fair timing claim.

`export_cocoa_halo.py` initializes the real LSST project, checks its
cosmology, and requires every linear and nonlinear CAMB sample to equal
the saved input bundle bitwise. It uses the native export's exact output
nodes and saves `cocoa-halo-ingredients-v1` with the common input hash,
the native manifest hash, core/project revisions and interface hash.

The existing OneCov refresh archive was inspected first. Its
`halo_cocoa_matched_800` result uses the current Wynn domain and the same
interface hash, and contains the requested moments. However, it lacks a
saved CAMB-table hash and a native-concentration NFW profile. A fresh
small export avoids claiming stronger provenance than that archive has.
The existing archive remains untouched.

## Units and axes

- CCL takes mass in Msun: divide the shared Msun/h mass by h.
- CCL takes comoving k in Mpc^-1: multiply the shared h/Mpc value by h.
- CCL's mass function returns dn/dlog10(M) in Mpc^-3. Divide by ln(10)
  and h^3 to obtain dn/dln(M) in (h/Mpc)^3.
- CCL's analytic Fourier NFW contains the mass M. Divide its output by M
  to obtain the dimensionless u(k|M) used by CoCoA.
- Its public profile normalization is the comoving matter density in
  Msun/Mpc^3. Divide I01/I11 by this once and I02/I12 by its square.
  Multiply the latter by h^3 to express them in (Mpc/h)^3.
- The saved density divides CCL's density by h^2, giving
  (Msun/h)/(Mpc/h)^3. CoCoA's own density constant is retained separately;
  it is not silently replaced to improve agreement.

`sigma`, `dndlnm`, `bias` and `concentration` have [redshift,mass] axes;
`linear`, `i01`, `i11`, `i02` and `i12` have [redshift,k] axes;
`profile` has [redshift,k,mass] axes. I01 is exported only by CCL: no
unavailable CoCoA quantity is filled with zero. `rho` has one value per z.

## Distinct native models

- **CCL as selected by TJPCov:** M200m, Tinker08 abundance, Tinker10 bias,
  Duffy08 concentration and analytic, truncated NFW. Its HMCalculator
  defaults to Simpson integration on 128 log10(M) nodes from 1e8 to
  1e16 Msun. A controlled `--mass-refinement 2` uses 255 nodes over the
  same interval; refinement 4 uses 509. These are output-independent
  internal integration nodes. Refinement is a public CCL ingredient
  diagnostic; TJPCov itself does not expose this control in its SSC YAML.
- **CoCoA:** its current bias-normalized Tinker10 multiplicity, fitted
  Tinker10 bias, Bhattacharya13 concentration, M200m NFW and production
  mass panels. Only I11 uses the Wynn tail prescription. The other halo
  moments remain direct integrals on their production finite domain.

CCL's `HMCalculator` computes both a missing mass coefficient and a
missing bias-weighted mass coefficient. Each is the target mean density
minus the resolved mass integral, divided by the minimum mass. Its
integration helpers add the appropriate coefficient times the integrand
at that minimum mass. Consequently this completion enters I01/I11 and
also I02/I12. The fitted abundance and halo bias are not divided by a
finite-range normalization. This is an additive prescription, distinct
from both OneCov's bias rescaling and CoCoA's extended I11 integration.

The native CCL sigma reader builds its sigma spline from the supplied
linear Pk2D through `Cosmology.compute_sigma`; this pilot leaves its
native numerical settings unchanged. Equal CAMB tables do not imply
equal interpolation, extrapolation, integration or growth conventions.
Inspect its sigma integration controls separately before assigning a
measured sigma difference to one cause.

## Diagnostics that isolate a fit choice

CoCoA exports both its native NFW profile and `profile_matched_c`, which
uses CCL's Duffy concentration with CoCoA's actual `u_nfw_c` routine.
The latter isolates profile/radius conventions from the concentration
fit; it must not be labeled a native CoCoA profile.

A further explicit diagnostic can match the radius constant: at fixed
mass R is proportional to rho^(-1/3), so evaluating CoCoA at
`k*(rho_cocoa/rho_ccl)^(1/3)` gives the same kR as CCL. This changes only
the diagnostic argument; native profiles and moments retain their own
physical constants. The coordinator tested this after finding that the
native density constants differ by 5.985e-5 fractionally. Matching both
concentration and radius lowers the maximum absolute profile difference
to 1.86e-9, 2.49e-9 and 3.97e-9 at z=0.1, 0.5 and 1, respectively.
This identifies the density/radius constant as the source of almost all
the discrepancy remaining after concentration alone was matched.

The two Tinker10 bias readers use slightly different collapse constants:
CCL uses 1.68647019984 (its EdS value), while CoCoA uses 1.686. The export
therefore saves two clearly named diagnostics:

- `bias_matched_sigma` gives CoCoA the CCL sigma values, retaining
  CoCoA's own 1.686/sigma peak definition.
- `bias_matched_nu` gives CoCoA the actual CCL peak height. Even then,
  a small difference can remain because delta_c also appears explicitly
  inside the Tinker10 bias formula's denominator.

Neither diagnostic modifies native spectra, abundance or halo moments.
Model differences are not pass/fail accuracy tolerances. They must later
be followed through the separated trispectra and covariance components.

## First runs

In the active TJPCov environment, from this repository:

```bash
python scripts/run_halo.py work/lsst_y1 --tjpcov ../TJPCov --output work/halo_native
```

Then in the CoCoA environment, from cocoa/Cocoa:

```bash
python ../../tjcovbenchmark/scripts/export_cocoa_halo.py ../../tjcovbenchmark/work/lsst_y1 ../../tjcovbenchmark/work/halo_native --cocoa . --output ../../tjcovbenchmark/work/halo_cocoa
```

## First numerical checkpoint

The coordinator executed both exports successfully, including the
bitwise CAMB-table assertion. CCL's mass grid was refined from 128 to
255 nodes; CoCoA was checked at integration levels 0, 1 and 2. The
reported maximum mass-refinement changes are 8.2e-6 for CCL and 8.9e-7
for CoCoA. These small changes do not imply agreement of the halo models.

The current comparison of the refined exports finds maximum fractional
sigma differences below 9.1e-5, abundance differences up to 5.44%, and
I12 differences up to 9.22% over the saved grids. The matched-concentration
NFW discrepancy is at most 1.32e-5 in absolute u; matching the radius
constant reduces this to below 4e-9. The same-peak Tinker bias discrepancy is at most
1.43e-5 fractionally, consistent with the remaining explicit collapse
constant difference requiring its own interpretation.

Final records are in `work/halo_comparison_final/comparison.json`, with
native `work/halo_native[_m255]` and CoCoA
`work/halo_cocoa_final_i0`, `work/halo_cocoa_radius_i1` and
`work/halo_cocoa_radius_i2` exports. Earlier exports remain preserved.
Python syntax, command-line help and diff whitespace checks also pass.
The coordinator owns the comparison/plot script, result publication and
final commit. None of these ingredient
checks certifies full trispectra, covariance matrices or Fisher forecasts.

The comparator requires the exact native manifest used by the CoCoA
export. Final CoCoA exports use the refined 255-node native manifest;
the copied diagnostic inputs are bitwise identical to the 128-node ones.
A deliberate mismatched-manifest call is rejected before publication.
