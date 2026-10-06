# Two-lens Gaussian extension — 2026-10-06

The user authorized small numerical TJPCov cases, up to six threads,
with timing comparisons deferred until correctness is understood. The
coordinator owns numerical runs. This entry records the source contracts
and the completed shared-input Gaussian checks.

## Scope and public API

`run_gaussian.py --include-lenses` adds the two lens populations exported
from the actual LSST configuration. It preserves the one-source default
and its existing `tjpcov-gaussian-shear-v1` output schema. The new case uses
`tjpcov-gaussian-fields-v1` so a shear-only SSC reader cannot consume the
three-field case silently.

Fields are ordered lens 1, lens 2, source 3 (original survey bin numbers,
with the actual source selection retained if different). The six
observables are g1g1, g1g2, g2g2, g1gamma, g2gamma and gammagamma, each
followed by its five increasing bands. Thus the matrix has all 30x30
entries, including cross-lens clustering and every covariance cross block.

The runner places those observables in SACC and checks its indices. It
uses native `FourierGaussianFsky.get_covariance()` for the three-field
case. That method selects SACC E/density blocks, places every native
block and mirrors the transpose. Do not pass `include_b_modes=False` to
this high-level method: its internal reshape expects the complete native
spin layout before selecting the SACC quantities.

The one-source default still calls the original single-block method with
`include_b_modes=False`. Both modes use the same five bands and native
ell-weighted operator audit described in the earlier source study.

## Exact zero-noise diagnostic without a monkeypatch

The inspected public TJPCov configuration accepts Ngal as a number and
converts arcmin^-2 to sr^-1. There is no finite-value rejection. Its
NumberCountsTracer is constructed from n(z) and bias, without Ngal;
the galaxy shot-noise term is subsequently `1/Ngal`.

- For noise factor 1, keep the physical survey density.
- For factor 2, halve Ngal so that its reciprocal doubles.
- For factor 0, use Ngal=infinity, the public zero-shot-noise limit.
  The runner asserts that TJPCov actually returns exactly zero noise.

The source still uses sigma_e multiplied by sqrt(factor), with its fixed
physical number density. The saved metadata describes the infinite-density
diagnostic but contains only the original finite survey densities. This
is not a change to the galaxies' redshift distributions, biases or native
spectra, and no cached tracer/noise arrays or functions are monkeypatched.

The same three native calls recover CC/CN/NN from the quadratic noise
dependence. Their time is diagnostic, not a production timing. Every
field pair's actual CCL spectrum is exported, with shot/shape noise kept
separate. CoCoA then consumes these arrays through its production C API.

## Reproducing the small cases

After the coordinator frees the numerical slot, in the TJPCov environment:

```bash
python scripts/run_gaussian.py work/lsst_y1 --tjpcov ../TJPCov --include-lenses --output work/gaussian_3x2_low
```

The high-multipole case uses the same LSST bundle:

```bash
python scripts/run_gaussian.py work/lsst_y1 --tjpcov ../TJPCov --include-lenses --ell-range 1500 1620 --output work/gaussian_3x2_high
```

From the active CoCoA terminal in cocoa/Cocoa, for the first case:

```bash
python ../../tjcovbenchmark/scripts/compare_gaussian.py ../../tjcovbenchmark/work/gaussian_3x2_low --cocoa . --output ../../tjcovbenchmark/work/comparison_3x2_low
```

Then plot the passing comparison in the TJPCov environment:

```bash
python scripts/plot_gaussian.py work/comparison_3x2_low --output figures/gaussian_3x2_low
```

Repeat the comparison and plot with high-case paths. The Gaussian runner
has a default 600-second Unix deadline and leaves partial outputs intact.
Use fresh output directories; the script never overwrites an existing run.

The comparator retains all entries and checks component differences,
total positivity and generalized variance ratios. Plotting groups all
five bands of each observable and never connects a component curve across
different spectra. Original five-band results remain readable.

## Validation of the final script

The coordinator reran both three-field cases and the original source-only
case after adding the deadline. Each final-script matrix agrees bitwise
with its earlier saved native calculation. The source-only schema remains
readable by the original comparison and SSC paths.

Both 30x30 CoCoA comparisons pass for sample variance, signal times noise,
pure noise and total, retaining every entry. The largest absolute
difference after division by the CoCoA total rms product is:

| Multipoles | Largest component difference | Largest mode change |
| --- | ---: | ---: |
| 30 to 149 | 6.9604e-16 | 5.1071e-15 |
| 1500 to 1619 | 6.5774e-16 | 1.9985e-15 |

Both total matrices are positive definite. The minimum eigenvalue after
diagonal scaling is about 0.04943 at low multipoles and 0.61911 at high
multipoles. Mode changes use generalized covariance eigenvalues. The
source-only compatibility comparison also passes. These are correctness
measurements; the three noise calls are not a fair production timing.

The records are `work/comparison_3x2_low_checked/comparison.json`,
`work/comparison_3x2_high_checked/comparison.json` and
`work/comparison_shear_checked/comparison.json`; their manifests preserve
the actual script, binary and source hashes. The low-multipole six-spectrum
plot was visually reviewed. Python syntax and command-line help checks
also pass.

No native-spectrum agreement, SSC/cNG, galaxy HOD or Fisher-convergence
claim follows from this Gaussian shared-input extension.
