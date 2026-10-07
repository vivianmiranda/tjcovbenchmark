# Global covariance power-table refresh

The adopted CoCoA covariance configuration subdivides each native CAMB
log-k interval eight times at global accuracy boost one. This turns 1,500
native samples into 11,993 installed samples. Natural cubic splines of
ln(P) versus log10(k) fill the linear, nonlinear and cb tables at fixed
redshift. All covariance readers receive them and continue using their
ordinary linear lookups. No new CAMB information or k support is added.

The older results in `results/` and `work/` retain their original inputs
and hashes. Do not describe their approximately 61.6% maximum native 4h
discrepancy as a result of this new configuration. Recompute before
replacing public results or scientific conclusions. The older power
diagnostic used SciPy's default not-a-knot cubic boundary; the adopted
production helper instead uses natural boundaries.

## Completed refresh

All 63 sequential stages completed in
`work/global_power_11993/refresh_20261007T002618891993Z.json`.
The input export confirms exact retention of all three native CAMB power
tables and exact reproduction of the installed 11,993-node tables by the
production helper. There are 140 redshift nodes. All argument parsers,
script syntax and saved report/file hashes were checked as well.

The current scientific results, figures and README were replaced only
after this successful campaign. The previous public targets and README
are preserved in `work/historical_before_global11993/`; earlier work
directories also remain intact. The publication staging tree is
`work/publication_global_power_11993/`. No performance claim or projected
cNG result is inferred from these 63 stages.

Key refreshed findings:

- Gaussian 5x5 and 30x30 cases pass; their largest component differences
  are 6.58e-16 and 7.67e-16 in total-Gaussian rms normalization.
- CoCoA low-ell SSC integration 0→1 and 1→2 changes are 0.00626% and
  0.00418% at fixed global boost one. CCL 50→99 time nodes still change
  low SSC by 31.5726% in finer-matrix normalization; 393→785 changes
  0.017998%. The main sampling issue remains its background variance.
- Refined native SSC diagonals differ by −26.61% to −10.51% at low ell
  and −20.56% to −20.35% high. Shared physical-input projection residuals
  remain 0.00204%/0.000251%. Both refined symmetric parts are positive;
  high SSC remains nearly rank deficient.
- Native 4h at K=.001,Q=.316 h/Mpc now differs by +0.486925%,+0.464792%,
  +0.430773% at z=.1,.5,1. Native maxima instead occur at high k:
  1.1853%,5.7215%,7.1214%. Their 4h fractions of the CCL sum are only
  0.00902%,0.00244%,0.00604%, respectively.
- At z = 1, K = Q = 10 h/Mpc, the native I11 weight ratio raised to the fourth power
  predicts −7.122249%, compared with actual 4h −7.121364%. This is a
  direct decomposition of saved native arrays, not a new model run.
- Shared CCL powers and moments with direct growth reproduce 4h within
  4.983e-6 fractionally. The fixed-CoCoA-moment power-reader swap changes
  4h by at most 0.4922%. The natural-cubic smooth-control residual is
  0.4875% at 11,993 and 0.3883% at 23,985 on the 36 interior pairs. These
  are distinct comparisons and neither certifies a full cNG matrix.

`export_lsst_y1.py` now initializes the actual forecast twice:

- `power_accuracyboost=1` saves `camb_native_tables.npz`, the native CAMB
  interchange arrays. Global boost stays one.
- The project's production setting saves `camb_tables.npz`, the installed
  tables, and the matching `inputs.npz` supplied to CCL.

The export checks exact equality to
`covariance.power.refine_power_tables(native_tables, refinement)` and
checks that every original node and all three log-power arrays survive
at their nested positions. The manifest records native and installed
node counts, the spline boundary, the helper hash and both archives.
All CoCoA ingredient exports check grid and linear/nonlinear/cb input
equality. No stale-array guard has been weakened.

## Running the prepared campaign

Use an activated CoCoA terminal. Work in `tjcovbenchmark/`. The runner
uses that Python for CoCoA stages and the supplied private Python for
TJPCov stages; never resolve the `.local/bin/python` symlink to its
Conda-base executable. It removes CoCoA's import/library overrides from
the TJPCov child. Existing environment installation is reused.

**Step 1:** inspect the command list without starting either code.

```bash
python scripts/refresh_comparison.py --cocoa ../cocoa/Cocoa --tjpcov ../TJPCov --tjpcov-python .local/bin/python --output work/global_power_11993 --list
```

**Step 2:** export inputs and run the Gaussian pilot cases.

```bash
python scripts/refresh_comparison.py --cocoa ../cocoa/Cocoa --tjpcov ../TJPCov --tjpcov-python .local/bin/python --output work/global_power_11993 --groups inputs gaussian
```

**Step 3:** run the SSC cases and response/window diagnostics.

```bash
python scripts/refresh_comparison.py --cocoa ../cocoa/Cocoa --tjpcov ../TJPCov --tjpcov-python .local/bin/python --output work/global_power_11993 --groups ssc
```

**Step 4:** run halo ingredients and their integration refinements.

```bash
python scripts/refresh_comparison.py --cocoa ../cocoa/Cocoa --tjpcov ../TJPCov --tjpcov-python .local/bin/python --output work/global_power_11993 --groups halo
```

**Step 5:** run separated trispectra and common-input/model diagnostics.

```bash
python scripts/refresh_comparison.py --cocoa ../cocoa/Cocoa --tjpcov ../TJPCov --tjpcov-python .local/bin/python --output work/global_power_11993 --groups trispectrum
```

Omitting `--groups` runs these five groups in order. Each child is limited
to 600 seconds, including scripts without their own alarm, and is fully
stopped before the next starts. OpenMP comes from the environment, at
most six threads, with `OPENBLAS_NUM_THREADS=1`. No timing comparison is
authorized during this correctness refresh: saved stage wall times are
operational diagnostics only.

The runner saves an execution journal and one log per stage. A failed
child, unmet scientific gate, changed script or invalid saved-file hash
stops the campaign and preserves its output. After diagnosis, archive a
failed output folder before an explicit `--start-at STAGE_NAME` restart.
Do not run two runners together or reuse an existing TJPCov cache folder.
No automatic deletion, reference refreeze or tolerance relaxation occurs.

## Scope and interpretation

| Group | Stages | Native output and diagnostics |
| --- | ---: | --- |
| Inputs | 1 | Native CAMB and installed production tables; actual LSST catalogs |
| Gaussian | 12 | Complete 5x5 shear and 30x30 two-lens/one-source cases at low/high ell; signal, mixed and pure noise; matched assembly and figures |
| SSC | 29 | Complete 5x5 shear matrices at low/high ell; native CCL a/k/integration checks, CoCoA integration 0/1/2, response/window swaps, common projection, galaxy-bias placement and figures |
| Halo | 8 | Sigma, abundance, bias, concentration, NFW and I11/I02/I12 at three redshifts; CCL 128/255 mass nodes, CoCoA integration 0/1/2; figures |
| Trispectrum | 13 | Five separated terms, all 45 k pairs at three redshifts; native fits, shared inputs, growth/fit switches, power refinement and figures |

SSC now holds CoCoA's global boost at one while changing the independent
integration level. Its old global-boost-two test cannot silently share
this bundle: global boost two installs 23,985 power nodes and requires a
separately exported, matching CCL input bundle. That is distinct from CCL
refining its internal response-a or response-k sampling on fixed supplied
power tables.

The power interpolation diagnostic now uses the production helper and
natural cubic control. Its default factors one and two are relative to
the installed 11,993-node grid, yielding 11,993 and 23,985 nodes. Both
linear and nonlinear power and cb power are kept compatible. Halo moments
and angular nodes stay fixed only in this diagnostic. The explicit
interior-pair metric and boundary-extrapolation caveat remain necessary.

The refreshed Gaussian and SSC figures show every entry of their small
matrices. Trispectrum figures remain three-dimensional matter ingredients,
not projected cNG. Native projected cNG, complete G/SSC/cNG/total Fourier
matrices and native real-space Gaussian are still unimplemented comparison
stages. TJPCov's inspected native real-space path has no SSC/cNG calculator.
Do not represent any missing component as zero or call Gaussian totals
full covariance totals. Its Fourier Gaussian band average also differs
from native SSC/cNG point sampling at effective ell.

After execution, inspect results and render the scientific figures before
updating the README. Preserve original runs in the historical skill record;
the README should describe the newly measured production configuration.

## Publication preparation

`prepare_publication.py` is a reporting/copying tool, not a numerical
runner. It requires completed records for all 63 stages, checks the last
attempt of each stage, verifies report/data hashes and checks that the
scripts still match the executed plan. It only creates a fresh staging
tree; it never overwrites the current README, figures or tracked results.

After the numerical campaign completes, in the TJPCov environment:

```bash
python scripts/prepare_publication.py work/global_power_11993 --output work/publication_global_power_11993
```

This command completed on the new results. All 101 staged files were
checked again after copying to the public paths; representative Gaussian,
SSC, halo and all four trispectrum figure layouts were visually inspected.
The publisher retains the README's existing URLs through this mapping:

| New campaign source | Public result/figure destination |
| --- | --- |
| `gaussian_shear_low/high`, `comparison_shear_low/high` | `results/gaussian_low/high`; corresponding Gaussian figure folders |
| `gaussian_3x2_low/high`, `comparison_3x2_low/high` | Existing `results/gaussian_3x2_low/high` and matching figure folders |
| `results/ssc/comparison.json`, `ssc_native.npz` | Existing `results/ssc_native.json` and `results/ssc_native.npz` |
| `ssc_sampling_a2` | New `results/ssc_sampling`; existing `figures/ssc_sampling` |
| `ssc_comparison_low/high_fine` | `results/ssc_low/high_refined`; existing `figures/ssc_low/high_refined` |
| `results/ssc_models` | Existing `results/ssc_models` and `figures/ssc_models` |
| `halo_comparison_i2` and five native exports | Existing `results/halo` and `figures/halo`; regenerated `refinements.json` |
| `results/trispectrum` and power diagnostic | Existing `results/trispectrum` and `figures/trispectrum` |
| `galaxy_bias_placement` | Existing `results/galaxy_bias` |

It also stages the base SSC comparisons, native execution journals, and
`results/refresh_summary.json`. That summary contains the actual Gaussian
component metrics, SSC refinements and model substitutions, halo rows,
native trispectrum metrics and the measured 4h maximum's location/weight.
Native CAMB and installed arrays remain intact in the campaign input
bundle; their configuration and verified hashes are embedded in the
summary, without duplicating dense inputs in every published result.

The entire old tracked result directories should be archived under a
historical work folder before replacing them. Do not merge old halo
subdirectories into the new publication or leave old result JSON beside
new figures. Preserve `results/environment`, which is not recalculated.
No commit or public-file replacement is performed by the publisher.

### README replacements after verified output

- **Overview/physical choices:** identify the installed 11,993-node
  natural-cubic production preparation and original CAMB anchors. Both
  codes receive the same installed tables. Retain native CCL background,
  power interpolation and halo modeling distinctions.
- **Gaussian:** replace both component tables, maximum generalized-mode
  changes, small-entry fractional maximums and source revisions from the
  four new reports. State that the totals are total Gaussian, not G+SSC+cNG.
- **SSC sampling:** replace all percentages using the new sampling report
  and collected refinements. Report CoCoA integration 0/1/2 at fixed global
  boost one. Remove old boost-two results and their reproduction commands:
  boost two now changes the shared power grid and needs its own bundle.
- **SSC models:** replace both native diagonal ranges, all swap/projection
  residuals and local-response examples. Recheck bitwise-reassembly and
  positivity claims against the new reports. The archive now contains 16
  native matrices, not the old 19; link the separate four-way sampling
  report explicitly.
- **Halo:** replace all eight ingredient rows, profile/bias matched-input
  values and integration checks. Keep the comparison at CoCoA integration
  two and CCL's 255-node mass grid, as the new staged plots do.
- **Trispectra:** replace native-order/sum tables, 4h maximum location and
  fractional contribution to the sum, fit-switch examples, common-input
  residuals and refinement values. Remove the old 1,500-to-23,985 table and
  the statement that production is unchanged. The new diagnostic compares
  the installed 11,993 nodes and 23,985-node refinement with the natural
  cubic control, holding moments fixed only in that diagnostic. A current
  native 4h discrepancy must be read from the new native comparison, not
  inferred from the old smooth-control percentages.
- **Galaxy bias:** replace the source revision and measured scaling
  residuals from the refreshed diagnostic; its interpretation remains
  conditional on the actual tested source.
- **Reproduction:** point new commands into one fresh campaign root or use
  the sequential runner's groups. Keep one command per numbered Step box.
  All old `_ab2_i2` and `_checked` private run names should disappear from
  the current recipe. Scientific figure URLs stay stable through the map.

No projected cNG, complete Fourier G/SSC/cNG/total or real-space comparison
becomes complete merely because these 63 ingredient/pilot stages pass.
Those pending scopes must remain explicit. The README does not need the
historical 1,500-node argument or internal campaign execution details.

## Quiet component timings, 2026-10-07

The user requested timing tables and scaling warnings comparable to the
OneCov README. The missing timing phase is now completed for the two
refined SSC pilots and separated matter trispectra. No full TJPCov
G+SSC+cNG or real-space measurement has been added.

`scripts/time_components.py` launches six existing cases three times in
fresh processes, with alternating code order between rounds. It uses
eight OpenMP threads, BLAS one, new output/cache directories and a
600-second hard limit per child. No other numerical/compiler job was
active; the M2 Pro/macOS 13.7.5 preflight found 90.41% idle CPU. The local
execution journal is `work/component_timings_20261007/timing.json`.

All 18 processes passed. `collect_component_timings.py` checked every
scientific setting and native source/binary fingerprint against the
accuracy archives (allowing new commit labels, thread counts and timer
metadata). Every saved scientific array also agrees bitwise in dtype,
shape and bytes. Updated CoCoA core/LSST commit labels are recorded;
the compiled interface fingerprint is unchanged.

| Scope | CoCoA mean ± sample scatter (s) | TJPCov/CCL (s) |
| --- | ---: | ---: |
| Low-ell 5x5 SSC | 2.798584 ± 0.014324 | 5.656082 ± 0.009809 |
| High-ell 5x5 SSC | 2.803663 ± 0.010699 | 5.667027 ± 0.029702 |
| Five matter terms, 45 pairs, three redshifts | 0.159405 ± 0.006791 | 0.164755 ± 0.000545 |

First-use halo tables are included. CoCoA SSC includes geometry, response,
cap variance and projection. TJPCov SSC times its unmodified native block
call, including its small compressed cache write. Benchmark exports and
extra diagnostics are excluded. The original runner's generic exclusion
of file writing was corrected in the publication scope after inspecting
that native cache write; no measurement was altered. The executed runner
is preserved and hash-checked in the work folder. Setup includes CAMB in
CoCoA but saved-table installation in CCL; do not compare those as equal
work. These native models retain their documented physical differences.

`results/component_timings_20261007.json` contains all samples, hashes,
settings and peak memory. The README now gives tables and a warning with
the separately measured full CoCoA LSST real-space time of 53.4 s versus
48.0 s for its source-bin pilot. That 1.11 ratio is a CoCoA reuse example,
not a TJPCov speed ratio. Do not copy OneCov's scaling behavior onto TJPCov.
The trispectrum timing is an ingredient calculation, not projected cNG.

Still pending: fair Gaussian timing (the existing native TJPCov path also
generates spectra while CoCoA receives saved spectra), equivalent halo
ingredient timing scopes, projected cNG and complete Fourier totals.
Native real-space SSC/cNG remain absent from the inspected checkout.
