# Native SSC pilot preparation — 2026-10-06

The user authorized small TJPCov runs on up to six OpenMP threads while
the other correctness checks finish. Timing comes later on a quiet
machine. The coordinator runs the jobs; this preparation was source-only.

## Prepared entry point

`scripts/run_ssc.py` reuses the completed Gaussian pilot's SACC file and
the exact LSST/CAMB bundle identified in its manifest. It calls the native
`FourierSSCHaloModelFsky.get_covariance_block`, with all four legs equal
to the selected source bin and `include_b_modes=False`. No library method
is replaced. A new output directory prevents a filename-only SSC cache
hit. It saves the complete 5x5 SSC matrix, its native multipoles, the disc
variance used by the native window prescription and source/input hashes.

Use the active TJPCov environment from the README. The first low-ell run:

```bash
python scripts/run_ssc.py work/lsst_y1 work/gaussian_low --tjpcov ../TJPCov --match-power-a-range --output work/ssc_low_qag
```

The independently selected native integration algorithm:

```bash
python scripts/run_ssc.py work/lsst_y1 work/gaussian_low --tjpcov ../TJPCov --match-power-a-range --integration-method spline --output work/ssc_low_spline
```

The combined table-density refinement, after separate effects are checked:

```bash
python scripts/run_ssc.py work/lsst_y1 work/gaussian_low --tjpcov ../TJPCov --match-power-a-range --a-refinement 2 --k-refinement 2 --output work/ssc_low_refined
```

Repeat with `work/gaussian_high` and fresh outputs only after reviewing
the first case. These commands were prepared, not executed by this agent.
They generate native SSC components, not CoCoA-versus-TJPCov comparisons.
Elapsed block times are diagnostic and are not quiet-machine benchmarks.
The script's `--timeout` defaults to 600 seconds. A Unix SIGALRM uses
the OS default action so it can stop compiled code as well as Python.
The supervising runner must preserve its log and recognize that signal;
run_settings.json and any native outputs are retained after interruption.

## CAMB support and the native response grid

The Gaussian run works with the shared CAMB power table covering
0 <= z <= 49.99. Its minimum scale factor is about 0.0196117. The runtime
CCL default scale-factor grid begins at 0.01 (z=99).

TJPCov truncates the sigma2_B sampling to the tracer redshift range, but
it does not pass that smaller grid to halomod_Tk3D_SSC_linear_bias. That
routine builds the halo response at every default CCL scale-factor node
and sets `extrap_pk=False`. Thus the untreated shared-table calculation
would evaluate power outside its supplied time support.

The explicit `--match-power-a-range` flag uses CCL's public control
`spline_params.A_SPLINE_MINLOG_PK`, before CosmologyCalculator creation.
It sets the response-grid endpoint to the actual supplied CAMB endpoint.
It retains A_SPLINE_MIN_PK, A_SPLINE_NLOG_PK and A_SPLINE_NA_PK, and saves
the complete resulting grid. It does not alter or monkeypatch TJPCov,
extrapolate CAMB power in redshift, or change the projected source range.
Without this flag, a guard explains the unsupported time range and stops.

This is a disclosed domain adaptation, not a convergence result. Check
that the source kernel's complete support lies inside the retained range.
Changing this control only redistributes the low-a logarithmic segment;
the higher-a linear segment stays unchanged. For a separate refinement
study, compare denser response grids and/or CAMB inputs explicitly extended
to z >= 99. Do not manufacture high-redshift tables by rescaling P with an
assumed growth factor and call them shared CAMB inputs.

## Supported native controls

- TJPCov exposes `SSC.integration_method`: `qag_quad` (default) or `spline`.
  The argument to get_covariance_block overrides the configuration entry.
- It hardcodes M200m, Tinker08 abundance, Tinker10 bias, Duffy08
  concentration, analytic NFW and a CCL HMCalculator with its defaults.
  In the inspected CCL 3.3.3 source these are log10M = 8..16, nM=128,
  Simpson integration. Mass is in solar masses, not solar masses/h.
- It does not expose hmc mass bounds/counts or response k/a arrays in its
  public configuration. CCL's documented global sampling controls can
  define a separately recorded refinement before cosmology construction;
  replacing the internal HMCalculator is a different diagnostic adapter.
- No MPI is requested. OpenMP comes only from OMP_NUM_THREADS; this pilot
  rejects more than six threads. BLAS remains single-threaded.
- Native blocks are cached under tracer names with no input hash. Never
  reuse their output directory for a new cosmology/control or timing run.

### Public sampling refinements to test next

CCL exposes `spline_params.A_SPLINE_NA_PK` and
`spline_params.A_SPLINE_NLOG_PK` for the linear and logarithmic pieces of
the scale-factor grid. Set them before constructing CosmologyCalculator.
The pieces share the transition point A_SPLINE_MIN_PK. For a refinement
factor r, use `r*(n-1)+1` for each count rather than `r*n`; this preserves
the old segment nodes. Verify that fact on the actual returned native
grid before claiming nested sampling. Save both response-a arrays.
The implemented `--a-refinement 2` option follows this interval rule and
asserts that taking every second refined node recovers the original grid
within 5e-14 relative rounding tolerance. This native control refines both
the response table and the a sampling of sigma2_B, because TJPCov derives
both from CCL's power-grid sampling.

The public wavenumber-density control is `spline_params.N_K` (points per
decade); bounds are `K_MIN` and `K_MAX`. CCL's pyutils reports
`nk = ceil(log10(K_MAX/K_MIN)*N_K)`. Doubling N_K generally does **not**
preserve old nodes when the endpoints are included: nested refinement
would need `nk_new = r*(nk_old-1)+1`, which an integer density parameter
may not represent. Record the actual old/new log-k arrays and distinguish
density refinement from nested refinement. The native TJPCov class does
not forward an explicit lk_arr to its CCL response helper, so forcing an
arbitrary nested k grid would require a separately labeled lower-level
CCL diagnostic, not an undisclosed native-run modification.
The implemented `--k-refinement 2` doubles N_K, leaves the bounds fixed
and saves the original and refined log-k arrays. No nested-k claim is made.

Keep K_MIN/K_MAX fixed while testing density. A cutoff change is a
different test. Neither density refinement establishes the accuracy of
CCL's fixed HMCalculator mass grid (128 Simpson nodes in this checkout).
Likewise, qag_quad versus spline tests the radial integration route, not
the response-table sampling. Test these controls one at a time on the
entire 5x5 matrix, then a combined refinement if each effect is understood.

## Scientific differences to resolve before a total comparison

1. **Estimator.** Gaussian band averages and SSC point samples are not
   the same observable estimator. The SSC uses SACC mean/effective ell
   directly and never calls bin_cov. CoCoA must first be evaluated at the
   same ell for a point-sampled SSC comparison. A common band-averaged
   diagnostic is a separate test. Do not silently add this native SSC
   matrix to the matched band-averaged Gaussian and label a common total.

2. **Window.** CCL integrates the transverse linear-power spectrum with
   a circular-disc Bessel window. CoCoA's forecast uses a spherical-cap
   harmonic mask. The same area does not equate those integrations.
   CCL uses R=chi*acos(1-2*fsky). Its sigma2_B has units of Mpc because
   k dk times the three-dimensional power has length units; it is not a
   dimensionless three-dimensional cell variance. Respect that unit when
   transferring it into CoCoA's internal length convention.

3. **Matter response.** TJPCov's selected CCL helper uses
   (47/21 - dlnP_linear/dlnk / 3)*P_linear + I12/normalization_squared.
   It does not include I11-squared in this two-halo response. CoCoA uses
   the slope of its two-halo power and transfers the fractional halo
   response to its nonlinear power. Shared CAMB tables do not equate
   these choices. In the one-source case no number-count correction
   applies. Galaxy local-mean subtraction requires its own later audit.

4. **Halo fits and finite integration.** The native Tinker08/Duffy08
   choices differ from CoCoA's Tinker10/Bhattacharya13 model and its
   extended-mass Wynn treatment of I11. Compare ingredients before
   interpreting their native SSC difference as a numerical error.

5. **Geometry and tracer interpolation.** CCL keeps its own background,
   growth and n(z) interpolation even when given the common CAMB P.
   The shear-spin factors are handled by the native tracer projection;
   do not multiply them again when exporting a projected response.

Next, separate native window, matter response and projection comparisons.
The prepared script has no cross-code convergence claim and does not
insert uncomputed SSC/cNG arrays as zeros.
