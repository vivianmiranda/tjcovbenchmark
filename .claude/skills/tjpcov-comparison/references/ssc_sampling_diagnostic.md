# Separating the SSC scale-factor tables

Executed on 2026-10-06. The coordinator ran all four cases; the numerical
findings and limits are recorded below.

The native TJPCov SSC calculator uses CCL's power-grid scale factors for
two different tables. The halo response is tabulated across the complete
scale-factor range. The transverse survey-window variance uses the subset
covering the source redshift support, including one bracketing point.
Increasing the native a-grid density therefore changes both quantities.

`scripts/diagnose_ssc_sampling.py` separates them through public CCL APIs.
It reads a coarse native run and an a-only refined native run, verifies
their manifests and files, and reuses their actual saved a arrays. The
coarse cosmology controls, shared CAMB power, SACC tracer, effective ell,
response k grid and radial integration algorithm stay fixed.

| Case | Halo-response a sampling | Disc-variance a sampling |
| --- | --- | --- |
| coarse_coarse | Coarse | Coarse |
| fine_coarse | Fine | Coarse |
| coarse_fine | Coarse | Fine |
| fine_fine | Fine | Fine |

The response is produced by `halomod_Tk3D_SSC_linear_bias`, with explicit
`a_arr` and `lk_arr`. The halo choices exactly follow TJPCov: M200m,
Tinker08 abundance, Tinker10 bias, Duffy08 concentration, analytic NFW and
the default HMCalculator. All four legs are shear, so the default unit
biases and absence of number-count subtraction are appropriate.

The variance is produced by `sigma2_B_disc`, with explicit `a_arr` and
the unchanged fsky. Each covariance comes from `angular_cl_cov_SSC` with
the actual TJPCov tracer from `get_tracer_info()`. No response, variance,
spin correction or projection formula is independently reimplemented.

The coarse/coarse covariance must first reproduce every entry of the
saved native coarse covariance to 1e-8 when the absolute difference is
divided by its diagonal rms product. This is an adapter-reproduction
gate, not a scientific accuracy target. The script saves that matrix and
the measured discrepancy before stopping on a failed gate. Do not infer
a cause from the remaining cases if the baseline gate fails.

For an authorized run, use the TJPCov environment from the README:

```bash
python scripts/diagnose_ssc_sampling.py work/lsst_y1 work/gaussian_low work/ssc_low_qag work/ssc_low_a2 --tjpcov ../TJPCov --output work/ssc_sampling_a2
```

The run requires OMP_NUM_THREADS between one and six and single-threaded
BLAS. The default wall-time deadline is 600 seconds; `--timeout` sets a
different positive deadline. A fresh output directory is required. An OS
SIGALRM interrupts compiled calls as well as Python, preserving partial
`report.json` and `sampling.npz` outputs for the supervising runner.

The output retains all four complete 5x5 SSC matrices, both sets of a
nodes, the common log-k nodes and both variance arrays. Reported case
durations include differing amounts of reused work and are diagnostic;
they are not fair speed comparisons. Native source and binary hashes are
checked against the input records before calculation.

CAMB k and P are converted from h/Mpc and (Mpc/h)^3 to CCL's Mpc^-1 and
Mpc^3. The disc variance has units Mpc, because its transverse integral
contains k dk times a three-dimensional power spectrum. CCL applies all
distance and spin factors during projection. No additional area, chi or
shear-spin factors belong in this diagnostic.

The first coarse native record predates saved response log-k arrays. In
that case the diagnostic obtains unchanged defaults from the same CCL
build, checks them against the refined run's saved k grid and additionally
requires full coarse-matrix reproduction. A k-refined input is rejected.

Compare fine/coarse and coarse/fine against coarse/coarse to locate the
effect. Their sum need not equal the combined change because projection
multiplies the response and variance. The report also compares fine/fine
with the saved native fine result: any residual there signals an effect
of changing CCL's internal cosmology tables or an incomplete recreation.
This test isolates two interpolation samplings. It does not establish
mass integration, physical response modeling, cross-code SSC agreement,
or a converged full covariance.

## Measured result

The coarse/coarse and fine/fine diagnostics reproduce the corresponding
native TJPCov matrices exactly. Relative to the coarse SSC diagonal rms
product, refining only the response a grid changes entries by 0.00311%;
refining only the disc-variance a grid changes them by 23.9959%. Refining
both gives 23.9963%. Doubling k sampling gives only 0.0000735%.

Native a-grid refinement retains shared nodes. Using the finer matrix's
SSC rms normalization, 50->99->197->393->785 response-grid nodes gives
successive changes 31.5726%, 2.6184%, 0.15742%, 0.0179984% at low ell.
The last change at high ell is 0.000496484%. At 393 nodes, qag_quad versus
spline changes low-ell SSC by 0.000509868%. The first 31.57% uses a different
normalization from the 24% diagnostic; this is not a numerical conflict.

CoCoA integration 0->1->2 changes the low-ell SSC by 0.0565306% and
0.0457384%, respectively. At integration 2, global boost 1->2 changes
entries by 0.00588654%. All runs use six or fewer OpenMP workers and one
BLAS worker. Their elapsed times are diagnostic only because project
regression validation was active concurrently.

The refined native comparison retains different halo fits, responses,
background/tracer interpolation and cap/disc windows. TJPCov shear SSC
diagonals are 10.48--26.62% below CoCoA at low ell, 20.34--20.55% below at
high ell. Do not attribute the whole residual to a single ingredient.
The high-ell component is almost rank deficient: the initial small
negative mode disappears after refinement, with no eigenvalue clipping.
SSC-only weak-mode ratios do not certify total covariance accuracy.

The native archive includes all 19 cases, all four diagnostics and raw
matrices. It is not a uniformly band-averaged total: native TJPCov SSC
samples band centers, whereas its Gaussian calculator averages bands.
