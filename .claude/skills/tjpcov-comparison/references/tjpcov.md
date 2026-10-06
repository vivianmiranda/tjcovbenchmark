# Local source study — 2026-10-06

Inspected TJPCov commit `2f59302af33607aec712185632d8274e59e6b33c`.
All paths below are relative to the sibling TJPCov checkout. During this initial source review, no source changes, imports, covariance
runs or TJPCov dependency installs were made. The later numerical campaign
is recorded in comparison_plan.md and ssc_sampling_diagnostic.md.
Also inspected the already installed CCL 3.3.3 Python source, read-only.
The proposed environment selects that release, not the data-vector study's
CCL PR #1296. Record the actual installed code before each campaign.

## Entry points and data conventions

- `run_tjpcov.py` uses CovarianceCalculator to write a SACC covariance.
  `tjpcov/covariance_calculator.py` exposes individual terms through
  get_covariance_terms(); create_sacc_cov saves separate terms by default.
- `tjpcov/covariance_io.py` accepts a YAML path or a dict. sacc_file may be
  a FITS path or a Sacc object. SACC supplies tracers, binning and ordering;
  its mean C_ell values do not replace spectra calculated by the f_sky API.
- `covariance_builder.py::get_cosmology` accepts a CCL Cosmology object.
  Check CosmologyCalculator compatibility in the installed dependency
  before using this route for shared CAMB tables. It also supports YAML
  parameters or a saved cosmology. No benchmark-specific path exists yet.
- get_tracer_info converts Ngal from arcmin^-2 to sr^-1. Shear noise is
  sigma_e^2/Ngal, with sigma_e per ellipticity component; galaxy noise is
  1/Ngal. Tracer names and quantity metadata determine source/lens type.
- Lens tracers use has_rsd=False, constant input bias and no magnification.
  IA=None omits IA. A supplied IA value feeds a constant amplitude array
  into CCL WeakLensingTracer; there is no TATT switch in this builder.
- Gaussian calls omit angular_cl's Limber options. In inspected CCL 3.3.3,
  cells.py defaults l_limber=-1: all nonnegative multipoles use Limber.
  Recheck if the dependency changes; TJPCov itself does not set this.

## Gaussian Fourier

`covariance_gaussian_fsky.py::FourierGaussianFsky` computes four CCL spectra
for the Wick contractions 13×24 and 14×23. Noise is added only for identical
tracers. It forms a diagonal ell covariance, divides by
(2ell+1) gradient(ell) fsky and then calls wigner_transform.py::bin_cov.

This binning routine weights by ell*dell; in real space it weights by
theta*dtheta. Neither equals CoCoA's (2ell+1) band weighting or sin(theta)
annular weighting automatically. Supply explicit SACC windows, inspect
the selected discrete ell nodes and export the actual weights. A common
operator diagnostic can use CoCoA's supplied-operator interface without
rewriting either code's contractions. Check bin endpoints explicitly.

The native f_sky function does not expose supplied C_ell arrays. A possible
shared-spectrum test exports the four spectra used by TJPCov and feeds
them to CoCoA. Do not monkeypatch ccl.angular_cl or claim that SACC means
are consumed. Report native TJPCov spectrum generation separately from
assembly; its block-call time includes both unless measured otherwise.

`covariance_fourier_gaussian_nmt.py` is a separate NaMaster mask-coupled
path. Its cache accepts cl13/cl24/cl14/cl23, noise arrays, masks, fields,
workspaces and bins. It uses the narrow-kernel approximation. Optional
pymaster>=2 is required. It is not the initial common-f_sky comparison.

## SSC

`covariance_fourier_ssc.py` fixes M200m, MassFuncTinker08, HaloBiasTinker10,
ConcentrationDuffy08, analytic NFW and CCL HMCalculator. It uses
halomod_Tk3D_SSC_linear_bias with the four biases and number-count flags.
The scale-factor table follows CCL's P(k) nodes over the tracer overlap.
Projection uses angular_cl_cov_SSC at effective multipoles, not the
Gaussian calculator's explicit band average. Its integration choices are
qag_quad (default) or spline.

- `covariance_fourier_ssc_fsky.py` calls sigma2_B_disc with fsky.
- The mask path computes the power of mask products m12 and m34, with
  (2ell+1)/(integral(m12) integral(m34)) normalization, then calls
  sigma2_B_from_mask. Equal area is insufficient to match these windows.

CCL 3.3.3 halos/pk_4pt.py::halomod_Tk3D_SSC_linear_bias combines
(47/21 - dlnP/dlnk/3)*P and I12, and applies number-count corrections when
requested. Check the actual selected P, profile normalizations, density
reference and bias placement before comparing with CoCoA's two-halo
slope and fractional-response transfer. The later native SSC comparisons
are recorded in ssc_sampling_diagnostic.md.

SSC blocks are reused from ssc_tr1_tr2_tr3_tr4.npz if present. Filenames
do not identify cosmology or accuracy settings. Always use a fresh outdir
for a new configuration or a full-construction timing.

## Connected non-Gaussian covariance

`covariance_fourier_cNG.py` uses the same M200m/Tinker08/Tinker10/Duffy08/NFW
choices. Its CCL calls expose 2h22, 2h13, 3h and 4h independently; the class
sums them. separable_growth=True is passed to 2h22, 3h and 4h. It multiplies
this higher-halo sum by the product of four configured linear biases.

The one-halo term instead selects HOD profiles for number-count legs and
NFW profiles for shear legs. HOD-HOD pairs use Profile2ptHOD; other pairs
use Profile2pt. Thus galaxy cNG is not CoCoA's bias-weighted matter model.
The constructor requires the complete HOD configuration even for shear;
copy all keys from its current example, do not invent them or count on
defaults. For shear, verify that the chosen HOD object does not enter the
selected profiles.

It packages the trispectrum in a Tk3D with is_logt=False and calls
angular_cl_cov_cNG at effective multipoles. Unlike the Gaussian path,
it does not apply bin_cov. Record this estimator difference and perform
a separate common-binning diagnostic before interpreting discrepancies.

Use the same public CCL functions for the separated halo-order comparison;
do not recode their algebra. Check both 1+3 permutations and K/Q symmetry.
The OneCov off-diagonal discrepancy is not evidence that CCL has the same
problem. Inspect CCL's own equal-pair limits and angular integration.
Saved cng_*.npz blocks have the same filename-only reuse as SSC.

## Real space: capability and resource gates

`tjpcov/__init__.py` exposes RealGaussianFsky only. That class in
covariance_gaussian_fsky.py calls FourierGaussianFsky with for_real=True;
the source explicitly selects the EE block. Audit BB pure-noise effects
in xi+/xi− and their cross-covariance before claiming a complete Gaussian
noise treatment. There are no exported native real-space SSC/cNG classes.

`covariance_builder.py::CovarianceProjectedReal` requires ProjectedReal.lmax.
An example's lmax=90 is a test setting, not a converged survey default.
The Fourier path allocates np.diag of an lmax+1 vector: one float64 matrix
at lmax=100000 requires about 80 GB, before intermediates. Do not launch
that case on this laptop or lower lmax and call the result converged.

The real-space builder assumes logarithmic angular bins when reconstructing
edges from centers. It samples theta with separate log and linear pieces.
The transform uses integer ell from 2 to lmax. wigner_transform.py uses
Jacobi-polynomial Wigner d below ell=10000 and a Bessel approximation above
that value. Do not label TJPCov as purely flat-sky or exact full-sky over
all ell. CoCoA uses full-sky annular transforms.

WignerTransform has an ncpu argument, but the high-level builder does not
pass it. wigner_d_parallel defaults to multiprocessing.cpu_count(), which
is independent of OMP_NUM_THREADS. Plan an explicit, documented use of the
lower-level API or a supported worker limit before running here. Do not
silently monkeypatch cpu_count and report an untouched native run.

## Environment status

pyproject.toml requires Python>=3.10, pyccl>=3.2.0, sacc>=0.12, SciPy,
NumPy, YAML, Jinja2, CAMB, healpy and h5py. It uses setuptools_scm at build
time. The README's older Python text is not the packaging contract.
NaMaster and MPI are optional, not part of the first f_sky environment.

The benchmark recipe pins CCL 3.3.3 and a separate Python 3.12 base. It
does not reuse Cocoa's Python 3.11 environment or the data-vector study's
CCL PR. Installation, imports and pip dependency checks subsequently
passed. Resolved package records are in results/environment/.
