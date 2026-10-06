# CoCoA vs TJPCov: covariance comparison

This repository compares [CoCoA/CosmoLike](https://github.com/CosmoLike/cocoa)
with [TJPCov](https://github.com/LSSTDESC/TJPCov), the LSST DESC covariance
code. **LSST Y1 is the benchmark survey.**

We follow the component comparisons of
[CoCoA vs OneCovariance](https://github.com/vivianmiranda/OneCov-benchmark-):
Gaussian covariance, super-sample covariance (SSC), halo ingredients,
separate halo trispectra, and connected non-Gaussian covariance (cNG).
Complete matrix comparisons retain all entries and show Gaussian, SSC,
cNG and total differences separately.

**Status:** the local TJPCov source study, installation recipes and first
Gaussian comparison scripts are prepared. Numerical comparisons and
timings have not yet been run.
The new TJPCov environment recipe has not yet been installed and tested.

## Contents

1. [Comparison scope](#scope)
2. [Physical choices](#physics)
3. [Installation and compilation](#installation)
4. [Starting and stopping](#sessions)
5. [CoCoA environment](#cocoa-environment)
6. [First Gaussian comparison](#gaussian)
7. [Results](#results)

## Comparison scope <a name="scope"></a>

The source study uses TJPCov revision
[`2f59302`](https://github.com/LSSTDESC/TJPCov/tree/2f59302af33607aec712185632d8274e59e6b33c).
TJPCov uses CCL for angular spectra and halo-model calculations, and SACC
to describe the survey tracers and the ordering of measured quantities.

| Comparison | TJPCov calculation | Initial scope |
| --- | --- | --- |
| Gaussian Fourier covariance | `FourierGaussianFsky` | One source bin; then two lens bins and one source bin. |
| SSC | `FourierSSCHaloModelFsky` | Background variance, matter response and projected shear covariance. |
| Halo ingredients | CCL routines called by TJPCov | Mass variance, abundance, bias, concentration and halo moments. |
| 1h, 2h, 3h, 4h trispectra | CCL routines called by `FouriercNGHaloModelFsky` | Equal and unequal wavenumbers at several redshifts. |
| cNG | `FouriercNGHaloModelFsky` | Shared-trispectrum projection and native halo-model results. |
| Complete Fourier matrix | Sum of the three Fourier components | One source bin, with all off-diagonal entries. |
| Real-space covariance | `RealGaussianFsky` | One source bin's xi+ and xi−, including their cross-covariance. |

TJPCov's public real-space calculator supplies **Gaussian covariance**.
It currently projects the EE block only. The source study therefore
requires a separate check of shear-noise contributions before a total
real-space Gaussian agreement claim. Native real-space SSC and cNG are
not exposed by the studied checkout.

Mask-coupled Gaussian covariance through NaMaster is a separate TJPCov
capability. The first comparison uses the common sky-fraction
approximation; the installation below does not include NaMaster or MPI.

## Physical choices <a name="physics"></a>

Start with the same LSST Y1 redshift distributions, survey area, number
densities, shape noise and CAMB tables used in the OneCovariance study.
The exporter reads the actual LSST project configuration and preserves
the selected bin's identity and original redshift samples.
The initial comparison uses massless neutrinos and zero intrinsic
alignment, magnification and redshift-space distortions.

Shared inputs isolate an individual calculation. Native-model results
then show the combined effect of each code's physical choices. These
are reported separately.

- **CoCoA:** the galaxy covariance uses linear-bias matter projections.
  The CoCoA baseline selected for this study includes the low-mass Wynn
  treatment in I11;
  higher halo moments use direct integrals.
- **TJPCov:** the connected galaxy covariance uses HOD profiles in the
  one-halo term, and bias-weighted matter terms at higher halo order.
  Cosmic shear avoids this galaxy-model mismatch in the first cNG test.

The native halo prescriptions also differ:

- **CoCoA:** Tinker (2010) multiplicity and bias, with its bias-consistency
  normalization, and the Bhattacharya (2013) concentration prescription.
- **TJPCov:** Tinker (2008) abundance, Tinker (2010) bias, Duffy (2008)
  concentration and NFW profiles, as selected in its SSC and cNG classes.

Equal sky area alone does not match the SSC window:

- **CoCoA:** the survey examples use a spherical-cap footprint.
- **TJPCov:** its sky-fraction SSC calculator calls CCL's circular-disc
  background variance.

These are modeling differences to measure, not evidence of a numerical
error. Numerical refinements must be checked within each code first.

## Installation and compilation <a name="installation"></a>

We assume Conda is installed and `TJPCov/` and `tjcovbenchmark/` are sibling
directories. Open a fresh Bash terminal in `tjcovbenchmark/`. Keep the
CoCoA environment in a separate terminal.

This follows the Cocoa setup/compile/start/stop sequence. Conda supplies
CCL and its compiled dependencies; the repository-private `.local`
environment contains the pinned TJPCov installation.

**Step :one:**: create the separate Conda base.

```bash
conda env create --file=environment.yml
```

**Step :two:**: activate that base.

```bash
conda activate tjcovbenchmark
```

**Step :three:**: inspect
[set_installation_options.sh](set_installation_options.sh) for the local
TJPCov path and revision.

**Step :four:**: create the private environment.

```bash
source setup_tjpcov.sh
```

**Step :five:**: install the local source and check calculator imports.

```bash
source compile_tjpcov.sh
```

**Step :six:**: activate the private environment.

```bash
source start_tjpcov.sh
```

| File | Purpose |
| --- | --- |
| [environment.yml](environment.yml) | Create the separate Conda base, including CCL 3.3.3. |
| [set_installation_options.sh](set_installation_options.sh) | Select the existing TJPCov source path and revision. |
| [setup_tjpcov.sh](setup_tjpcov.sh) | Create `.local` using the active Conda base. |
| [compile_tjpcov.sh](compile_tjpcov.sh) | Install TJPCov offline and check its calculator imports. |
| [start_tjpcov.sh](start_tjpcov.sh) | Activate `.local`. |
| [stop_tjpcov.sh](stop_tjpcov.sh) | Restore the shell's previous Python environment. |

TJPCov itself is Python; the compile step installs it without rebuilding
CCL or changing Cocoa. The recipe pins Python, CCL and the principal
numerical packages. It is not yet a platform lock file: the remaining
dependency versions will be recorded after the first successful install.

## Starting and stopping <a name="sessions"></a>

After installation, open a fresh Bash terminal in `tjcovbenchmark/`.
Setup and compilation are not needed for each calculation.

**Step :one:**: activate the Conda base.

```bash
conda activate tjcovbenchmark
```

**Step :two:**: activate the private environment.

```bash
source start_tjpcov.sh
```

**Step :three:**: after the calculation, leave the private environment.

```bash
source stop_tjpcov.sh
```

## CoCoA environment <a name="cocoa-environment"></a>

Use the [main Cocoa installation recipe](https://github.com/CosmoLike/cocoa#required_packages_conda)
and the [LSST Y1 covariance instructions](https://github.com/CosmoLike/cocoa_lsst_y1#computing_covariances).
The LSST Y1 interface must be compiled with covariance support enabled.

We assume Cocoa is already installed. Open a **separate** Bash terminal
in its `Cocoa/` directory. Replace `cocoa` below if the installed Conda
base has a different name.

**Step :one:**: activate Cocoa's Conda base.

```bash
conda activate cocoa
```

**Step :two:**: load Cocoa's runtime environment.

```bash
source start_cocoa.sh
```

**Step :three:**: select the benchmark's OpenMP thread count.

```bash
export OMP_NUM_THREADS=8
```

**Step :four:**: after the calculation, leave Cocoa's runtime environment.

```bash
source stop_cocoa.sh
```

## First Gaussian comparison <a name="gaussian"></a>

The prepared first case uses **LSST Y1 source bin 3**, five Fourier bands
over `30 <= ell < 150`, and its complete 5 × 5 Gaussian covariance.
The source density and ellipticity noise come from the project settings.

This first test asks whether the two codes assemble the same covariance
from the same angular spectra:

- **TJPCov:** reads the LSST distribution through SACC and imports CoCoA's
  CAMB power tables through CCL. Its native Gaussian calculator computes
  the shear spectrum and covariance. The runner exports those same CCL
  spectrum samples and the bin operator used by TJPCov.
- **CoCoA:** receives the exported spectra, noise and operator through
  its production covariance interface. Its C kernels compute the matrix.

TJPCov's bin weights are proportional to `ell` on this integer grid.
The scripts check its reconstructed edges and match those weights in
CoCoA. The final endpoint is present to define the last edge and has zero
weight in the covariance. SACC's stored mean values are not used as
supplied spectra by this TJPCov calculator.

Sample variance, signal–noise and pure-noise pieces are saved separately.
TJPCov does not expose that split directly: three public calls with noise
power multiplied by zero, one and two recover the coefficients of its
quadratic dependence on noise. CoCoA obtains the same pieces by turning
the supplied signal or noise off.

This isolates **Gaussian assembly**. It does not compare native spectra,
validate the halo model, or generate SSC/cNG. CCL uses the supplied CAMB
power but its own background and interpolation; a later native-spectrum
comparison must measure those differences. The three noise calls also
do not constitute a production timing measurement.

### Run the prepared case

These commands assume the sibling directory layout described above.
Finish each step before starting the next; use the two activated terminals
from the environment instructions. Use new output names for a new case.

**Step :one:**: in the **CoCoA terminal**, from `cocoa/Cocoa/`, export the
actual survey inputs and CAMB tables. This does not compute a covariance.

```bash
python ../../tjcovbenchmark/scripts/export_lsst_y1.py \
  --cocoa . --output ../../tjcovbenchmark/work/lsst_y1
```

**Step :two:**: in the **TJPCov terminal**, from `tjcovbenchmark/`, select
the OpenMP allocation.

```bash
export OMP_NUM_THREADS=8
```

**Step :three:**: keep the numerical-array library single-threaded in
that terminal.

```bash
export OPENBLAS_NUM_THREADS=1
```

**Step :four:**: compute TJPCov's native Gaussian components and export
the angular spectra it uses.

```bash
python scripts/run_gaussian.py work/lsst_y1 --tjpcov ../TJPCov \
  --output work/gaussian_low
```

**Step :five:**: return to the **CoCoA terminal**, still in `cocoa/Cocoa/`,
and compare its Gaussian matrix using those exact arrays.

```bash
python ../../tjcovbenchmark/scripts/compare_gaussian.py \
  ../../tjcovbenchmark/work/gaussian_low --cocoa . \
  --output ../../tjcovbenchmark/work/comparison_low
```

**Step :six:**: in the **TJPCov terminal**, plot the passing comparison.

```bash
python scripts/plot_gaussian.py work/comparison_low \
  --output figures/gaussian_low
```

The matrix figure shows each code's correlations and the signed
TJPCov-minus-CoCoA difference, normalized by CoCoA's total rms product.
A second figure separates the three Gaussian variance contributions.
Every matrix entry is retained; the report also checks positivity and
generalized variance ratios. No eigenvalues are clipped.

For the matching high-multipole case, use `--ell-range 1500 1620` in
Step 4 and fresh output paths in Steps 4–6. Reuse the exported LSST bundle;
the runners verify its hashes before using it.

| Script | Purpose |
| --- | --- |
| [export_lsst_y1.py](scripts/export_lsst_y1.py) | Save the project's actual catalog, noise, cosmology and CAMB arrays. |
| [run_gaussian.py](scripts/run_gaussian.py) | Run native TJPCov shear Gaussian blocks and export the shared spectra/operator. |
| [compare_gaussian.py](scripts/compare_gaussian.py) | Call CoCoA's production C kernels and compare all Gaussian entries. |
| [plot_gaussian.py](scripts/plot_gaussian.py) | Make correlation, residual and variance-component panels in the OneCov comparison style. |

Manifests record the input hashes, survey identities, units, code
revisions, installed TJPCov source hashes and package versions. These
scripts have passed syntax and source/API review; their first numerical
execution still awaits installation of the separate TJPCov environment.

## Results <a name="results"></a>

The first numerical result will compare the Gaussian Fourier covariance
for LSST Y1 source bin 3. Subsequent results will include component
difference plots and execution-time tables, with the input conventions,
tested revisions and numerical refinements stated beside them.

No TJPCov runtime or cross-code accuracy has been measured yet.
