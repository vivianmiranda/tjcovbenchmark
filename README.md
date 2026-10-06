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

**Status:** the local TJPCov source study and installation recipes are
prepared. Numerical comparisons and timings have not yet been run.
The new TJPCov environment recipe has not yet been installed and tested.

## Contents

1. [Comparison scope](#scope)
2. [Physical choices](#physics)
3. [Installation and compilation](#installation)
4. [Starting and stopping](#sessions)
5. [CoCoA environment](#cocoa-environment)
6. [Results](#results)

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

## Results <a name="results"></a>

The first numerical result will compare the Gaussian Fourier covariance
for LSST Y1 source bin 3. Subsequent results will include component
difference plots and execution-time tables, with the input conventions,
tested revisions and numerical refinements stated beside them.

No TJPCov runtime or cross-code accuracy has been measured yet.
