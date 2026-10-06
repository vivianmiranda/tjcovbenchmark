---
name: tjpcov-comparison
description: Run and document the LSST Y1 covariance comparison between CoCoA and the local TJPCov checkout, including Gaussian, SSC, halo ingredients, separated trispectra, cNG and complete matrix checks.
---

# CoCoA versus TJPCov

Use the same physical cases and reporting contract as the OneCovariance
comparison. Study the local TJPCov checkout; do not replace it with another
clone or infer implemented physics from a paper or class name alone.
Do not optimize or rewrite TJPCov. Commit coherent progress locally;
never push.

Read [the source study](references/tjpcov.md) before choosing a calculator,
input convention or numerical control. Read [the comparison plan](references/comparison_plan.md)
before adding or timing a case. That plan maps the OneCov cases onto actual
TJPCov capabilities and identifies unsupported native comparisons.

## Comparison contract

- LSST Y1 is the baseline. Export the actual CoCoA survey configuration,
  distributions and CAMB tables; do not recreate approximate distributions.
- Use the validated, committed CoCoA Wynn implementation for new results.
  Do not relabel earlier cutoff results as current production results.
- Begin with massless neutrinos, zero IA/RSD/magnification. Check the
  installed CCL angular-spectrum defaults explicitly before calling a
  native TJPCov spectrum Limber or non-Limber.
- Separate shared-input diagnostics from native-model comparisons.
  SACC's mean vector is not a supplied-C_ell input to the f_sky calculator.
- Keep Gaussian signal, mixed noise and pure noise distinct where the API
  permits their extraction. Keep G, SSC, cNG and their total distinct.
- Match ordering, estimator weights, noise units, mass definitions and
  survey windows. Equal bin edges or survey areas do not establish this.
- Ingredient agreement is preliminary. Carry changes through complete
  matrices with every entry retained. Plot all four component differences,
  test total positivity and generalized variance ratios, and report cuts.
  Do not clip eigenvalues or repair a matrix silently. A cNG component need
  not itself be positive definite.
- Use each code's real routines. No independent reference implementation
  in the public comparison. No silent monkeypatches or changes to improve
  agreement. Label any diagnostic adapter and state its scope.
- Do not claim native real-space SSC/cNG: the inspected dispatcher exports
  only RealGaussianFsky. Its EE-only projection also needs a noise audit.

## Installation and runs

Follow Cocoa's schema: Conda base, repository-private .local, sourced
setup/compile/start/stop scripts, choices in set_installation_options.sh.
Keep Cocoa and TJPCov in separate terminals. Follow the Cocoa/LSST README
for Cocoa activation, compilation and thread settings. Never copy a full
Cocoa environment dump; it can contain credentials.

The new TJPCov recipe is not installed or numerically validated yet.
Before publishing results, save resolved versions and a platform-specific
Conda export, check imports and run a minimal source-bin case. Dependency
installation requires authorization; preparing a recipe is not proof that
the environment works.

Run only one numerical job at a time. Use at most eight total workers on
this laptop, BLAS one thread, and derive OpenMP count from OMP_NUM_THREADS.
TJPCov's Wigner multiprocessing pool has a separate worker count; an
OpenMP variable alone does not limit it. Resolve this before a real-space
run. Do not run timing jobs alongside tests, compilation or another code.

Use fresh output directories: TJPCov reuses saved SSC/cNG blocks by name
without a cosmology/settings hash. Never time that disk-cache hit as a
new covariance. Record initialization, halo preparation, projection and
assembly separately. Full-construction timings include first-use tables.
Projection-only timings must be labeled as such.

Start with bounded pilots and record peak memory, elapsed time, revisions,
settings and input hashes. Small output matrices still require converged
physical integration. Never force ell_max=100000 into the native dense
real-space path before assessing its allocation cost. Never lower a cutoff
only to fit memory and then label the result converged.

## Documentation

The README is for current CoCoA-versus-TJPCov results and reproducible
commands. Keep internal studies, historical runs and unadopted proposals
in references. Do not send human readers to Claude skills or private test
folders. Use bullet points for contrasts between the codes; separate short
paragraphs by physical question. Timing differences belong in tables, not
plots. Avoid inline equations in table cells when plain labels suffice.

Use the main Cocoa README style: numbered contents with explicit anchors,
numbered Step blocks, one command per box. Review each completed component
for physics and didactics before proceeding. Keep numerical scripts small,
document units and array axes, and state the scalar meaning of SIMD blocks.
