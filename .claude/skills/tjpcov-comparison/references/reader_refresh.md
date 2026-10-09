# Log-domain reader refresh of the complete pilot — 2026-10-09

CoCoA adopted a log-domain linear-power reader for the connected
covariance, a per-call z slice of the bilinear read, and a blocked
fusion feeding the tree averages directly from that reader. The merged
core is `23b0127`; interface SHA-256
`3c361f80d018377cda141fd710e75adf27befe89ab945d371bd23c26cb27f1fb`.

The complete-pilot CoCoA cases were rerun against the unchanged
archived input bundle (`work/global_power_11993/lsst_y1`, manifest
`50ec3809...`) at six OpenMP threads: `cocoa_i0` and `cocoa_i1` under
`work/reader_refresh_20261009/`. Elapsed accuracy-run walls dropped
from about 31 s to 19 s (i0) — these are not quiet timings and must
not be published as such. Against the archived 2026-10-07 runs,
Gaussian, SSC, ell, edges, both signals, spectra, noise power and
operators are bitwise equal; cNG moves at most 2.3e-36 in the saved
units (entries are of order 1e-20). The collector was rerun with the
new CoCoA cases and the unchanged archived natives
(`native_k32_a2`, `native_k96_a2`, refinements k16/k32/k64); every
published diagnostic in `results/complete_fourier/report.json` is
identical at displayed precision, and the record now carries core
`23b0127` provenance. The collector requires `OPENBLAS_NUM_THREADS`,
`MKL_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` all set to 1.

The shared-input cNG projection, time-grid control and assembly
verification records were not rerun: their C kernels are untouched by
the reader work, and their records accurately state the core at their
measurement time. The component timing tables (2026-10-07) time SSC
and halo-trispectrum kernels the reader work does not touch; they are
retained as dated measurements. The cross-repo full-LSST claim in the
README now quotes the OneCov-benchmark reader-refresh record
(29.5 s full versus 23.5 s real pilot, 1.25x), measured 2026-10-08 on
the same merged core. Quiet TJPCov complete-pilot timings remain
pending, unchanged by this refresh.
