#!/bin/bash
# Install the pinned local TJPCov source into .local without downloads.
# TJPCov is Python; Conda already supplies CCL's compiled libraries.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source compile_tjpcov.sh" >&2
  exit 1
fi

(
  cd "$(dirname "${BASH_SOURCE[0]}")" || return 1
  source ./set_installation_options.sh || return 1
  if [[ -z "${CONDA_PREFIX:-}" || ! -x .local/bin/python ||
        -n "${VIRTUAL_ENV:-}" || -n "${ROOTDIR:-}" ]]; then
    echo "Activate only the Conda base and run setup_tjpcov.sh first." >&2
    return 1
  fi
  if [[ "$(git -C "${TJPCOV_PATH}" rev-parse HEAD)" != \
        "${TJPCOV_GIT_COMMIT}" ]]; then
    echo "TJPCov revision differs from set_installation_options.sh." >&2
    return 1
  fi

  .local/bin/python -m pip install "${TJPCOV_PATH}" \
    --no-dependencies --no-index --no-build-isolation || return 1

  # Import each supported calculator without constructing a covariance.
  # This checks the real library dependencies, not just package metadata.
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    .local/bin/python - "${CCL_VERSION}" <<'PY' || return 1
import sys
import pyccl
import sacc
from tjpcov import covariance_from_name

if pyccl.__version__ != sys.argv[1]:
    raise RuntimeError("CCL version differs from installation options")
for name in (
    "FourierGaussianFsky", "FourierSSCHaloModelFsky",
    "FouriercNGHaloModelFsky", "RealGaussianFsky",
):
    print(name, covariance_from_name(name).__module__)
print("CCL:", pyccl.__version__, pyccl.__file__)
print("SACC:", sacc.__file__)
PY
  echo "Installation complete. Next: source start_tjpcov.sh"
) || return 1
return 0
