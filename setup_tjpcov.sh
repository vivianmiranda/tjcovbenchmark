#!/bin/bash
# Source after creating and activating the benchmark's Conda environment.
# Conda supplies dependencies; this creates the repository-private layer.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source setup_tjpcov.sh" >&2
  exit 1
fi

(
  cd "$(dirname "${BASH_SOURCE[0]}")" || return 1
  source ./set_installation_options.sh || return 1

  if [[ -z "${CONDA_PREFIX:-}" || -n "${VIRTUAL_ENV:-}" ||
        -n "${ROOTDIR:-}" ]]; then
    echo "Activate only the tjcovbenchmark Conda base in a fresh shell." >&2
    return 1
  fi

  PYTHON="${CONDA_PREFIX}/bin/python"
  if [[ "$("${PYTHON}" -c \
    'import sys; print("%d.%d" % sys.version_info[:2])')" != \
    "${PYTHON_VERSION}" ]]; then
    echo "The Conda base must use Python ${PYTHON_VERSION}." >&2
    return 1
  fi
  if [[ "$(git -C "${TJPCOV_PATH}" rev-parse HEAD)" != \
        "${TJPCOV_GIT_COMMIT}" ]]; then
    echo "TJPCov revision differs from set_installation_options.sh." >&2
    return 1
  fi

  "${PYTHON}" -m venv --system-site-packages .local || return 1
  echo "Setup complete. Next: source compile_tjpcov.sh"
) || return 1
return 0
