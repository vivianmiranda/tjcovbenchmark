"""Publish the saved native cNG time-grid control; no covariance is rerun."""
import argparse
import shutil
from pathlib import Path

import numpy as np
from common import load_bundle, sha256, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coarse', type=Path, required=True)
    p.add_argument('--fine', type=Path, required=True)
    p.add_argument('--assembly', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    records = [load_bundle(x, 'tjpcov-complete-shear-v1') for x in (a.coarse, a.fine)]
    assert all(r['status'] == 'completed' for r in records)
    assert records[0]['input_manifest_sha256'] == records[1]['input_manifest_sha256']
    settings = [r['components']['cng'] for r in records]
    for key in ('N_K', 'k_count', 'k_range_Mpc_inverse', 'integration_method'):
        assert settings[0][key] == settings[1][key]
    with np.load(a.coarse/'covariance.npz') as c, np.load(a.fine/'covariance.npz') as f:
        np.testing.assert_array_equal(c['ell'], f['ell'])
        coarse, fine, ell = c['cng'], f['cng'], f['ell']
    assert coarse.shape == fine.shape == (100, 100)
    assert np.isfinite(coarse).all() and np.isfinite(fine).all()
    assert np.all(np.diag(fine) > 0)
    scale = np.sqrt(np.diag(fine)[:, None]*np.diag(fine)[None, :])
    residual = (coarse-fine)/scale
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'time_grid.npz', coarse=coarse, fine=fine,
                        ell=ell, scaled_difference=residual)
    shutil.copy2(a.assembly, a.output/'assembly_verification.json')
    write_json(a.output/'manifest.json', {
        'schema': 'cng-time-grid-publication-v1', 'status': 'completed',
        'scope': 'Native TJPCov cNG at fixed 234-node wavenumber grid; every entry retained',
        'maximum_difference_percent_of_fine_cng_rms': float(100*np.max(abs(residual))),
        'normalization': 'sqrt(fine cNG diagonal_i times fine cNG diagonal_j)',
        'source_manifests': records,
        'source_manifest_sha256': [sha256(x/'manifest.json') for x in (a.coarse, a.fine)],
        'script_sha256': sha256(__file__),
        'files': {n: sha256(a.output/n) for n in ('time_grid.npz', 'assembly_verification.json')},
    })


if __name__ == '__main__':
    main()
