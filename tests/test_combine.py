"""Tests for multi-instrument, multi-seed CCF map combination."""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pytest
from astropy.io import fits

from exoatmosim.combine import combine_maps


def make_input(root, instrument, seed, species, value, kp=None, snr=None):
    run = root / 'WASP-189_b' / 'simulations' / instrument / f'seed_{seed:04d}'
    out = run / 'cc' / 'sysrem_5' / species
    out.mkdir(parents=True)
    k = np.array([100., 200., 300.]) if kp is None else np.asarray(kp)
    vsys = np.array([-5., 0., 5.])
    for filename, array in [('kp_vsys.fits', np.full((len(k), 3), value)),
                            ('kp_axis.fits', k), ('vsys_axis.fits', vsys)]:
        fits.writeto(out / filename, array)
    if snr is not None:
        fits.writeto(run / 'snr_res_0.fits', np.full((4, 10), snr))
    return run



def fake_poseidon_plot(kp, vsys, *, CCF_Kp_Vsys, **kwargs):
    """Stub visualization: arithmetic tests intentionally use constant maps."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.imshow(CCF_Kp_Vsys, origin='lower',
              extent=[vsys[0], vsys[-1], kp[0], kp[-1]], aspect='auto')
    return fig, ax


def test_average_four_instrument_seed_maps():
    with TemporaryDirectory() as temp:
        root = Path(temp)
        for inst, seed, value in [('VIOLA_J', 1, 1.), ('VIOLA_J', 2, 2.),
                                   ('VIOLA_H', 1, 3.), ('VIOLA_H', 2, 4.)]:
            make_input(root, inst, seed, 'O2', value)
        with patch('exoatmosim.cc.expected_planet_kp', return_value=191.862), \
             patch('POSEIDON.high_res.plot_CCF_Kp_Vsys',
                   side_effect=fake_poseidon_plot):
            result = combine_maps('WASP-189 b',
                                  instruments=('VIOLA_J', 'VIOLA_H'), seeds=(1, 2),
                                  species=('O2',), data_dir=root)
        np.testing.assert_allclose(fits.getdata(result / 'O2' / 'kp_vsys.fits'), 2.5)
        assert (result / 'O2' / 'kp_vsys.png').is_file()


def test_rejects_different_kp_grid():
    with TemporaryDirectory() as temp:
        root = Path(temp)
        make_input(root, 'VIOLA_J', 1, 'O2', 1.)
        make_input(root, 'VIOLA_H', 1, 'O2', 2., kp=[100., 201., 300.])
        with patch('exoatmosim.cc.expected_planet_kp', return_value=191.862):
            with pytest.raises(ValueError, match='Incompatible'):
                combine_maps('WASP-189 b', instruments=('VIOLA_J', 'VIOLA_H'),
                             seeds=(1,), species=('O2',), data_dir=root)


def test_combine_cli():
    from exoatmosim.cli import main
    with patch('exoatmosim.combine.combine_maps', return_value=Path('combined')) as fn:
        assert main(['combine', '--planet', 'WASP-189 b',
                     '--instruments', 'VIOLA_J', 'VIOLA_H', 'VIOLA_K',
                     '--seeds', '1', '2', '3', '4', '--species', 'O2']) == 0
        assert fn.call_args.kwargs['seeds'] == (1, 2, 3, 4)
        assert fn.call_args.kwargs['instruments'] == ('VIOLA_J', 'VIOLA_H', 'VIOLA_K')


def test_snr_squared_weighted_average():
    with TemporaryDirectory() as temp:
        root = Path(temp)
        make_input(root, 'VIOLA_J', 1, 'O2', 10., snr=10.)
        make_input(root, 'VIOLA_H', 1, 'O2', 20., snr=20.)
        with patch('exoatmosim.cc.expected_planet_kp', return_value=191.862), \
             patch('POSEIDON.high_res.plot_CCF_Kp_Vsys',
                   side_effect=fake_poseidon_plot):
            result = combine_maps('WASP-189 b', instruments=('VIOLA_J', 'VIOLA_H'),
                                  seeds=(1,), species=('O2',), data_dir=root)
        np.testing.assert_allclose(fits.getdata(result / 'O2' / 'kp_vsys.fits'), 18.)
        import json
        metadata = json.loads((result / 'O2' / 'combination.json').read_text())
        assert metadata['aggregation'] == 'weighted_mean_of_raw_ccf_maps'
        assert metadata['weights']['VIOLA_J/seed_0001']['weight'] == 100.
        assert metadata['weights']['VIOLA_H/seed_0001']['weight'] == 400.
