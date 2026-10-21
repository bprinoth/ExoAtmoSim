"""Smoke tests for native CC visualizations (no POSEIDON dependency)."""
import numpy as np
import pytest
from exoatmosim.cc_plots import plot_kp_vsys, plot_phase_rv, standardize_ccf


def test_native_cc_figures(tmp_path):
    kp = np.arange(-100., 301., 2.)
    vsys = np.arange(-200., 201., 2.)
    kgrid, vgrid = np.meshgrid(kp, vsys, indexing='ij')
    raw = .15 * np.sin(vgrid / 12. + kgrid / 17.)
    raw += 5 * np.exp(-((kgrid - 190.) / 8.)**2 - (vgrid / 5.)**2)
    image = tmp_path / 'kp_vsys.png'
    plot_kp_vsys(kp, vsys, raw, image, molecule='H2O', kp_expected=190.)
    assert image.is_file() and image.stat().st_size > 1000

    phi = np.linspace(1.4, 1.8, 31)
    rv = np.arange(-500., 501., 2.)
    grid = np.sin(rv[None, :] / 19. + phi[:, None] * 14.)
    image = tmp_path / 'phase_rv.png'
    plot_phase_rv(phi, rv, grid, image)
    assert image.is_file() and image.stat().st_size > 1000


def test_flat_ccf_cannot_be_called_sigma():
    with pytest.raises(ValueError, match='flat'):
        standardize_ccf(np.full((3, 5), 5.), np.array([True] * 5))
