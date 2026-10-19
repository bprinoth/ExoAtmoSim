"""Diagnostic figures for native VIOLA simulations, without recomputing spectra."""
from __future__ import annotations

from pathlib import Path
import numpy as np


def plot_stellar_spectrum_comparison(output_dir, original_wave_nm, original_flux,
                                     broadened_wave_nm, broadened_flux, vsini_kms,
                                     *, instrument_resolution=None):
    """Compare unbroadened stellar flux with actual broadened stellar spectrum."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from ..plotting import publication_style

    original_wave_nm = np.asarray(original_wave_nm, float)
    original_flux = np.asarray(original_flux, float)
    broadened_wave_nm = np.asarray(broadened_wave_nm, float)
    broadened_flux = np.asarray(broadened_flux, float)
    if (original_wave_nm.ndim != 1 or broadened_wave_nm.ndim != 1 or
            original_wave_nm.shape != original_flux.shape or
            broadened_wave_nm.shape != broadened_flux.shape or
            not all(np.isfinite(a).all() for a in
                    (original_wave_nm, original_flux, broadened_wave_nm, broadened_flux))):
        raise ValueError('Invalid stellar rotation diagnostic arrays')
    out = Path(output_dir) / 'diagnostics' / 'stellar_spectrum_comparison.png'
    out.parent.mkdir(parents=True, exist_ok=True)
    if instrument_resolution is None:
        label = f'Rotationally broadened (vsini = {vsini_kms:g} km/s)'
    else:
        label = f'Broadened (vsini = {vsini_kms:g} km/s, R = {instrument_resolution:g})'
    with publication_style():
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(original_wave_nm, original_flux, lw=.75, alpha=.7, label='PHOENIX (unbroadened)')
        ax.plot(broadened_wave_nm, broadened_flux, lw=.75,
                label=label)
        ax.set(xlabel='Air wavelength [nm]',
               ylabel=r'Stellar flux [erg s$^{-1}$ cm$^{-2}$ nm$^{-1}$]',
               title='Unbroadened versus broadened stellar spectrum')
        ax.legend(loc='best')
        fig.tight_layout()
        fig.savefig(out, dpi=180)
        plt.close(fig)
    return out
