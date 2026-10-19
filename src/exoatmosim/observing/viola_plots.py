"""Diagnostic figures for native VIOLA simulations, without recomputing spectra."""
from __future__ import annotations

from pathlib import Path
import numpy as np


def plot_viola_diagnostics(output_dir, wavelength_nm, stellar_e, model_e,
                           planet_wavelength_nm, planet_transmission, transit_weights):
    """Write observed-spectrum comparison from finished arrays.

    Uses deterministic, noise-free electron expectations; plots are diagnostics,
    not a claim of validated detector SNR or telluric PWV calibration.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from ..plotting import publication_style

    out = Path(output_dir) / 'diagnostics'
    out.mkdir(parents=True, exist_ok=True)
    w = np.asarray(wavelength_nm, dtype=float)
    star = np.asarray(stellar_e, dtype=float)
    model = np.asarray(model_e, dtype=float)
    weights = np.asarray(transit_weights, dtype=float)
    if (w.ndim != 1 or star.ndim != 2 or star.shape != model.shape or
        star.shape[1] != len(w) or star.shape[0] != len(weights) or
        np.any(np.diff(w) <= 0) or
        not all(np.isfinite(a).all() for a in (w, star, model, weights))):
        raise ValueError('Invalid VIOLA arrays for diagnostics')
    i = int(np.argmax(weights))
    with publication_style():
        fig, (ax, residual) = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                                       gridspec_kw={'height_ratios': [3, 1]})
        ax.plot(w, model[i], label='Star + tellurics + planet', lw=.75, alpha=.8)
        ax.plot(w, star[i], label='Star + tellurics (no planet)', lw=.75)
        ax.set(ylabel='Expected electrons / resolution element',
               title=f'VIOLA exposure {i}, transit weight = {weights[i]:.3f}')
        ax.legend(loc='best')
        ratio = np.divide(model[i], star[i], out=np.full(len(w), np.nan), where=star[i] > 0)
        residual.plot(w, ratio - 1, lw=.75)
        residual.axhline(0, lw=.6, color='grey')
        residual.set(xlabel='Air wavelength [nm]', ylabel='Planet / star − 1')
        fig.tight_layout()
        fig.savefig(out / 'viola_spectral_comparison.png', dpi=160)
        plt.close(fig)
    return out / 'viola_spectral_comparison.png'
