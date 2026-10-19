"""Diagnostic figures for native VIOLA simulations, without recomputing spectra."""
from __future__ import annotations

from pathlib import Path
import numpy as np


def plot_viola_diagnostics(output_dir, wavelength_nm, stellar_e, model_e,
                           planet_wavelength_nm, planet_transmission, transit_weights):
    """Write input-planet and observed-spectrum comparisons from finished arrays.

    Uses deterministic, noise-free electron expectations; plots are diagnostics,
    not a claim of validated detector SNR or telluric PWV calibration.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    out = Path(output_dir) / 'diagnostics'
    out.mkdir(parents=True, exist_ok=True)
    w = np.asarray(wavelength_nm, dtype=float)
    star = np.asarray(stellar_e, dtype=float)
    model = np.asarray(model_e, dtype=float)
    pw = np.asarray(planet_wavelength_nm, dtype=float)
    pf = np.asarray(planet_transmission, dtype=float)
    weights = np.asarray(transit_weights, dtype=float)
    if (w.ndim != 1 or star.ndim != 2 or star.shape != model.shape or
        star.shape[1] != len(w) or star.shape[0] != len(weights) or
        pw.ndim != 1 or pf.shape != pw.shape or len(pw) < 2 or
        np.any(np.diff(w) <= 0) or np.any(np.diff(pw) <= 0) or
        not all(np.isfinite(a).all() for a in (w, star, model, pw, pf, weights))):
        raise ValueError('Invalid VIOLA arrays for diagnostics')
    with plt.rc_context({'figure.figsize': (10, 4.5)}):
        fig, ax = plt.subplots()
        ax.plot(pw, pf, linewidth=.75)
        ax.set(xlabel='Air wavelength [nm]', ylabel='Planet transmission factor',
               title='Input planetary transmission spectrum (before observation)')
        fig.tight_layout()
        fig.savefig(out / 'planet_transmission.png', dpi=160)
        plt.close(fig)

    i = int(np.argmax(weights))
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
    return (out / 'planet_transmission.png', out / 'viola_spectral_comparison.png')
