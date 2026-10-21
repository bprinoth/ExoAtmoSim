"""Native, publication-style cross-correlation diagnostic figures.

FITS CCF arrays are unchanged; only display coordinates are normalized.
The plotted sigma-like score is the background-standardized CCF, not a
calibrated false-alarm significance. Background excludes |Vsys| < 50 km/s.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np


def _check_grid(x, y, array):
    x, y, array = np.asarray(x, float), np.asarray(y, float), np.asarray(array, float)
    if (x.ndim != 1 or y.ndim != 1 or array.shape != (len(y), len(x)) or
            len(x) < 3 or len(y) < 2 or not np.all(np.isfinite(x)) or
            not np.all(np.isfinite(y)) or
            np.any(np.diff(x) <= 0) or np.any(np.diff(y) <= 0)):
        raise ValueError('Invalid CCF map or velocity coordinates')
    return x, y, array


def standardize_ccf(data, background_mask):
    """Standardize with a single off-trail background mean and standard deviation.

    A global background preserves the relative contrasts across Kp. A nearly
    constant map is scientifically undefined as a significance estimate.
    """
    data = np.asarray(data, float)
    background = data[:, background_mask]
    finite = background[np.isfinite(background)]
    if finite.size < 4:
        raise ValueError('Insufficient finite off-trail CCF samples')
    scale = float(np.std(finite, ddof=1))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps * max(1., abs(float(np.mean(finite)))):
        raise ValueError('Cannot normalize a flat CCF map: zero background scatter')
    return (data - float(np.mean(finite))) / scale


def plot_kp_vsys(kp, vsys, ccf, destination, *, molecule=None,
                 kp_expected=None, vsys_expected=0.):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .plotting import publication_style

    vsys, kp, ccf = _check_grid(vsys, kp, ccf)
    keep = kp >= 0
    if np.count_nonzero(keep) < 2:
        raise ValueError('Kp axis must contain at least two nonnegative samples')
    noise_region = np.abs(vsys - vsys_expected) >= 50.
    if np.count_nonzero(noise_region) < 2:
        raise ValueError('Need at least two Vsys samples outside ±50 km/s for CCF noise')
    sigma = standardize_ccf(ccf, noise_region)
    positive = sigma[keep]
    kp_positive = kp[keep]
    slice_kp = kp_expected if kp_expected is not None else kp_positive[np.nanargmax(np.nanmax(positive, axis=1))]
    index = int(np.argmin(abs(kp_positive - slice_kp)))
    curve = positive[index]
    finite = positive[np.isfinite(positive)]
    if not finite.size:
        raise ValueError('No finite standardized CCF map')
    limit = max(3., float(np.percentile(np.abs(finite), 99.5)))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with publication_style():
        fig = plt.figure(figsize=(10, 7), layout='constrained')
        gs = fig.add_gridspec(2, 2, width_ratios=[28, 1], height_ratios=[3, 1],
                             hspace=.03, wspace=.025)
        ax = fig.add_subplot(gs[0, 0])
        slice_ax = fig.add_subplot(gs[1, 0], sharex=ax)
        cax = fig.add_subplot(gs[0, 1])
        image = ax.pcolormesh(vsys, kp_positive, positive, cmap='magma',
                              shading='auto', vmin=-limit, vmax=limit, rasterized=True)
        fig.colorbar(image, cax=cax, label='CCF / background scatter')
        if kp_expected is not None:
            ax.axhline(kp_expected, lw=1.2, color='white', ls='--', alpha=.85)
            ax.plot(vsys_expected, kp_expected, marker='x', color='white', ms=8, mew=1.8)
        ax.axvline(vsys_expected, color='white', lw=1., ls=':', alpha=.8)
        ax.set(ylabel=r'$K_p$ (km/s)', ylim=(0, float(kp_positive[-1])))
        if molecule:
            ax.text(.02, .97, molecule, transform=ax.transAxes, ha='left',
                    va='top', color='white', fontsize=14, weight='bold')
        ax.tick_params(labelbottom=False)
        slice_ax.plot(vsys, curve, color='#bd2943', lw=1.5)
        slice_ax.axhline(0, lw=.8, color='0.5')
        slice_ax.axvline(vsys_expected, color='0.4', lw=1., ls=':')
        slice_ax.set(xlabel=r'$V_{\\rm sys}$ (km/s)', ylabel='CCF / noise',
                     title=rf'Slice at $K_p={kp_positive[index]:.1f}$ km/s')
        slice_ax.title.set_fontsize(11)
        slice_ax.tick_params(top=False)
        fig.savefig(destination, dpi=190, bbox_inches='tight')
        plt.close(fig)
    return destination


def plot_phase_rv(phi, rv, ccf, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .plotting import publication_style
    rv = np.asarray(rv, float)
    phi = np.asarray(phi, float)
    ccf = np.asarray(ccf, float)
    if (phi.ndim != 1 or rv.ndim != 1 or ccf.shape != (len(phi), len(rv))
            or not np.isfinite(phi).all() or not np.isfinite(rv).all()
            or np.any(np.diff(rv) <= 0)):
        raise ValueError('Invalid phase–RV data shape')
    # Phase need not be monotonically increasing after eccentric geometry.
    order = np.argsort(phi)
    phi, ccf = phi[order], ccf[order]
    noise = abs(rv) >= 250.
    if np.count_nonzero(noise) < 2:
        noise = np.ones(rv.shape, bool)
    scaled = standardize_ccf(ccf, noise)
    bound = max(3., float(np.nanpercentile(np.abs(scaled), 99.5)))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with publication_style():
        fig, ax = plt.subplots(figsize=(10, 4.7), layout='constrained')
        image = ax.pcolormesh(rv, phi, scaled, shading='auto', cmap='magma',
                              vmin=-bound, vmax=bound, rasterized=True)
        fig.colorbar(image, ax=ax, pad=.015, fraction=.028,
                     label='CCF / background scatter')
        ax.set(xlabel='Radial velocity (km/s)', ylabel='True anomaly (rad)')
        fig.savefig(destination, dpi=190, bbox_inches='tight')
        plt.close(fig)
    return destination
