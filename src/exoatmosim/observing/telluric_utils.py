"""Wavelength conversion and telluric-diagnostic plotting helpers.

Telluric transmission is loaded exclusively from external HDF5 grids.
"""
from __future__ import annotations
import numpy as np

def vacuum_to_air_nm(wavelength_vac_nm):
    """Convert vacuum wavelengths to standard-air nm (Morton/Edlen relation).

    The expression uses vacuum wavenumber in inverse micrometres. It is
    intended for optical/near-infrared wavelengths, including CRIRES+ K band.
    """
    vac = np.asarray(wavelength_vac_nm, dtype=float)
    if not np.all(np.isfinite(vac)) or np.any(vac <= 0):
        raise ValueError('Vacuum wavelengths must be positive and finite')
    sigma2 = (1000. / vac)**2
    refractive_index = 1. + 1e-8 * (8342.13 + 2406030. / (130. - sigma2)
                                      + 15997. / (38.9 - sigma2))
    return vac / refractive_index


def air_to_vacuum_nm(wavelength_air_nm):
    """Numerically invert the standard-air dispersion relation, in nm."""
    air = np.asarray(wavelength_air_nm, dtype=float)
    if not np.all(np.isfinite(air)) or np.any(air < 200.):
        raise ValueError('Air wavelengths must be finite and >= 200 nm')
    vac = air * 1.00028
    for _ in range(6):
        vac = vac * air / vacuum_to_air_nm(vac)
    return vac


def plot_tellurics(wavelength_nm, transmission, schedule, outfile):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from exoatmosim.plotting import publication_style
    from pathlib import Path
    path = Path(outfile)
    path.parent.mkdir(parents=True, exist_ok=True)
    transmission = np.asarray(transmission)
    with publication_style():
        fig, axes = plt.subplots(2, 1, figsize=(12, 7), constrained_layout=True)
        for i in sorted(set((0, len(transmission)//2, len(transmission)-1))):
            axes[0].plot(wavelength_nm, transmission[i], lw=0.8, label=f'Exposure {i}, X={schedule.airmass[i]:.2f}')
        axes[0].set(xlabel='Air wavelength [nm]', ylabel='Transmission', ylim=(-0.02, 1.05))
        axes[0].legend()
        image = axes[1].imshow(transmission, aspect='auto', origin='lower', extent=(wavelength_nm[0], wavelength_nm[-1], 0, len(transmission)))
        axes[1].set(xlabel='Air wavelength [nm]', ylabel='Exposure index')
        fig.colorbar(image, ax=axes[1], label='Transmission')
        fig.savefig(path, dpi=180)
        plt.close(fig)
    return path
