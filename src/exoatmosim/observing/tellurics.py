"""Optional independent TelFit interface, with explicit AIR-nm wavelengths.

TelFit/LBLRTM itself is a user-installed optional dependency. For physical
transmission, evaluate at every exposure's actual target zenith angle.
No assumed constant airmass, no extrapolation and no silent transmission=1.
"""
from __future__ import annotations
import numpy as np


def zenith_angles(altitude_deg):
    alt = np.asarray(altitude_deg, dtype=float)
    if alt.ndim != 1 or not np.isfinite(alt).all() or np.any((alt <= 0) | (alt > 90)):
        raise ValueError('TelFit requires a visible target (0 < altitude <= 90 deg) at every exposure')
    return 90.0 - alt


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


def telfit_transmission(schedule, wavelength_air_nm, *, pressure_hpa=600.,
                        temperature_k=273., latitude_deg=-24.627,
                        elevation_km=2.635, resolving_power=100000., modeler_factory=None, conditions=None,
                        ringing_tolerance=2e-4):
    """Return (air-nm model wavelengths, Nexp x Nwave transmission).

    Negative transmission is always clipped to zero with a warning,
    regardless of its magnitude. Non-finite values remain an error. Values
    above unity remain subject to ringing_tolerance (default 2e-4).
    Clipping is a numerical recovery, not proof of physical absorption.

    ``modeler_factory`` supports isolated tests with a fake modeler.
    TelFit's ``MakeModel(..., vac2air=True)`` output is treated as *AIR nm*;
    this must be confirmed with the local TelFit install before publication.
    Requested instrument wavelengths are AIR nm. TelFit accepts *vacuum*
    wavenumber bounds even when ``vac2air=True`` is requested for its output.
    Convert bounds to vacuum first and include a small coverage margin.
    """
    if not np.isfinite(ringing_tolerance) or not 0 <= ringing_tolerance <= 2e-4:
        raise ValueError('ringing_tolerance must be between 0 and 2e-4')
    wavelengths = np.asarray(wavelength_air_nm, float)
    if (wavelengths.ndim != 1 or wavelengths.size < 2
            or not np.isfinite(wavelengths).all() or np.any(np.diff(wavelengths) <= 0)
            or wavelengths[0] <= 0):
        raise ValueError('Wavelengths must be finite increasing air nm')
    if not np.isfinite([pressure_hpa, temperature_k, latitude_deg, elevation_km, resolving_power]).all() or pressure_hpa <= 0 or temperature_k <= 0 or resolving_power <= 0:
        raise ValueError('Invalid TelFit atmospheric or instrument settings')
    angles = zenith_angles(schedule.altitude_deg)
    if modeler_factory is None:
        try:
            from telfit import Modeler
        except ImportError as exc:
            raise ImportError('Tellurics requested but TelFit could not be imported. Install/configure TelFit and LBLRTM first.') from exc
        modeler_factory = Modeler
    # TelFit interprets wavenumber limits in vacuum cm^-1. Convert the
    # desired AIR-nm limits to VACUUM-nm *before* specifying those limits.
    # Padding avoids interpolation/boundary rounding excluding edge pixels.
    pad_nm = max(0.02, 2.0 * float(np.median(np.diff(wavelengths))))
    bounds_vac_nm = air_to_vacuum_nm([wavelengths[0] - pad_nm,
                                     wavelengths[-1] + pad_nm])
    lowfreq_cm1 = 1e7 / bounds_vac_nm[1]
    highfreq_cm1 = 1e7 / bounds_vac_nm[0]
    if conditions is not None:
        if len(conditions) != len(angles):
            raise ValueError('Atmospheric conditions must have one row per exposure')
        if not np.allclose(np.asarray(conditions['zenith_angle_deg'], float), angles, atol=1e-7):
            raise ValueError('Atmospheric conditions zenith angles disagree with scheduler')
    spectra = []
    reference = wavelengths.copy()
    for i, angle in enumerate(angles):
        modeler = modeler_factory()
        exposure_state = ({'pressure': pressure_hpa, 'temperature': temperature_k} if conditions is None else {
            'humidity': float(conditions['humidity_percent'][i]),
            'pressure': float(conditions['pressure_hpa'][i]),
            'temperature': float(conditions['temperature_k'][i]),
        })
        model = modeler.MakeModel(**exposure_state, lowfreq=lowfreq_cm1,
                                  highfreq=highfreq_cm1,
                                  angle=float(angle), resolution=resolving_power,
                                  vac2air=True, lat=latitude_deg, alt=elevation_km)
        x = np.asarray(model.x, float)
        y = np.asarray(model.y, float)
        if x.ndim != 1 or x.size < 2 or x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError('TelFit returned invalid wavelength or transmission arrays')
        if np.any(y > 1+ringing_tolerance):
            above = int(np.count_nonzero(y > 1.0+ringing_tolerance))
            raise ValueError(
                f'TelFit exposure {i}: transmission outside [0,1] beyond ringing_tolerance={ringing_tolerance:g}; '
                f'min={float(np.min(y)):.9e}, max={float(np.max(y)):.9e}, '
                f'above_tolerance={above}. '
                'Inspect LBLRTM output, wavelength grid, and model conditions; '
                'check the unmodified LBLRTM spectrum before interpreting clipped regions.')
        if x[0] > x[-1]:
            x, y = x[::-1], y[::-1]
        if np.any(np.diff(x) <= 0) or x[0] > wavelengths[0] or x[-1] < wavelengths[-1]:
            raise ValueError('TelFit AIR nm model fails to cover instrument grid; verify units, vac2air and buffer')
        # TelFit/LBLRTM may return slightly different sample grids for each
        # exposure. Resample every model onto the caller's requested AIR-nm grid.
        # Coverage is checked above; np.interp never extrapolates here.
        if x[0] > reference[0] or x[-1] < reference[-1]:
            raise ValueError(f'TelFit exposure {i}: output grid does not cover reference grid')
        if np.any((y < 0.) | (y > 1.)):
            import warnings
            warnings.warn(
                f'TelFit exposure {i}: clamping {np.count_nonzero(y < 0.)} negative '
                f'transmission pixels (minimum {np.min(y):.6g}) and '
                f'{np.count_nonzero(y > 1.)} pixels above unity into [0, 1]. '
                'This is numerical post-processing; inspect the telluric model if excursions are substantial.',
                RuntimeWarning)
        spectra.append(np.interp(reference, x, np.clip(y, 0., 1.)))
    return reference, np.asarray(spectra)


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
