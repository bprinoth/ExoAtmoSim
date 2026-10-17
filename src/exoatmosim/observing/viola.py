"""VIOLA ETC stellar spectra already in the observatory rest frame.

VIOLA output contains telluric absorption and receives a *per-exposure*
observer-frame stellar radial velocity. Never shift or telluric-multiply it
again during order synthesis. ``signal_noise.csv`` flux is electrons/resel,
not PHOENIX physical surface flux; no cross-source flux normalization is implied.
"""
from __future__ import annotations

import csv
from pathlib import Path
import numpy as np


def read_viola_csv(path, *, wavelength_medium='air'):
    """Read VIOLA signal_noise.csv; require explicit AIR wavelength convention.

    VIOLA's documented CSV wavelength_um does not encode air/vacuum metadata.
    Do not guess: caller must confirm the VIOLA run uses air wavelengths.
    """
    if str(wavelength_medium).lower() != 'air':
        raise ValueError('VIOLA wavelength medium is not confirmed AIR; convert explicitly before using air-nm orders')
    path = Path(path)
    with path.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    required = {'wavelength_um', 'stellar_e_per_resel', 'snr_per_resel'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f'{path}: expected VIOLA CSV columns {sorted(required)}')
    wl = np.array([float(r['wavelength_um']) * 1000.0 for r in rows])
    flux = np.array([float(r['stellar_e_per_resel']) for r in rows])
    snr = np.array([float(r['snr_per_resel']) for r in rows])
    if (not np.isfinite(wl).all() or not np.isfinite(flux).all()
            or not np.isfinite(snr).all() or np.any(np.diff(wl) <= 0)
            or np.any(wl <= 0) or np.any(flux < 0) or np.any(snr < 0)):
        raise ValueError(f'{path}: invalid or unordered VIOLA spectrum')
    return wl, flux, snr


def viola_observer_spectra(schedule, observer_rv_kms, cfg, output_root, *,
                           run_viola=None):
    """Produce exact per-exposure VIOLA spectra and SNR without time interpolation.

    ``run_viola`` is injectable to test without installed external VIOLA.
    The default imports ``run_viola_etc_cli`` from the user-installed package.
    One run per exposure, never grouped by rounded altitude or RV.
    """
    rv = np.asarray(observer_rv_kms, dtype=float)
    altitude = np.asarray(schedule.altitude_deg, dtype=float)
    if rv.shape != altitude.shape or not np.isfinite(rv).all():
        raise ValueError('VIOLA needs one finite observer-frame stellar RV per exposure')
    if np.any(altitude <= 0):
        raise ValueError('VIOLA requires target above horizon at every exposure')
    band = str(cfg.get('viola_band', 'K')).upper()
    if band not in ('Y', 'J', 'H', 'K'):
        raise ValueError('VIOLA band must be Y, J, H, or K')
    medium = str(cfg.get('viola_wavelength_medium', '')).lower()
    if medium != 'air':
        raise ValueError('Set viola_wavelength_medium air only after verifying VIOLA wavelength_um output uses air wavelengths')
    # Y-band magnitude may not be present in NASA Archive defaults.
    if band == 'Y' and cfg.get('Ymag') is None:
        raise ValueError('VIOLA_Y requires an explicit calibrated Ymag in config; do not substitute Jmag')
    mag = cfg.get_float(f'{band}mag')
    teff = cfg.get_float('Teff')
    texp = cfg.get_float('exposure_time')
    pwv_mean = cfg.get_float('viola_pwv_mm', 2.5)
    pwv_sigma = cfg.get_float('viola_pwv_sigma_mm', 0.0)
    pwv_timescale = cfg.get_float('viola_pwv_timescale_minutes', 30.0)
    pwv_seed = cfg.get_int('viola_pwv_seed', 42)
    seeing = cfg.get_float('viola_seeing_arcsec', 1.0)
    if not np.isfinite([mag, teff, texp, pwv_mean, pwv_sigma, pwv_timescale, seeing]).all() or texp <= 0 or pwv_mean <= 0 or pwv_sigma < 0 or pwv_timescale <= 0 or seeing <= 0:
        raise ValueError('Invalid VIOLA configuration')
    if run_viola is None:
        try:
            from run_viola_etc_cli import run_viola_etc as run_viola
        except ImportError as exc:
            raise ImportError('VIOLA not importable: add viola-etc to PYTHONPATH or install it') from exc
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    if hasattr(schedule, 'midpoint_jd_utc'):
        times_min = (np.asarray(schedule.midpoint_jd_utc, float) -
                     float(schedule.midpoint_jd_utc[0])) * 1440.0
    else:
        times_min = np.arange(len(rv), dtype=float) * (texp / 60.0)
    if not np.isfinite(times_min).all() or np.any(np.diff(times_min) <= 0):
        raise ValueError('VIOLA requires strictly increasing exposure midpoint times')
    # A correlated *PWV* process (mm), distinct from TelFit's RH percentage.
    rng = np.random.default_rng(pwv_seed)
    standardized = np.zeros(len(rv), dtype=float)
    if len(rv):
        standardized[0] = rng.normal()
    for j in range(1, len(rv)):
        rho = np.exp(-(times_min[j] - times_min[j-1]) / pwv_timescale)
        standardized[j] = rho * standardized[j-1] + np.sqrt(1-rho*rho)*rng.normal()
    pwv_per_exposure = np.maximum(0.001, pwv_mean + pwv_sigma * standardized)
    with (output_root / 'viola_conditions.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['exposure_id', 'midpoint_minutes', 'altitude_deg',
                         'stellar_observer_rv_kms', 'pwv_mm'])
        for j, (minute, alt, vel, water) in enumerate(zip(times_min, altitude, rv, pwv_per_exposure)):
            writer.writerow([j, f'{minute:.9f}', f'{alt:.8f}', f'{vel:.8f}', f'{water:.8f}'])
    waves, fluxes, snrs = None, [], []
    for i, (alt, velocity, pwv) in enumerate(zip(altitude, rv, pwv_per_exposure)):
        folder = output_root / f'exposure_{i:04d}' / band
        folder.mkdir(parents=True, exist_ok=True)
        result = run_viola(
            band=band, magnitude=mag, t_star_k=teff,
            target_altitude_deg=float(alt), exposure_time_s=texp,
            seeing_arcsec=seeing, pwv_mm=pwv, enable_molecules=False,
            radial_velocity_km_s=float(velocity), output_folder=str(folder),
        )
        result_folder = Path(result) if result is not None else folder
        wave, flux, snr = read_viola_csv(result_folder / 'signal_noise.csv',
                                         wavelength_medium=medium)
        if waves is None:
            waves = wave
        elif wave.shape != waves.shape or not np.allclose(wave, waves, rtol=0, atol=1e-8):
            raise ValueError(f'VIOLA exposure {i}: output wavelength grid differs; cannot blindly stack; align explicitly')
        fluxes.append(flux)
        snrs.append(snr)
    return waves, np.asarray(fluxes), np.asarray(snrs)
