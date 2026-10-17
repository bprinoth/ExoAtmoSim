"""VIOLA detector-electron simulation on a fixed, externally produced wavelength grid.

VIOLA is an instrument, not a PHOENIX input to another spectrograph. The ETC
returns already Doppler-shifted and telluric-absorbed stellar electrons/pixel.
No second stellar RV, stellar LSF, or telluric multiplication is permitted.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import numpy as np

C_KMS = 299792.458


@dataclass(frozen=True)
class ViolaDetectorInput:
    wavelength_nm: np.ndarray
    stellar_e_pixel: np.ndarray
    stellar_e_resel: np.ndarray
    snr_resel: np.ndarray
    total_noise_e_resel: np.ndarray
    valid_noise: np.ndarray | None = None


def read_viola_detector_csv(path, *, wavelength_medium='air') -> ViolaDetectorInput:
    """Parse exact, native VIOLA arrays without changing their detector sampling."""
    if wavelength_medium.lower() != 'air':
        raise ValueError('VIOLA wavelength_um medium must first be confirmed AIR')
    with Path(path).open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        required = ('wavelength_um', 'stellar_e_per_pixel', 'stellar_e_per_resel',
                    'snr_per_resel', 'total_noise_e_per_resel')
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f'VIOLA signal_noise.csv must contain {required}')
        rows = list(reader)
    if len(rows) < 2:
        raise ValueError('VIOLA requires at least two detector samples')
    arrays = [np.array([float(row[col]) for row in rows], dtype=float) for col in required]
    wl, pix, res, snr, noise = arrays
    # A zero-SNR sample may have NaN *derived* resolution-element noise.
    # Do not discard its wavelength or manufacture an uncertainty.
    if (not all(np.isfinite(a).all() for a in (wl, pix, res, snr)) or
            np.any(np.diff(wl) <= 0) or np.any(wl <= 0) or
            np.any(pix < 0) or np.any(res < 0) or np.any(snr < 0)):
        raise ValueError(f'{path}: invalid physical detector values or wavelengths')
    usable = np.isfinite(noise) & (noise >= 0) & (snr > 0)
    # Conversely, nonfinite / negative noise is NOT permitted at positive SNR.
    if np.any((snr > 0) & ~usable) or np.any((snr == 0) & np.isfinite(noise) & (noise < 0)):
        raise ValueError(f'{path}: invalid resolution-element noise at positive SNR')
    # Preserve nonzero noise reported where the rounded SNR is exactly zero,
    # but mask those samples until their ETC noise provenance is resolved.
    cleaned_noise = noise.copy()
    cleaned_noise[~usable] = np.nan
    return ViolaDetectorInput(wl * 1000., pix, res, snr, cleaned_noise, usable)


def pixel_variance(etc: ViolaDetectorInput):
    """Approximate per-pixel ETC variance; requires resel-to-pixel sampling.

    Photon noise from stellar electrons/pixel is explicit. Remaining ETC
    resolution-element variance (sky, thermal, dark, read, etc.) is spread
    uniformly over N_eff = stellar_e_per_resel / stellar_e_per_pixel.
    This is an approximation pending access to VIOLA's detector-pixel variance.
    """
    pix, res = etc.stellar_e_pixel, etc.stellar_e_resel
    valid = (etc.valid_noise if etc.valid_noise is not None else
             np.isfinite(etc.total_noise_e_resel))
    valid = np.asarray(valid, bool)
    if valid.shape != pix.shape:
        raise ValueError('VIOLA noise quality mask differs from detector grid')
    result = np.full(pix.shape, np.nan, dtype=float)
    if not valid.any():
        return result
    if np.any((pix[valid] <= 0) | (res[valid] <= 0)):
        raise ValueError('Cannot infer VIOLA per-pixel noise where stellar counts vanish')
    n_eff = res[valid] / pix[valid]
    if np.any(~np.isfinite(n_eff)) or np.any(n_eff < 1):
        raise ValueError('VIOLA resel-to-pixel sampling is inconsistent')
    res_variance = etc.total_noise_e_resel[valid]**2
    if np.any(res_variance + 1.e-7 < res[valid]):
        raise ValueError('VIOLA resolution-element noise is below stellar shot noise')
    result[valid] = pix[valid] + np.maximum(0., res_variance - res[valid]) / n_eff
    return result


def inject_viola_planet(etc_rows, planet_wave_nm, planet_factor, planet_rv_kms,
                        transit_weight, *, seed=42, add_noise=True):
    """Return native per-exposure star, planet-injected mean, error and sampled counts.

    `planet_factor` is dimensionless transmission, unity outside atmosphere;
    no additional stellar shift or telluric absorption is performed.
    """
    rows = list(etc_rows)
    rv = np.asarray(planet_rv_kms, float)
    weight = np.asarray(transit_weight, float)
    pw = np.asarray(planet_wave_nm, float)
    pf = np.asarray(planet_factor, float)
    if not rows or rv.shape != (len(rows),) or weight.shape != rv.shape:
        raise ValueError('One VIOLA ETC spectrum, RV and transit weight per exposure required')
    if (pw.ndim != 1 or pf.shape != pw.shape or len(pw) < 2 or
            np.any(np.diff(pw) <= 0) or not np.isfinite(pf).all() or
            not np.isfinite(rv).all() or not np.isfinite(weight).all() or
            np.any((weight < 0) | (weight > 1))):
        raise ValueError('Invalid planet spectrum, velocity, or occultation weights')
    grid = rows[0].wavelength_nm
    stellar, modeled, errors, noisy = [], [], [], []
    rng = np.random.default_rng(seed)
    for i, etc in enumerate(rows):
        if etc.wavelength_nm.shape != grid.shape or not np.allclose(etc.wavelength_nm, grid, rtol=0, atol=1e-8):
            raise ValueError(f'VIOLA exposure {i} has different sampling: cannot stack without validated resampling')
        at = grid / (1. + rv[i]/C_KMS)
        if at[0] < pw[0] or at[-1] > pw[-1]:
            raise ValueError(f'Planetary model does not cover Doppler-shifted VIOLA exposure {i}')
        star = etc.stellar_e_pixel
        planet = np.interp(at, pw, pf)
        expected = star * (1. + weight[i]*(planet - 1.))
        if np.any(expected < 0):
            raise ValueError('Negative planetary-injected electron expectation')
        base_var = pixel_variance(etc)
        # Stellar photon variance follows the planet-modulated expectation.
        variance = np.maximum(0., base_var - star) + expected
        sigma = np.sqrt(variance)
        observed = expected + rng.normal(size=len(grid))*sigma if add_noise else expected.copy()
        stellar.append(star.copy()); modeled.append(expected); errors.append(sigma); noisy.append(observed)
    return (grid.copy(), np.asarray(stellar), np.asarray(modeled),
            np.asarray(errors), np.asarray(noisy))


def make_viola_observation(schedule, star_rv_observer_kms, planet_rv_observer_kms,
                           transit_weight, cfg, planet_wave_nm, planet_factor,
                           output_dir, *, seed=42, overwrite=False,
                           run_viola=None, resume_from=None):
    """Execute independently for each exposure; keep science products, not ETC dumps.

    Writes native air-nm order 0, observer-frame star electrons/pixel,
    model electrons/pixel, and predicted noise electrons/pixel. No synthetic
    detector order boundaries are fabricated. VIOLA CLI is run in a temporary
    directory; each output is read and discarded before the next exposure.
    """
    import json
    import os
    import tempfile

    band = str(cfg.get('instrument')).upper().removeprefix('VIOLA_')
    if band not in ('K', 'H', 'Y'):
        raise ValueError('Select instrument VIOLA_K, VIOLA_H, or VIOLA_Y')
    if cfg.get_bool('tellurics', False):
        raise ValueError('VIOLA stellar output already contains tellurics: set tellurics False')
    if cfg.get('viola_wavelength_medium', '').lower() != 'air':
        raise ValueError('Confirm VIOLA wavelengths are AIR in viola_wavelength_medium')
    if band == 'Y' and cfg.get('Ymag') is None:
        raise ValueError('VIOLA_Y requires explicit calibrated Ymag')
    if run_viola is None and resume_from is None:
        try:
            from run_viola_etc_cli import run_viola_etc as run_viola
        except ImportError as exc:
            raise ImportError('Could not import run_viola_etc_cli; add its directory to PYTHONPATH') from exc
    rv_star = np.asarray(star_rv_observer_kms, float)
    rv_planet = np.asarray(planet_rv_observer_kms, float)
    altitude = np.asarray(schedule.altitude_deg, float)
    if (rv_star.shape != altitude.shape or rv_planet.shape != altitude.shape or
            not np.isfinite(rv_star).all() or np.any(altitude <= 0)):
        raise ValueError('Invalid exposure geometry or stellar RV inputs')
    time_minutes = ((np.asarray(schedule.midpoint_jd_utc, float) - schedule.midpoint_jd_utc[0]) * 1440.
                    if hasattr(schedule, 'midpoint_jd_utc') else
                    np.arange(len(altitude))*cfg.get_float('exposure_time')/60.)
    if len(time_minutes) > 1 and np.any(np.diff(time_minutes) <= 0):
        raise ValueError('Exposure times must strictly increase')
    mean = cfg.get_float('viola_pwv_mm', 2.5)
    scatter = cfg.get_float('viola_pwv_sigma_mm', 0.)
    tau = cfg.get_float('viola_pwv_timescale_minutes', 30.)
    if not np.isfinite([mean, scatter, tau]).all() or mean <= 0 or scatter < 0 or tau <= 0:
        raise ValueError('Invalid VIOLA PWV configuration')
    rng = np.random.default_rng(cfg.get_int('viola_pwv_seed', seed))
    pwv_noise = np.empty(len(altitude))
    pwv_noise[0] = rng.normal()
    for i in range(1, len(altitude)):
        rho = np.exp(-(time_minutes[i]-time_minutes[i-1])/tau)
        pwv_noise[i] = rho*pwv_noise[i-1] + np.sqrt(max(0., 1.-rho**2))*rng.normal()
    pwv = np.maximum(0.001, mean + scatter*pwv_noise)
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    planned = [output_dir / n for n in ('wave_0.fits', 'order_0.fits',
                  'stellar_flux_0.fits', 'model_flux_0.fits', 'noise_0.fits',
                  'viola_conditions.csv', 'noise_quality_0.fits', f'VIOLA_{band}_orders_nm.npy', 'viola_manifest.json')]
    if not overwrite and any(p.exists() for p in planned):
        raise FileExistsError('VIOLA outputs already exist; pass --overwrite')
    kwargs = dict(band=band, magnitude=cfg.get_float(f'{band}mag'),
                  t_star_k=cfg.get_float('Teff'),
                  exposure_time_s=cfg.get_float('exposure_time'),
                  seeing_arcsec=cfg.get_float('viola_seeing_arcsec', 1.),
                  enable_molecules=False)
    from astropy.io import fits
    rows = []
    if resume_from is not None:
        # Read the real, previously completed ETC exposures. No VIOLA calls.
        source = Path(resume_from).resolve()
        if not source.is_dir():
            raise FileNotFoundError(f'VIOLA recovery directory does not exist: {source}')
        for i in range(len(altitude)):
            csv_path = source / f'exposure_{i:04d}' / band / 'signal_noise.csv'
            if not csv_path.is_file():
                raise FileNotFoundError(f'Missing VIOLA exposure {i}: {csv_path}')
            rows.append(read_viola_detector_csv(csv_path))
        # The original VIOLA conditions are necessary for physical provenance.
        # Do NOT pretend regenerated PWV values match prior ETC calls.
        original_conditions = source / 'viola_conditions.csv'
        if original_conditions.is_file():
            with original_conditions.open(newline='', encoding='utf8') as fh:
                prior = list(csv.DictReader(fh))
            if len(prior) != len(rows):
                raise ValueError('Recovered VIOLA condition count differs from exposure count')
            for i, record in enumerate(prior):
                if int(record['exposure_id']) != i:
                    raise ValueError(f'Recovered VIOLA exposure IDs mismatch at {i}')
                if not np.isclose(float(record['altitude_deg']), altitude[i], atol=1e-4):
                    raise ValueError(f'Recovered VIOLA altitude mismatch at exposure {i}')
                if not np.isclose(float(record['stellar_observer_rv_kms']), rv_star[i], atol=1e-4):
                    raise ValueError(f'Recovered stellar RV mismatch at exposure {i}')
            pwv = np.array([float(r['pwv_mm']) for r in prior])
        else:
            raise FileNotFoundError('Cannot scientifically validate recovery without the original '
                                    f'VIOLA conditions: {original_conditions}')
    else:
        # The ETC changes CWD. Restore it after every call. Temporary files are
        # discarded only after their required arrays have been extracted.
        with tempfile.TemporaryDirectory(prefix='exoatmosim_viola_') as td:
            for i in range(len(altitude)):
                folder = Path(td) / f'exp_{i:04d}'
                prev = Path.cwd()
                try:
                    result = run_viola(**kwargs, n_exposures=1,
                               target_altitude_deg=float(altitude[i]),
                               radial_velocity_km_s=float(rv_star[i]),
                               pwv_mm=float(pwv[i]), output_folder=str(folder))
                finally:
                    os.chdir(prev)
                actual_dir = Path(result) if result is not None else folder
                rows.append(read_viola_detector_csv(actual_dir / 'signal_noise.csv'))
    grid, star, model, sigma, observed = inject_viola_planet(
         rows, planet_wave_nm, planet_factor, rv_planet, transit_weight, seed=seed)
    wave = np.broadcast_to(grid, observed.shape)
    quality = np.asarray([np.isfinite(row.total_noise_e_resel) & (row.snr_resel > 0)
                          for row in rows], dtype=np.uint8)
    for filename, data in (('noise_quality_0.fits', quality), ('wave_0.fits', wave), ('order_0.fits', observed),
                 ('stellar_flux_0.fits', star), ('model_flux_0.fits', model),
                 ('noise_0.fits', sigma)):
        fits.writeto(output_dir / filename, data, overwrite=overwrite,
                     header=fits.Header({'BUNIT': ('1' if filename == 'noise_quality_0.fits' else
                                                  'nm' if filename == 'wave_0.fits' else 'electron'),
                                        'INSTRUME': f'VIOLA_{band}', 'WAVEMED': 'AIR'}))
    # Store exact measured samples locally, but do not mutate the installed package.
    np.save(output_dir / f'VIOLA_{band}_orders_nm.npy', grid[None, :])
    with (output_dir / 'viola_conditions.csv').open('w', newline='', encoding='utf8') as f:
        writer = csv.writer(f)
        writer.writerow(('exposure_id','midpoint_jd_utc','altitude_deg','stellar_rv_observer_kms','pwv_mm'))
        for i in range(len(altitude)):
            writer.writerow((i, float(schedule.midpoint_jd_utc[i]),float(altitude[i]),
                             float(rv_star[i]),float(pwv[i])))
    manifest = {'instrument': f'VIOLA_{band}', 'wavelength_medium': 'air',
        'wavelength_unit': 'nm', 'flux_unit': 'electron_per_pixel',
        'n_exposures': len(altitude), 'n_pixels': len(grid),
        'star_doppler_and_tellurics_included_by_viola': True,
        'additional_stellar_rv_or_telfit': False,
        'noise_model': 'Gaussian approximation for valid SNR pixels; zero-SNR pixels have NaN uncertainty and NaN noisy flux, with noise_quality_0.fits=0; resel-to-pixel allocation is approximate',
        'n_masked_noise_samples': int(np.size(quality) - np.count_nonzero(quality)),
        'wavelength_grid_provenance': 'VIOLA signal_noise.csv; each exposure checked against first',
        'recovered_from_existing_etc': resume_from is not None,
        'pwv_model': {'mean_mm': mean, 'sigma_mm': scatter, 'timescale_minutes':tau,
                      'seed':cfg.get_int('viola_pwv_seed',seed)}}
    (output_dir / 'viola_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return [output_dir / 'order_0.fits', output_dir / 'wave_0.fits']
