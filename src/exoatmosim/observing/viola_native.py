"""Experimental independent VIOLA J/H/K backend: PHOENIX + HDF5 tellurics + ETC-inspired noise.

Output samples are supplied VIOLA wavelength positions; each reports a
*resolution-element* electron estimate (lambda/R), even if oversampled.
Adjacent samples are therefore not statistically independent; noise draws
are independent as a first approximation, not a detector covariance model.
"""
from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from ..instruments import load_instrument_orders
from ..objects import Star
from ..stellar_spectra import PhoenixSpectrumProvider
from .spectral_processing import broaden_spectrum
from .telluric_utils import vacuum_to_air_nm

C = 299792458.0
H = 6.62607015e-34
K_B = 1.380649e-23
FNU0_VEGA = {'J': 1.594e-23, 'H': 1.024e-23, 'K': 0.6668e-23}
THROUGHPUT = (0.95, 0.95, 0.90, 0.96, 0.90, 0.90, 0.80, 0.90, 0.95)


def viola_etc_variance(signal_res_e, sky_res_e, thermal_res_e, dark_res_e,
                       *, nod=True, read_noise_e=20., n_reads=1, n_detector_pixels=4):
    """Replicate VIOLA ETC runner.py variance for a resolution element."""
    s, sky, thermal, dark = (np.asarray(a, float) for a in
                             (signal_res_e, sky_res_e, thermal_res_e, dark_res_e))
    if any(np.any(~np.isfinite(a)) or np.any(a < 0) for a in (s, sky, thermal, dark)):
        raise ValueError('VIOLA electron terms must be finite and nonnegative')
    if read_noise_e < 0 or n_reads < 1 or n_detector_pixels < 1:
        raise ValueError('Invalid VIOLA detector readout parameters')
    return s + (2. if nod else 1.)*(sky + thermal + dark) + read_noise_e**2*n_reads*n_detector_pixels


def calibrate_phoenix_flux(wave_air_nm, phoenix_surface_flux, band_magnitude, band,
                           *, magnitude_system='Vega'):
    """Anchor PHOENIX *shape* to ETC-style Fnu at band reference wavelength.

    This is a monochromatic approximation, not synthetic bandpass photometry.
    """
    wave = np.asarray(wave_air_nm, float)
    flux = np.asarray(phoenix_surface_flux, float)
    if wave.ndim != 1 or flux.shape != wave.shape or np.any(flux < 0) or not np.all(np.isfinite(flux)):
        raise ValueError('Invalid PHOENIX flux')
    reference_nm = 0.5*(wave[0]+wave[-1])
    reference_flux = float(np.interp(reference_nm, wave, flux))
    if reference_flux <= 0 or not np.isfinite(band_magnitude):
        raise ValueError('PHOENIX reference flux and band magnitude must be valid')
    system = magnitude_system.upper()
    if system == 'VEGA':
        fnu0 = FNU0_VEGA[band]
    elif system == 'AB':
        fnu0 = 3631e-26
    else:
        raise ValueError('VIOLA magnitude system must be Vega or AB')
    fnu = fnu0*10**(-0.4*band_magnitude)
    flambda_w_m2_um = C/(reference_nm*1e-9)**2 * fnu * 1e-6
    return flux*(flambda_w_m2_um/reference_flux)


def sky_emission_on_air_grid(path, wavelength_air_nm):
    """Load VIOLA ETC SkyCalc emission components, excluding flux_tie.

    SkyCalc FITS lam is VACUUM nm; interpolate sky brightness on AIR nm.
    """
    from astropy.io import fits
    with fits.open(path) as hdul:
        table = hdul[1].data
        cols = ('flux_zl', 'flux_ael', 'flux_tme', 'flux_sml', 'flux_ssl', 'flux_arc')
        missing = [name for name in ('lam', *cols) if name not in table.names]
        if missing:
            raise ValueError(f'SkyCalc file missing {missing}')
        wave = vacuum_to_air_nm(np.asarray(table['lam'], float))
        emission = sum(np.asarray(table[col], float) for col in cols)
    target = np.asarray(wavelength_air_nm, float)
    if wave[0] > target[0] or wave[-1] < target[-1]:
        raise ValueError('SkyCalc emission grid does not cover VIOLA band')
    if np.any(emission < 0) or not np.isfinite(emission).all():
        raise ValueError('Invalid SkyCalc emission')
    return np.interp(target, wave, emission)


def electron_budget(wave_nm, flux_w_m2_um, sky_phi, cfg, *, exposure_s):
    """VIOLA v0.5 equations on resolution-element centers (not per sample)."""
    wave = np.asarray(wave_nm, float)
    flux = np.asarray(flux_w_m2_um, float)
    sky = np.asarray(sky_phi, float)
    if wave.ndim != 1 or wave.shape != flux.shape or wave.shape != sky.shape:
        raise ValueError('VIOLA electron budget arrays must match')
    if np.any(flux < 0) or np.any(sky < 0) or not np.isfinite(flux).all() or not np.isfinite(sky).all():
        raise ValueError('Negative or nonfinite VIOLA photon flux')
    R = cfg.get_float('resolution', 300000.)
    diameter = 8.0
    seeing = cfg.get_float('viola_seeing_arcsec', 1.)
    slit = cfg.get_float('viola_slit_width_arcsec', 1.6)
    tau_opt = cfg.get_float('viola_optics_throughput', float(np.prod(THROUGHPUT)))
    nx = cfg.get_int('viola_nx', 2)
    ny = cfg.get_int('viola_ny', 2)
    spatial_scale = cfg.get_float('viola_spatial_scale_arcsec', 0.8)
    temperature = cfg.get_float('viola_ambient_temperature_k', 275.)
    emissivity = cfg.get_float('viola_thermal_emissivity', 0.10)
    dark = cfg.get_float('viola_dark_rate_e_s_pixel', 0.01)
    read = cfg.get_float('viola_read_noise_e_pixel', 20.)
    nod = cfg.get_bool('viola_nod_subtraction', True)
    if (min(R, diameter, seeing, slit, spatial_scale, temperature, exposure_s) <= 0
            or nx < 1 or ny < 1 or not 0 <= tau_opt <= 1 or not 0 <= emissivity <= 1
            or dark < 0 or read < 0):
        raise ValueError('Invalid VIOLA ETC instrument parameters')
    import math
    tau_slit = math.erf(math.sqrt(math.log(2.))*slit/seeing)
    area = np.pi*(diameter/2.)**2
    wave_um = wave*1e-3
    width_um = wave_um/R
    photon_energy = H*C/(wave*1e-9)
    star = flux/photon_energy*width_um*area*exposure_s*tau_opt*tau_slit
    omega_arcsec2 = slit*spatial_scale*ny
    sky_e = area*omega_arcsec2*tau_opt*width_um*sky*exposure_s
    omega_sr = omega_arcsec2*(np.pi/(180*3600))**2
    nu = C/(wave*1e-9)
    bnu = (2*H*nu**3/C**2)/np.expm1(H*nu/(K_B*temperature))
    thermal_e = area*omega_sr/(H*R)*emissivity*bnu*exposure_s
    dark_e = np.full_like(star, dark*exposure_s*nx*ny)
    variance = viola_etc_variance(star, sky_e, thermal_e, dark_e, nod=nod,
                                  read_noise_e=read, n_reads=1,
                                  n_detector_pixels=nx*ny)
    return star, variance, {'sky_e':sky_e, 'thermal_e':thermal_e, 'dark_e':dark_e,
                             'tau_slit':tau_slit, 'tau_optics':tau_opt}


def sample_viola_pwv(schedule, *, mean_mm=2.5, sigma_mm=0.5,
                     timescale_minutes=30., seed=42):
    """Stationary exponentially correlated PWV process (mm), per exposure."""
    times = np.asarray(schedule.midpoint_jd_utc, dtype=float)
    if (times.ndim != 1 or len(times) == 0 or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0) or not np.isfinite([mean_mm, sigma_mm, timescale_minutes]).all()
            or mean_mm <= 0 or sigma_mm < 0 or timescale_minutes <= 0):
        raise ValueError('Invalid VIOLA PWV schedule or parameters')
    from .viola_plots import plot_stellar_spectrum_comparison
    stellar_diagnostic = plot_stellar_spectrum_comparison(
        output_dir, star.wavelength_nm, original_stellar_flux,
        star_wave, star_shape, vsini)
    print(f'Stellar spectrum comparison: {stellar_diagnostic}')
    # Remove legacy duplicate only when regenerating an existing observation.
    if overwrite:
        for old_name in ('planet_transmission.png', 'viola_spectral_comparison.png'):
            old_plot = output_dir / 'diagnostics' / old_name
            if old_plot.exists():
                old_plot.unlink()
    rng = np.random.default_rng(seed)
    z = np.empty(len(times))
    z[0] = rng.normal() if sigma_mm else 0.
    for i in range(1, len(times)):
        rho = np.exp(-(times[i]-times[i-1])*1440./timescale_minutes)
        z[i] = rho*z[i-1] + np.sqrt(max(0., 1.-rho*rho))*rng.normal() if sigma_mm else 0.
    return np.maximum(0.001, mean_mm + sigma_mm*z)


def select_nearest_skycalc(candidates, altitudes_deg, pwv_mm):
    """Nearest available model per exposure; both axes separately in native units."""
    altitudes = np.asarray(altitudes_deg, float)
    pwv = np.asarray(pwv_mm, float)
    if (not candidates or altitudes.ndim != 1 or pwv.shape != altitudes.shape
            or not np.isfinite(altitudes).all() or not np.isfinite(pwv).all()
            or np.any((altitudes <= 0) | (altitudes > 90)) or np.any(pwv <= 0)):
        raise ValueError('Invalid VIOLA SkyCalc candidates, altitude or PWV')
    # The published VIOLA grid is rectangular. Pick closest altitude and PWV
    # nodes independently, then select the corresponding physical model.
    alt_nodes = sorted({a for _, a, _ in candidates})
    pwv_nodes = sorted({w for _, _, w in candidates})
    lookup = {(a, w): path for path, a, w in candidates}
    if len(lookup) != len(candidates):
        raise ValueError('Duplicate SkyCalc grid coordinates')
    selected = []
    for alt, water in zip(altitudes, pwv):
        a = min(alt_nodes, key=lambda x: abs(x-alt))
        w = min(pwv_nodes, key=lambda x: abs(x-water))
        if (a, w) not in lookup:
            raise ValueError(f'Missing SkyCalc altitude={a}, PWV={w} grid cell')
        selected.append((lookup[a,w], a, w))
    return selected


def resolve_viola_sky_dir(cfg, config_path):
    """Prefer explicit science config, else private site_paths.conf; never bundle FITS."""
    value = cfg.get('viola_sky_dir')
    if value:
        path = Path(value).expanduser()
        return path if path.is_absolute() else Path(config_path).resolve().parent/path
    from ..backend_paths import discover_site_config
    from ..config import Config
    site_path = discover_site_config()
    if site_path is None:
        return None
    site_cfg = Config.from_file(site_path)
    value = site_cfg.get('viola_skycalc_dir')
    if not value:
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else site_path.parent/path


def make_native_viola_observation(schedule, rv_star, rv_planet, transit_weight,
                                  cfg, config_path, planet_wave_nm, planet_factor,
                                  output_dir, *, seed=42, overwrite=False,
                                  phoenix_provider=None, telluric_provider=None):
    """Generate native-grid PHOENIX+HDF5 telluric VIOLA spectra and ETC electron noise."""
    from astropy.io import fits
    instrument = str(cfg.get('instrument')).upper()
    if instrument not in ('VIOLA_J', 'VIOLA_H', 'VIOLA_K'):
        raise ValueError('Native VIOLA supports J/H/K')
    band = instrument[-1]
    grid = load_instrument_orders(instrument)[0]
    if not cfg.get_bool('air', True):
        raise ValueError('VIOLA native requires air True for PHOENIX')
    resolution = cfg.get_float('resolution', 300000.)
    vsini = cfg.get_float('vsini', 0.)
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError('VIOLA instrumental resolution must be positive and finite')
    if not np.isfinite(vsini) or vsini < 0:
        raise ValueError('Stellar vsini must be finite and nonnegative')
    magnitude = cfg.get_float(f'{band}mag')
    n_exp = len(schedule.altitude_deg)
    rv_star = np.asarray(rv_star, float)
    rv_planet = np.asarray(rv_planet, float)
    weight = np.asarray(transit_weight, float)
    pw = np.asarray(planet_wave_nm, float)
    pf = np.asarray(planet_factor, float)
    if (rv_star.shape != (n_exp,) or rv_planet.shape != (n_exp,) or weight.shape != (n_exp,)
            or not np.isfinite(np.r_[rv_star,rv_planet,weight]).all()
            or np.any((weight < 0) | (weight > 1)) or np.any(np.diff(pw) <= 0)
            or pf.shape != pw.shape or not np.isfinite(pf).all()):
        raise ValueError('Invalid VIOLA orbital arrays')
    # Dense log wavelength grid resolves the configured LSF and narrow tellurics.
    margin = 5.0 + 2.*max(np.max(np.abs(rv_star)), np.max(np.abs(rv_planet)))*grid[-1]/299792.458
    r_sample = cfg.get_float('viola_internal_sampling_r', 900000.)
    if not np.isfinite(r_sample) or r_sample < max(900000., 3.*resolution):
        raise ValueError('VIOLA internal sampling must be finite and >= max(900000, 3*R)')
    n_hi = int(np.ceil(np.log((grid[-1]+margin)/(grid[0]-margin))*r_sample))+1
    high = np.geomspace(grid[0]-margin, grid[-1]+margin, n_hi)
    provider = phoenix_provider or PhoenixSpectrumProvider.from_config(cfg, config_path=config_path)
    provider = replace(provider, wavelength_min_um=(high[0]-1.)/1000.,
                       wavelength_max_um=(high[-1]+1.)/1000.)
    star = provider.get_spectrum(Star.from_config(cfg))
    if star.wavelength_medium != 'air':
        raise ValueError('PHOENIX must be in AIR nm')
    if high[0] < star.wavelength_nm[0] or high[-1] > star.wavelength_nm[-1]:
        raise ValueError('PHOENIX does not cover VIOLA high-resolution model grid')
    star_shape = calibrate_phoenix_flux(star.wavelength_nm, star.flux_erg_s_cm2_nm,
                                       magnitude, band,
                                       magnitude_system='Vega')
    # Ensure normalization at VIOLA band midpoint, not extended provider bounds.
    center = 0.5*(grid[0]+grid[-1])
    fnu0 = FNU0_VEGA[band]
    reference = C/(center*1e-9)**2 * fnu0*10**(-0.4*magnitude)*1e-6
    star_shape *= reference/np.interp(center, star.wavelength_nm, star_shape)
    original_stellar_flux = star_shape.copy()
    # Apply stellar rotation ONCE, in the stellar rest frame, before Doppler
    # shifting and before multiplying planet/telluric transmission. Never
    # convolve the star separately with the instrument LSF here.
    if vsini > 0:
        star_wave, star_shape = broaden_spectrum(
            star.wavelength_nm, star_shape, resolution=None, vsini_kms=vsini)
    else:
        star_wave = star.wavelength_nm
    # PWV is sampled once per exposure and used for the SkyCalc background.
    # The external HDF5 library was generated with its own atmospheric profile;
    # PWV (mm) is NOT interchangeable with RH (%).
    from .conditions import sample_pwv
    pwv_per_exp = sample_pwv(schedule,mean_mm=cfg.get_float('pwv',2.5),
                            sigma_mm=cfg.get_float('pwv_sigma',0.5),
                            timescale_minutes=cfg.get_float('pwv_timescale_minutes',30.),
                            seed=cfg.get_int('pwv_seed',seed))
    # VIOLA sky emission is deliberately out of scope for this release.
    # The zero-sky approximation must remain explicit in science products.
    sky_nodes = None
    sky_by_exp = [np.zeros_like(grid)] * n_exp
    sky_source = 'ZERO SKY (UNMODELLED; SNR IS AN OPTIMISTIC UPPER LIMIT)'
    if cfg.get_bool('tellurics', False):
        if telluric_provider is None:
            from .telluric_grid import grid_transmission, resolve_grid_file
            tw, transmission = grid_transmission(
                schedule,high,pwv_mm=pwv_per_exp,
                grid_file=resolve_grid_file(cfg.get('observatory','paranal'), config_path=config_path))
        else:
            tw, transmission = telluric_provider(schedule,high)
        if tw.shape != high.shape or not np.allclose(tw, high, atol=1e-8, rtol=0):
            raise ValueError('Telluric model grid differs from high-resolution grid')
        if transmission.shape != (n_exp,len(high)) or np.any((transmission < 0)|(transmission > 1)):
            raise ValueError('Invalid telluric transmission')
    else:
        transmission = np.ones((n_exp,len(high)))
    # Only samples within the instrumental LSF footprint contribute to the
    # delivered spectrum. The 5 nm high-grid padding is for the star/tellurics,
    # not a requirement for the planetary FITS.
    kernel_nm = 6.*grid[-1]/resolution
    required_low = min((grid[0]-kernel_nm)/(1.+v/299792.458) for v in rv_planet)
    required_high = max((grid[-1]+kernel_nm)/(1.+v/299792.458) for v in rv_planet)
    if pw[0] > required_low or pw[-1] < required_high:
        raise ValueError(f'Planetary FITS covers {pw[0]:.3f}-{pw[-1]:.3f} nm, '
                         f'but VIOLA {band} needs at least {required_low:.3f}-{required_high:.3f} nm. '
                         'Re-run exoatmosim spectrum with instrument VIOLA_'+band)
    if np.any(np.diff(star.wavelength_nm) <= 0):
        raise ValueError('PHOENIX wavelength grid not monotonic')
    output_dir = Path(output_dir).resolve()
    planned = [output_dir/name for name in ('wave_0.fits','order_0.fits','stellar_flux_0.fits',
               'model_flux_0.fits','noise_0.fits','snr_res_0.fits','viola_native_manifest.json')]
    if not overwrite and any(p.exists() for p in planned):
        raise FileExistsError('VIOLA native outputs exist; use --overwrite')
    rng = np.random.default_rng(seed)
    stars, models, sigmas, draws, snrs = [], [], [], [], []
    for i in range(n_exp):
        rest_star = high/(1.+rv_star[i]/299792.458)
        rest_planet = high/(1.+rv_planet[i]/299792.458)
        if rest_star[0] < star_wave[0] or rest_star[-1] > star_wave[-1]:
            raise ValueError('PHOENIX RV shifted beyond model coverage')
        stellar_high = np.interp(rest_star, star_wave, star_shape)
        planet_high = 1.+weight[i]*(np.interp(rest_planet,pw,pf,left=1.,right=1.)-1.)
        # Convolve the *combined* star x planet x atmosphere once.
        flux_star_high = stellar_high*transmission[i]
        flux_model_high = flux_star_high*planet_high
        lsf_star_wave, star_blurred = broaden_spectrum(high,flux_star_high,resolution=resolution)
        lsf_model_wave, model_blurred = broaden_spectrum(high,flux_model_high,resolution=resolution)
        # broaden_spectrum resamples to a log-wavelength grid. Use its returned
        # wavelength axis, never the original `high` axis: the lengths can differ.
        star_grid = np.interp(grid, lsf_star_wave, star_blurred)
        model_grid = np.interp(grid, lsf_model_wave, model_blurred)
        star_e, var, _ = electron_budget(grid,star_grid,sky_by_exp[i],cfg,
                                         exposure_s=cfg.get_float('exposure_time'))
        model_e, _, _ = electron_budget(grid,model_grid,sky_by_exp[i],cfg,
                                        exposure_s=cfg.get_float('exposure_time'))
        # Model-dependent stellar shot variance, same nonstellar variance.
        model_var = np.maximum(0.,var-star_e)+model_e
        sigma = np.sqrt(model_var)
        stars.append(star_e); models.append(model_e); sigmas.append(sigma)
        draws.append(model_e+rng.normal(size=len(grid))*sigma)
        snrs.append(np.divide(star_e,np.sqrt(var),out=np.zeros_like(star_e),where=var>0))
    data = {'wave_0.fits':np.broadcast_to(grid,(n_exp,len(grid))),
            'order_0.fits':np.asarray(draws), 'stellar_flux_0.fits':np.asarray(stars),
            'model_flux_0.fits':np.asarray(models), 'noise_0.fits':np.asarray(sigmas),
            'snr_res_0.fits':np.asarray(snrs)}
    output_dir.mkdir(parents=True,exist_ok=True)
    # One row per exposure, even when there is no sky model (for reproducibility).
    with (output_dir/'viola_sky_conditions.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['exposure_id','altitude_deg','sampled_pwv_mm',
                         'selected_sky_file','selected_sky_altitude_deg','selected_sky_pwv_mm'])
        for i, (alt, water) in enumerate(zip(schedule.altitude_deg,pwv_per_exp)):
            name, a, w = sky_nodes[i] if sky_nodes is not None else ('', '', '')
            writer.writerow([i, float(alt), float(water), name, a, w])
    for filename, array in data.items():
        hdr=fits.Header({'INSTRUME':instrument,'WAVEMED':'AIR',
                         'BUNIT':'nm' if filename.startswith('wave') else
                                 'dimensionless' if filename.startswith('snr') else 'electron/resel'})
        fits.writeto(output_dir/filename,array,header=hdr,overwrite=overwrite)
    # No redundant VIOLA spectral comparison; POSEIDON already has its plot.
    metadata={'backend':'VIOLA_NATIVE_PHOENIX_HDF5_ETC_V0.5', 'instrument':instrument,
              'telescope_diameter_m':8.0,
              'resolving_power':resolution,'stellar_vsini_kms':vsini,'magnitude':magnitude,'magnitude_system':'Vega',
              'sky_source':sky_source,
              'pwv_model':'correlated Gaussian PWV (mm), interpolated external HDF5 telluric grid',
              'pwv_calibration_source':'external_hdf5_grid_metadata',
              'pwv_note':'External HDF5 grids may use an approximate PWV-to-RH mapping; scientific calibration pending',
              'sky_conditions_csv':'viola_sky_conditions.csv',
              'tellurics':cfg.get_bool('tellurics',False),
              'sampling':'native supplied centers, resolution-element photon counts; overlapping samples correlated in reality',
              'noise':'independent Gaussian samples are approximation; ETC nod/read variance conventions',
              'stellar_normalization':'monochromatic PHOENIX anchor; synthetic filter photometry pending',
              'n_exposures':n_exp,'n_samples':len(grid)}
    (output_dir/'viola_native_manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return [output_dir/'order_0.fits',output_dir/'wave_0.fits']
