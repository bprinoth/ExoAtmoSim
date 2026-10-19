"""Independent constant-SNR CRIRES_K mock (experimental, no previous external mock implementation).

Run from project root after exporting the existing planetary FITS:
  python examples/run_independent_crires.py data/TOI-1842_b/config --snr 100

By default, retrieve a PHOENIX stellar spectrum using the existing provider.
VIOLA_J/H/K use the separate native PHOENIX + TelFit + ETC-noise backend.
An optional --stellar-npz path can override this with wavelength_nm and flux
arrays on the same AIR wavelength convention as the CRIRES_K order grid.
Tellurics use optional TelFit, evaluated at every exposure altitude.
"""
from __future__ import annotations

import argparse
import contextlib
import io
from pathlib import Path
import sys
import tempfile
from dataclasses import replace
import numpy as np

from exoatmosim.config import Config
from exoatmosim.observing import make_observation_schedule, save_schedule
from exoatmosim.observing.order_mock import simulate_orders, write_spectral_orders
from exoatmosim.observing.spectral_processing import observer_frame_velocities, broaden_spectrum
from exoatmosim.observing.obs_times import write_observation_times
from exoatmosim.observing.conditions import (observing_conditions, write_observing_conditions,
                                             read_observing_conditions, plot_observing_conditions)
from exoatmosim.instruments import instrument_resource
from exoatmosim.plotting import publication_style
from exoatmosim.stellar_spectra import PhoenixSpectrumProvider
from exoatmosim.objects import Star


def velocity_values_kms(velocity):
    """Return numeric km/s velocities from the legacy RVTraceEstimator.

    RVTraceEstimator.RVstar is returned by PyAstronomy as bare numeric km/s;
    RVplanet multiplies it by an Astropy dimensionless mass ratio, making it
    a dimensionless Quantity *even though its values are still km/s*.
    A physical speed Quantity is properly converted to km/s, while the
    legacy dimensionless quantity's numeric values are preserved.
    """
    try:
        import astropy.units as u
    except ImportError:
        if hasattr(velocity, "unit"):
            raise ImportError("Astropy is required for Quantity velocities")
        result = np.asarray(velocity, dtype=float)
        if result.ndim != 1:
            raise ValueError("RV estimator must return a finite 1-D velocity series")
        return result
    if isinstance(velocity, u.Quantity):
        if velocity.unit.is_equivalent(u.km / u.s):
            values = velocity.to_value(u.km / u.s)
        elif velocity.unit.is_equivalent(u.dimensionless_unscaled):
            values = velocity.to_value(u.dimensionless_unscaled)
        else:
            raise u.UnitConversionError(
                f'RV estimator returned an incompatible unit: {velocity.unit}')
    else:
        values = velocity
    result = np.asarray(values, dtype=float)
    if result.ndim != 1:
        raise ValueError('RV estimator must return a 1-D velocity series')
    return result


def validate_observation_inputs(rv_star, rv_planet, light_curve, cfg):
    """Explain invalid orbital inputs before entering expensive order synthesis.

    We deliberately do not replace unknown masses/RVs with zero: doing so would
    silently change the high-resolution planet Doppler trace.
    """
    for name, values in [('stellar RV', rv_star), ('planet RV', rv_planet),
                         ('transit light curve', light_curve)]:
        arr = np.asarray(values, dtype=float)
        bad = np.flatnonzero(~np.isfinite(arr))
        if bad.size:
            details = ''
            if name == 'planet RV':
                details = (f" Check Mp={cfg.get('Mp', '<missing>')} and "
                           f"Ms={cfg.get('Ms', '<missing>')}; the legacy estimator "
                           'scales RVstar by Ms/Mp.')
            elif name == 'stellar RV':
                details = f" Check K={cfg.get('K', '<missing>')} and orbital elements."
            raise ValueError(f'Non-finite {name} at {bad.size}/{arr.size} exposures '
                             f'(first index {int(bad[0])}, value {arr[bad[0]]}).{details}')
    curve = np.asarray(light_curve, float)
    if np.any((curve < -1e-8) | (curve > 1.0+1e-8)):
        raise ValueError('Transit light curve outside physical transmission range [0,1]. '
                         'Check occultation geometry in RVTraceEstimator.')


def staged_orbital_config(cfg, schedule):
    """Return a temporary config with an idealized transit at culmination.

    The original ephemeris is untouched. In max-altitude mode we intentionally
    interpret the geometric UTC JD as a synthetic Tc in the estimator's legacy
    JD convention; this is NOT a physical UTC -> BJD_TDB conversion.
    """
    staged = Config.from_text(cfg.to_text())
    if schedule.mode == 'max_altitude':
        staged.set('Tc', f'{schedule.center_jd_utc:.12f}')
    return staged


def run(config_path, stellar_npz, *, output=None, snr=None, seed=42, overwrite=False, resume_viola_from=None, model_path=None):
    from astropy.io import fits
    from exoatmosim.rvtrace import RVTraceEstimator

    config_path = Path(config_path).resolve()
    cfg = Config.from_file(config_path)
    instrument = str(cfg.get('instrument', 'CRIRES_K')).upper()
    if instrument not in ('CRIRES_K', 'ANDES', 'ANDES_K', 'IGRINS_HK', 'IGRINS2_HK', 'SPIROU', 'ESPRESSO', 'VIOLA_K', 'VIOLA_H', 'VIOLA_J'):
        raise ValueError(f'Unsupported instrument: {instrument}')
    if output is not None:
        out = Path(output).resolve()
    elif config_path.parent.name.startswith('seed_') and config_path.parent.parent.name == instrument:
        out = config_path.parent
    else:
        out = (config_path.parent / 'simulations' / instrument / f'seed_{seed:04d}').resolve()
    import hashlib
    import json
    cfg_digest = hashlib.sha256(config_path.read_bytes()).hexdigest()
    run_manifest_path = out / 'run_manifest.json'
    if run_manifest_path.exists():
        previous = json.loads(run_manifest_path.read_text())
        if (previous.get('config_sha256') != cfg_digest or
            previous.get('instrument') != instrument or
            previous.get('noise_seed') != seed):
            raise ValueError('Existing simulation run has different instrument, seed or config. '
                             'Select another --output directory or --seed; refusing to mix results.')
    out.mkdir(parents=True, exist_ok=True)
    snapshot_path = out / 'config'
    original = config_path.read_bytes()
    if snapshot_path.exists() and snapshot_path.read_bytes() != original and not overwrite:
        raise FileExistsError(f'Different config already saved for seed {seed}: {snapshot_path}')
    if snapshot_path != config_path:
        snapshot_path.write_bytes(original)
    diagnostic = out / 'diagnostics'
    diagnostic.mkdir(parents=True, exist_ok=True)
    schedule = make_observation_schedule(config_path)
    save_schedule(schedule, diagnostic)
    conditions = observing_conditions(
        schedule,
        humidity_percent=cfg.get_float('telfit_humidity_percent', 50.0),
        humidity_sigma_percent=cfg.get_float('telfit_humidity_sigma_percent', 0.0),
        humidity_timescale_minutes=cfg.get_float('telfit_humidity_timescale_minutes', 30.0),
        pressure_hpa=cfg.get_float('telfit_pressure_hpa', 600.0),
        temperature_k=cfg.get_float('telfit_temperature_k', 273.0),
        seed=cfg.get_int('telfit_seed', seed),
    )
    from exoatmosim.observing.conditions import sample_pwv
    pwv = sample_pwv(schedule, mean_mm=cfg.get_float('pwv',2.5),
                     sigma_mm=cfg.get_float('pwv_sigma',0.5),
                     timescale_minutes=cfg.get_float('pwv_timescale_minutes',30.),
                     seed=cfg.get_int('pwv_seed',seed))
    conditions['pwv_mm'] = pwv
    conditions['pwv_mm'].unit = 'mm'
    write_observing_conditions(out / 'observing_conditions.ecsv', conditions)
    conditions = read_observing_conditions(out / 'observing_conditions.ecsv', schedule)
    plot_observing_conditions(conditions, diagnostic / 'observing_conditions.png')
    if schedule.mode == 'max_altitude':
        print('[INFO] Idealized observation: the synthetic transit center is placed at '
              f'maximum altitude (JD UTC {schedule.center_jd_utc:.9f}); '
              'the historical Tc in the source config is unchanged.')
    # The integrated legacy estimator reads MJD obs_times and config in one folder.
    # Copy to a temporary directory to avoid overwriting the real planet files.
    with tempfile.TemporaryDirectory(prefix='exoatmosim_rv_') as tempdir:
        stage = Path(tempdir)
        (stage / 'config').write_text(staged_orbital_config(cfg, schedule).to_text())
        write_observation_times(stage / 'obs_times', schedule)
        # The legacy RV estimator prints many low-level progress messages.
        # Preserve them for debugging while keeping the interactive run readable.
        rvtrace_log = io.StringIO()
        try:
            with contextlib.redirect_stdout(rvtrace_log):
                estimator = RVTraceEstimator(str(stage), cfg.get('observatory', 'paranal'))
                rv_star = velocity_values_kms(estimator.RVstar)
                rv_planet = velocity_values_kms(estimator.RVplanet)
        finally:
            (diagnostic / 'rvtrace_backend.log').write_text(rvtrace_log.getvalue())
        # Save the original orbital-velocity diagnostic outside the temporary RV directory.
        import matplotlib.pyplot as plt
        with publication_style():
            estimator.plot_doppler_trace(RF='system', save=True)
            plt.close('all')
        trace_file = stage / 'RVTraceEstimator_system.png'
        if not trace_file.is_file():
            raise RuntimeError('RVTraceEstimator did not produce its system-frame plot')
        import shutil
        shutil.copy2(trace_file, diagnostic / trace_file.name)
        print(f'RV trace plot: {diagnostic / trace_file.name}')
        print('  Orbital solution and BERV calculated (details: diagnostics/rvtrace_backend.log)')
        # Fail with the specific quantity instead of an ambiguous order-mock error.
        validate_observation_inputs(rv_star, rv_planet, estimator.transit, cfg)
        # Legacy estimator returns a light curve near 1; use relative occultation.
        transit = np.asarray(estimator.transit, float)
        depth = np.clip(1.0 - transit, 0.0, 1.0)
        weight = depth / np.max(depth) if np.max(depth) > 0 else np.zeros_like(depth)
        # The estimator stores barycentric correction separately, and its RVstar
        # Kepler model uses c0=0. Follow frozen previous mock framework observer-frame shifts.
        rv_star_orbit, rv_planet_orbit = rv_star.copy(), rv_planet.copy()
        berv = velocity_values_kms(estimator.berv)
        rv_star, rv_planet = observer_frame_velocities(
            rv_star_orbit, rv_planet_orbit, cfg.get_float('vsys', 0.), berv)
        print(f'RV ranges (observer frame): star {rv_star.min():.3f} to {rv_star.max():.3f} km/s, '
              f'planet {rv_planet.min():.3f} to {rv_planet.max():.3f} km/s')
        print(f'Transit exposures: {np.count_nonzero(weight > 0)} / {len(weight)}')
        if not np.any(weight > 0):
            print('[WARNING] No in-transit exposures despite the selected center. '
                  'Check transit geometry, inclination, and T14.')
    if instrument in ('VIOLA_J', 'VIOLA_H', 'VIOLA_K'):
        if stellar_npz is not None or resume_viola_from is not None or snr is not None:
            raise ValueError('Native VIOLA uses PHOENIX and ETC noise; no stellar-npz, resume-viola-from or --snr')
        from exoatmosim.observing.viola_native import make_native_viola_observation
        model_file = (Path(model_path).resolve() if model_path is not None else
                      out / f'{config_path.parent.name}.fits')
        if not model_file.is_file():
            raise FileNotFoundError(f'Planet FITS not found: {model_file}; run spectrum first')
        with fits.open(model_file) as hdul:
            pw = np.asarray(hdul[1].data, float)
            pf = np.asarray(hdul[2].data, float)
        written = make_native_viola_observation(schedule, rv_star, rv_planet, weight,
                  cfg, config_path, pw, pf, out, seed=seed, overwrite=overwrite)
        write_observation_times(out / 'obs_times', schedule)
        run_manifest_path.write_text(json.dumps({
            'instrument':instrument,'noise_seed':seed,'config_sha256':cfg_digest,
            'config_snapshot':'config','n_exposures':schedule.n_exposures,
            'backend':'VIOLA_NATIVE_PHOENIX_HDF5_ETC_V0.5'},indent=2)+'\n')
        return written
    stellar_source = 'PHOENIX'
    if cfg.get('stellar_source', 'PHOENIX').upper() != 'PHOENIX':
        raise ValueError('Only PHOENIX stellar_source is supported')
    if stellar_npz is None:
        phoenix = PhoenixSpectrumProvider.from_config(cfg, config_path=config_path)
        # Provider's normal wavelength crop ends at the instrument boundaries.
        # Allow for Doppler motion and the finite convolution kernel BEFORE retrieval.
        if phoenix.wavelength_min_um is not None:
            phoenix = replace(phoenix,
                              wavelength_min_um=max(0.001, phoenix.wavelength_min_um - 0.010),
                              wavelength_max_um=phoenix.wavelength_max_um + 0.010)
        stellar = phoenix.get_spectrum(Star.from_config(cfg))
        if stellar.wavelength_medium != 'air':
            raise ValueError('Instrument order grid requires an AIR-wavelength stellar spectrum; set air True in config')
        sw = stellar.wavelength_nm
        sf = stellar.flux_erg_s_cm2_nm
        with publication_style():
            stellar.plot(diagnostic / 'stellar_phoenix.png')
        print(f'Loaded PHOENIX stellar spectrum: {len(sw)} pixels ({sw[0]:.2f}–{sw[-1]:.2f} nm air)')
    else:
        with np.load(stellar_npz, allow_pickle=False) as star:
            sw = np.asarray(star['wavelength_nm'], float)
            sf = np.asarray(star['flux'], float)
        print(f'Loaded supplied stellar spectrum: {stellar_npz}')
    model_file = (Path(model_path).resolve() if model_path is not None else
                  out / f'{config_path.parent.name}.fits')
    if not model_file.is_file():
        raise FileNotFoundError(f'Seed-specific planetary FITS not found: {model_file}. '
                                'Run the full pipeline to compute it first.')
    with fits.open(model_file) as hdul:
        pw = np.asarray(hdul[1].data, float)
        pf = np.asarray(hdul[2].data, float)
    from exoatmosim.instruments import load_instrument_orders
    orders = load_instrument_orders(instrument)
    if instrument in ('ANDES', 'ANDES_K'):
        import warnings
        warnings.warn(f'{instrument} is EXPERIMENTAL: supplied wavelength samples in computational segments, not calibrated physical orders. Constant SNR is not ETC SNR.', RuntimeWarning)
    if not all(np.all(np.diff(w) > 0) for w in orders):
        bad = [i for i,w in enumerate(orders) if not np.all(np.diff(w)>0)]
        raise ValueError(f'{instrument} calibration has non-monotonic wavelengths in order(s) '
                         f'{bad}; order data must be validated before simulation')
    # Crop PHOENIX before applying rotational broadening; the original model
    # can span millions of pixels. Keep enough margin for Doppler + broadening.
    margin_nm = max(3., np.max(np.abs(rv_star)) / 299792.458 * max(w.max() for w in orders) + 1.)
    selected = (sw >= min(w.min() for w in orders) - margin_nm) & (sw <= max(w.max() for w in orders) + margin_nm)
    if np.count_nonzero(selected) < 8:
        raise ValueError('Stellar spectrum does not cover instrument orders')
    sw = sw[selected]
    sf = sf[selected]
    resolution = cfg.get_float('resolution')
    vsini = cfg.get_float('vsini', 0.)
    print(f'Stellar rotational broadening: vsini={vsini:.3f} km/s; instrument R={resolution:.0f}')
    sw, sf = broaden_spectrum(sw, sf, resolution=resolution, vsini_kms=vsini)
    # Broaden the POSEIDON transmission to instrument resolution (without planetary rotation).
    pw, pf = broaden_spectrum(pw, pf, resolution=resolution)
    use_snr = float(snr) if snr is not None else cfg.get_float('SNR', 100.0)
    extra = {}
    if cfg.get_bool('tellurics', False):
        from exoatmosim.observing.telluric_grid import grid_transmission, resolve_grid_file
        from exoatmosim.observing.tellurics import plot_tellurics
        # Evaluate optical-depth interpolation at the instrument's native air-nm
        # sample wavelengths; no coarse 2000-point telluric sampling.
        grid = np.unique(np.concatenate(orders))
        tw, transmission = grid_transmission(
            schedule, grid, pwv_mm=np.asarray(conditions['pwv_mm'],float),
            grid_file=resolve_grid_file(cfg.get('observatory','paranal'), config_path=config_path))
        plot_tellurics(tw, transmission, schedule, diagnostic / 'telluric_transmission.png')
        extra = dict(telluric_wavelength_nm=tw, telluric_transmission=transmission)
    results = simulate_orders(orders, sw, sf, pw, pf, rv_star, rv_planet,
                              weight, snr=use_snr, seed=seed,
                              stellar_observer_frame=False, **extra)
    written = write_spectral_orders(results, out, overwrite=overwrite)
    # Mandatory quick-look diagnostic: verify that an actual order was synthesized.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    wave0, flux0, _ = results[0]
    with publication_style():
        fig, ax = plt.subplots(figsize=(10, 4))
        for j in sorted(set([0, len(flux0)//2, len(flux0)-1])):
            ax.plot(wave0[j], flux0[j], lw=0.9, label=f'Exposure {j}')
        ax.set(xlabel='Wavelength [nm]', ylabel='Flux', title=f'Independent {instrument} order 0')
        ax.legend()
        fig.tight_layout()
        fig.savefig(diagnostic / 'order_0_first_middle_last.png', dpi=180)
        plt.close(fig)
    write_observation_times(out / 'obs_times', schedule)
    np.savetxt(diagnostic / 'rv_kms.csv',
               np.column_stack((schedule.midpoint_jd_utc, rv_star_orbit, rv_planet_orbit,
                                np.full(len(rv_star), cfg.get_float('vsys', 0.)), berv,
                                rv_star, rv_planet, weight)),
               header='mid_jd_utc,rvstar_orbit_kms,rvplanet_orbit_kms,vsys_kms,berv_kms,rvstar_observer_kms,rvplanet_observer_kms,transit_weight',
               delimiter=',', comments='')
    run_manifest_path.write_text(json.dumps({
        'instrument': instrument, 'noise_seed': seed, 'config_sha256': cfg_digest,
        'config_snapshot': 'config', 'n_exposures': schedule.n_exposures, 'constant_snr': use_snr,
        'observer_frame_stellar_input': False,
    }, indent=2) + '\n')
    print(f'Stellar source: {stellar_source}')
    print(f'Wrote {len(results)} {instrument} orders for {schedule.n_exposures} exposures, SNR={use_snr:g}')
    print(f'Observation window: {schedule.window_hours:.3f} h ({schedule.mode})')
    print(f'Output: {out} ({len(written)} verified FITS files)')
    return written


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    parser.add_argument('--stellar-npz', type=Path, help='Optional wavelength_nm and flux arrays; default: PHOENIX')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--snr', type=float)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--model-path', type=Path, help='Explicit planetary model FITS (otherwise seed directory)')
    args = parser.parse_args()
    run(args.config, args.stellar_npz, output=args.output, snr=args.snr,
        seed=args.seed, overwrite=args.overwrite,
        model_path=args.model_path)
