"""High-resolution POSEIDON/SYSREM cross-correlation for existing observations.

No observation or template files are modified. POSEIDON and tayph are optional
runtime dependencies: importing exoatmosim does not import them.
"""
from __future__ import annotations

import contextlib
import csv
import inspect
import io
import json
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def true_anomaly(jd, tc_jd, period_days, eccentricity, omega_rad):
    """Solve Kepler's equation from the inferior-conjunction epoch.

    Uses the edge-on f(Tc) + omega = pi/2 convention; this should be
    compared against the RVTraceEstimator orbit for non-edge-on systems.
    """
    jd = np.asarray(jd, float)
    e = float(eccentricity)
    if (jd.ndim != 1 or not np.isfinite(jd).all() or
            not np.isfinite([tc_jd, period_days, e, omega_rad]).all()
            or period_days <= 0 or not 0 <= e < 1):
        raise ValueError('Invalid ephemeris for true anomaly')
    f0 = np.pi / 2 - omega_rad
    E0 = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(f0 / 2),
                         np.sqrt(1 + e) * np.cos(f0 / 2))
    M0 = E0 - e * np.sin(E0)
    M = (M0 + 2 * np.pi * (jd - tc_jd) / period_days + np.pi) % (2*np.pi) - np.pi
    E = M.copy()
    for _ in range(80):
        delta = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= delta
        if np.max(np.abs(delta)) < 1e-12:
            break
    else:
        raise RuntimeError('Kepler equation did not converge')
    f = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    # Preserve a continuous anomaly across periastron branch cuts. The RV
    # equation uses sin/cos, so adding 2*pi does not change physical speeds.
    # The observation schedule supplies times in chronological order.
    return np.unwrap(f)


def observation_orbit(run, n_exposures):
    """Reproduce the synthetic orbit configuration and transit weight."""
    from astropy.time import Time
    from .config import Config
    from .independent_crires import staged_orbital_config
    from .rvtrace import RVTraceEstimator

    cfg = Config.from_file(run / 'config')
    with (run / 'diagnostics' / 'observation_schedule.csv').open(newline='') as fh:
        schedule = list(csv.DictReader(fh))
    if len(schedule) != n_exposures:
        raise ValueError('Schedule exposure count does not match FITS')
    center = .5 * (float(schedule[0]['start_jd_utc']) +
                    float(schedule[-1]['end_jd_utc']))
    mode = 'explicit' if cfg.get('Tobs') else 'max_altitude'
    staged = staged_orbital_config(cfg, SimpleNamespace(mode=mode, center_jd_utc=center))
    rows = [line.split() for line in (run / 'obs_times').read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    if len(rows) != n_exposures:
        raise ValueError('obs_times exposure count does not match FITS')
    jd = Time(np.array([float(row[0]) for row in rows]), format='mjd', scale='utc').jd
    berv = np.array([float(row[4]) for row in rows])
    if not np.isfinite(berv).all():
        raise ValueError('Non-finite BERV')
    tc = staged.get_float('Tc')
    eccentricity = staged.get_float('ecc', 0.)
    omega_deg = staged.get_float('omega', 90.)
    phi = true_anomaly(jd, tc, staged.get_float('P'), eccentricity,
                       np.deg2rad(omega_deg))
    with tempfile.TemporaryDirectory(prefix='exoatmosim_cc_') as td:
        stage = Path(td)
        (stage / 'config').write_text(staged.to_text())
        shutil.copyfile(run / 'obs_times', stage / 'obs_times')
        with contextlib.redirect_stdout(io.StringIO()):
            estimator = RVTraceEstimator(str(stage), cfg.get('observatory', 'paranal'))
            _ = estimator.RVstar
            _ = estimator.RVplanet
            transit = np.asarray(estimator.transit, float)
    if transit.shape != (n_exposures,) or not np.isfinite(transit).all():
        raise ValueError('Invalid transit light curve')
    depth = np.clip(1 - transit, 0, 1)
    weights = depth / np.max(depth) if np.max(depth) > 0 else np.zeros_like(depth)
    if not np.any(weights > 0):
        raise ValueError('No in-transit exposures')
    return phi, berv, weights, eccentricity, np.deg2rad(omega_deg)


def expected_planet_kp(config):
    """Planet RV semi-amplitude in km/s for a general eccentric Kepler orbit.

    Config masses are in solar and Jupiter masses, period in days,
    inclination in degrees. Includes both masses rather than M_p << M_star.
    """
    import astropy.units as u
    from astropy.constants import G
    from .config import Config

    cfg = Config.from_file(config)
    mstar = cfg.get_float('Ms') * u.Msun
    mplanet = cfg.get_float('Mp') * u.Mjup
    period = cfg.get_float('P') * u.day
    incl = np.deg2rad(cfg.get_float('inclination'))
    ecc = cfg.get_float('ecc', 0.)
    if (not np.isfinite([mstar.value, mplanet.value, period.value,
                         incl, ecc]).all() or mstar <= 0*u.kg or
            mplanet <= 0*u.kg or period <= 0*u.s or
            not 0 <= ecc < 1):
        raise ValueError('Invalid orbital inputs for expected planetary Kp')
    kp = ((2*np.pi*G/period)**(1/3) * mstar *
          np.sin(incl) / ((mstar + mplanet)**(2/3) *
                          np.sqrt(1 - ecc*ecc)))
    return float(kp.to_value(u.km/u.s))


def read_observation_berv(run):
    """Read the exposure-midpoint BERV (km/s) from ExoAtmoSim obs_times."""
    path = Path(run) / 'obs_times'
    rows = [line.split() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    if not rows or any(len(row) < 5 for row in rows):
        raise ValueError(f'Missing five-column obs_times BERV: {path}')
    values = np.array([float(row[4]) for row in rows], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f'Non-finite BERV in {path}')
    return values


def order_arrays(run, index, n_exposures):
    from astropy.io import fits
    with fits.open(run / f'wave_{index}.fits') as hdul:
        wave = np.asarray(hdul[0].data, float)
        medium = hdul[0].header.get('WAVEMED', '').upper()
    flux = np.asarray(fits.getdata(run / f'order_{index}.fits'), float)
    if not medium:
        # Older ExoAtmoSim writers did not include wavelength-medium headers.
        # Recover the convention from the immutable seed-specific config.
        from .config import Config
        import warnings
        cfg = Config.from_file(run / 'config')
        medium = 'AIR' if cfg.get_bool('air', True) else 'VACUUM'
        warnings.warn(
            f'Order {index}: wave FITS lacks WAVEMED; using {medium} from '
            f'{run / "config"}', RuntimeWarning)
    if medium not in ('AIR', 'VAC', 'VACUUM'):
        raise ValueError(f'Order {index}: unrecognized WAVEMED={medium!r}')
    if medium != 'AIR':
        raise ValueError(f'Order {index}: expected AIR wavelengths in nm; '
                         'vacuum input requires a separate conversion path')
    if wave.ndim == 2:
        if wave.shape[0] != n_exposures or not np.allclose(wave, wave[0]):
            raise ValueError('Exposure-dependent grids are not supported yet')
        wave = wave[0]
    if wave.ndim != 1 or flux.shape != (n_exposures, len(wave)):
        raise ValueError(f'Order {index}: incompatible wavelength/flux shapes')
    try:
        from tayph.operations import airtovac
    except ImportError as exc:
        raise ImportError('CC requires tayph for air-to-vacuum conversion') from exc
    vac = np.asarray(airtovac(wave), float)
    valid = np.isfinite(vac) & np.all(np.isfinite(flux) & (flux > 0), axis=0)
    if np.count_nonzero(valid) < 20:
        raise ValueError(f'Order {index}: fewer than 20 usable wavelengths ({np.count_nonzero(valid)}/{len(wave)})')
    return vac[valid] / 1000, flux[:, valid]


def _write_raw(path, wave_um, flux, phi, berv, weight):
    import h5py
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, 'w') as handle:
        for name, value in {
            'flux': flux[np.newaxis, :, :],
            'wl_grid': wave_um[np.newaxis, :],
            'phi': phi, 'V_bary': berv, 'transit_weight': weight,
        }.items():
            handle.create_dataset(name, data=value)


def _template(path):
    from astropy.io import fits
    with fits.open(path) as hdul:
        wave = np.asarray(hdul[1].data, float).ravel() / 1000
        depth = np.asarray(hdul[2].data, float).ravel()
    if (len(wave) < 3 or depth.shape != wave.shape or
            not np.isfinite(wave).all() or not np.isfinite(depth).all() or
            np.any(np.diff(wave) <= 0)):
        raise ValueError(f'Invalid species template {path}')
    # ExoAtmoSim writes species-template wavelengths in AIR nm, even though
    # the underlying POSEIDON model is computed on a vacuum-micron grid.
    # The observed-order wavelengths are already air->vacuum converted by
    # order_arrays(); make the templates use the very same conversion.
    from tayph.operations import airtovac
    vacuum_um = np.asarray(airtovac(wave * 1000.), float) / 1000.
    if not np.isfinite(vacuum_um).all() or np.any(np.diff(vacuum_um) <= 0):
        raise ValueError(f'Invalid vacuum-converted species template: {path}')
    return vacuum_um, depth - np.median(depth)


def cross_correlate_run(config, *, seed=42, sysrem=5, species=('O2','O3','H2O'),
                        prepare_only=False, cc_only=False, overwrite=False):
    """Prepare existing order FITS and/or compute POSEIDON CCF maps.

    The project-local `cc/` tree contains all products. Templates are read
    from the seed's existing `templates/` folder without regeneration.
    """
    from .workflow import _locations
    from astropy.io import fits
    from POSEIDON.high_res import (remove_outliers, fit_uncertainties,
                                  prepare_high_res_data, read_high_res_data,
                                  cross_correlate)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if prepare_only and cc_only:
        raise ValueError('--prepare-only and --cc-only are mutually exclusive')
    if not isinstance(sysrem, int) or sysrem < 0:
        raise ValueError('SYSREM iterations must be a nonnegative integer')
    if not species or any(not x or '/' in x or '\\' in x for x in species):
        raise ValueError('Provide one or more valid molecule names')
    _, run, _ = _locations(config, seed)
    if not run.is_dir():
        raise FileNotFoundError(f'No simulation run: {run}')
    out = run / 'cc'
    prepared = out / f'sysrem_{sysrem}' / 'prepared'
    order_ids = sorted(int(p.stem.split('_')[-1]) for p in run.glob('order_*.fits')
                       if p.stem.split('_')[-1].isdigit())
    if not order_ids:
        raise FileNotFoundError(f'No observed spectral orders under {run}')
    names = [f'order_{i}' for i in order_ids]
    if not cc_only:
        if prepared.exists() and not overwrite:
            raise FileExistsError(f'Prepared products already exist: {prepared}; use --overwrite or --cc-only')
        from astropy.io import fits
        n_exposures = np.asarray(fits.getdata(run / f'order_{order_ids[0]}.fits')).shape[0]
        phi, berv, weight, ecc, omega_rad = observation_orbit(run, n_exposures)
        usable_names = []
        skipped = {}
        for i, name in zip(order_ids, names):
            try:
                wave_um, flux = order_arrays(run, i, n_exposures)
            except ValueError as exc:
                if 'fewer than 20 usable wavelengths' not in str(exc):
                    raise
                skipped[name] = str(exc)
                print(f'[CC] Skipping {name}: {exc}')
                continue
            usable_names.append(name)
            wave_by_order = wave_um[np.newaxis, :]
            flux_by_order = flux[np.newaxis, :, :]
            _write_raw(prepared / name / 'data_raw.hdf5', wave_um, flux, phi, berv, weight)
            clean = remove_outliers(wave_by_order, flux_by_order)
            uncertainty = fit_uncertainties(clean, n_components=5)
            prepare_high_res_data(str(prepared), name, 'transmission', 'sysrem',
                                  clean, wave_by_order, phi, uncertainty,
                                  sysrem_niter=sysrem, transit_weight=weight)
        if not usable_names:
            raise ValueError('No spectral segments have at least 20 usable wavelengths; '
                             'ANDES computational segments are too short for cross-correlation')
        names = usable_names
        (prepared / 'preparation.json').write_text(json.dumps({
            'sysrem': sysrem, 'phase_kind': 'true_anomaly_radians',
            'eccentricity': ecc, 'omega_rad': omega_rad, 'orders': names,
            'skipped_short_orders': skipped,
        }, indent=2))
        if prepare_only:
            return prepared
    metadata_path = prepared / 'preparation.json'
    if not metadata_path.is_file():
        raise FileNotFoundError(f'Missing preparation metadata: {metadata_path}')
    metadata = json.loads(metadata_path.read_text())
    if cc_only:
        names = metadata.get('orders', [])
        if not names or any(name not in [f'order_{i}' for i in order_ids] for name in names):
            raise ValueError('Prepared order selection does not match observed FITS')
    if (metadata.get('sysrem') != sysrem or
            metadata.get('phase_kind') != 'true_anomaly_radians' or
            metadata.get('orders') != names):
        raise ValueError('Prepared data were produced with incompatible settings; rerun preparation')
    data = read_high_res_data(str(prepared), names=names)
    # POSEIDON's prepare/read functions do not necessarily retain V_bary from
    # the raw HDF5. Inject the authoritative observer-frame correction directly
    # into the per-order dict supplied to cross_correlate.
    berv = read_observation_berv(run)
    print(f'ExoAtmoSim obs_times BERV: {berv.min():.3f} to '
          f'{berv.max():.3f} km/s')
    for name in names:
        if 'phi' not in data[name]:
            raise ValueError(f'Prepared order {name} lacks phi')
        if len(np.asarray(data[name]['phi']).ravel()) != len(berv):
            raise ValueError(f'BERV exposure count differs from prepared {name}')
        old = data[name].get('V_bary')
        if old is None:
            print(f'{name}: prepared V_bary absent; restoring from obs_times')
        else:
            old = np.asarray(old, float)
            print(f'{name}: prepared V_bary range: '
                  f'{old.min():.3f} to {old.max():.3f} km/s')
        data[name]['V_bary'] = berv.copy()
    kp_expected = expected_planet_kp(run / 'config')
    print(f'Expected planetary Kp: {kp_expected:.3f} km/s')
    kp = np.arange(-300., 301., 1.)
    vsys = np.arange(-200., 201., 1.)
    rv = np.arange(-500., 501., 1.)
    # Support both keyword-driven CCF implementations and POSEIDON forks
    # that read orbital elements from the first 'dp' (run directory) argument.
    params = inspect.signature(cross_correlate).parameters
    print(f'POSEIDON cross_correlate: {inspect.getfile(cross_correlate)}')
    print(f'POSEIDON signature: {inspect.signature(cross_correlate)}')
    kwargs = {}
    if abs(metadata['eccentricity']) > 1e-12:
        if {'ecc', 'omega_rad'} <= set(params):
            kwargs = {'ecc': metadata['eccentricity'],
                      'omega_rad': metadata['omega_rad']}
        elif 'dp' in params:
            # The Mantis POSEIDON fork reads ecc/omega from dp/config
            # internally. Do not pass unsupported keyword arguments.
            try:
                source = inspect.getsource(cross_correlate)
            except (OSError, TypeError) as exc:
                raise RuntimeError(
                    'Cannot inspect the eccentric correction in POSEIDON; '
                    'verify the implementation before running') from exc
            if not all(token in source for token in ('ecc', 'omega_rad', 'np.cos')):
                raise RuntimeError(
                    'POSEIDON dp interface found, but no recognizable eccentric '
                    'true-anomaly correction in cross_correlate source')
            print('Using POSEIDON eccentric velocity correction from dp/config')
        else:
            raise RuntimeError(
                'Installed POSEIDON cross_correlate lacks support for eccentric '
                'velocity corrections (ecc/omega_rad or dp configuration)')
    with plt.rc_context({'figure.constrained_layout.use': False,
                         'figure.autolayout': False}):
        for molecule in species:
            template = run / 'templates' / f'{molecule}.fits'
            if not template.is_file():
                raise FileNotFoundError(f'Missing template: {template}')
            wave, depth = _template(template)
            ccf_per_order = []
            maps = []
            for name in names:
                map_order, ccf_order = cross_correlate(str(run), kp, vsys, rv,
                                                      wave, depth, data[name], **kwargs)
                maps.append(map_order)
                ccf_per_order.append(ccf_order)
            output = out / f'sysrem_{sysrem}' / molecule
            output.mkdir(parents=True, exist_ok=True)
            phase_rv = np.nanmean(np.asarray(ccf_per_order), axis=0)
            kpmap = np.nansum(np.asarray(maps), axis=0)
            for filename, array in [('phase_rv.fits', phase_rv),
                                    ('kp_vsys.fits', kpmap), ('rv_axis.fits', rv),
                                    ('kp_axis.fits', kp), ('vsys_axis.fits', vsys)]:
                fits.writeto(output / filename, array, overwrite=overwrite)
            phi = np.asarray(data[names[0]]['phi'], float)
            from .cc_plots import plot_phase_rv, plot_kp_vsys
            plot_phase_rv(phi, rv, phase_rv, output / 'phase_rv.png')
            plot_kp_vsys(kp, vsys, kpmap, output / 'kp_vsys.png',
                         molecule=molecule, kp_expected=kp_expected)
    return out / f'sysrem_{sysrem}'
