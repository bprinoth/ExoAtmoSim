"""Reproducible exposure-level atmospheric conditions, stored as Astropy ECSV.

The first five columns of ``obs_times`` remain unchanged. This table is the
diagnostic record of observing conditions, and is also produced when
``tellurics False`` for a fully reproducible observing plan.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np


def observing_conditions(schedule, *, humidity_percent=50.0,
                         humidity_sigma_percent=0.0, humidity_timescale_minutes=30.0,
                         pressure_hpa=600.0, temperature_k=273.0, seed=42):
    """Return per-exposure conditions, with stationary correlated RH perturbations.

    ``humidity_sigma_percent`` is the stationary *absolute percentage-point*
    standard deviation, not a PWV or fractional humidity perturbation.
    The exponential correlation is exp(-delta_time / timescale); values are
    clipped to [0, 100] and clipping is recorded. Pressure and temperature
    remain fixed until an explicitly parameterized variability model exists.
    """
    from astropy.table import Table
    from astropy.time import Time
    n = len(schedule.midpoint_jd_utc)
    mids = np.asarray(schedule.midpoint_jd_utc, float)
    alt = np.asarray(schedule.altitude_deg, float)
    am = np.asarray(schedule.airmass, float)
    berv = None if schedule.berv_kms is None else np.asarray(schedule.berv_kms, float)
    if berv is None or any(x.shape != (n,) for x in (mids, alt, am, berv)) or n < 1:
        raise ValueError('A complete schedule with matching BERV values is required')
    if not all(np.isfinite(x).all() for x in (mids, alt, am, berv)) or np.any(np.diff(mids) <= 0):
        raise ValueError('Schedule must be finite and strictly increasing')
    if np.any(alt <= 0) or np.any(alt > 90) or np.any(am < 1):
        raise ValueError('Invalid target altitude or airmass')
    values = (humidity_percent, humidity_sigma_percent, humidity_timescale_minutes,
              pressure_hpa, temperature_k, schedule.duration_seconds)
    if not np.all(np.isfinite(values)) or not (0 <= humidity_percent <= 100) or humidity_sigma_percent < 0 or humidity_timescale_minutes <= 0 or pressure_hpa <= 0 or temperature_k <= 0 or schedule.duration_seconds <= 0:
        raise ValueError('Invalid atmospheric-state parameters')
    seed = int(seed)
    rng = np.random.default_rng(seed)
    z = np.empty(n)
    z[0] = rng.normal() if humidity_sigma_percent > 0 else 0.0
    for k in range(1, n):
        if humidity_sigma_percent > 0:
            dt_min = (mids[k] - mids[k-1]) * 1440.0
            rho = np.exp(-dt_min / humidity_timescale_minutes)
            z[k] = rho*z[k-1] + np.sqrt(max(0., 1. - rho*rho))*rng.normal()
        else:
            z[k] = 0.0
    humidity_raw = humidity_percent + humidity_sigma_percent*z
    humidity = np.clip(humidity_raw, 0., 100.)
    times = Time(mids, format='jd', scale='utc')
    table = Table()
    table['exposure_id'] = np.arange(n, dtype=int)
    table['midpoint_jd_utc'] = mids
    table['mjd_utc'] = times.mjd
    table['utc'] = times.isot.astype(str)
    table['exposure_seconds'] = np.full(n, schedule.duration_seconds)
    table['airmass'] = am
    table['berv_kms'] = berv
    table['altitude_deg'] = alt
    table['zenith_angle_deg'] = 90. - alt
    table['humidity_percent'] = humidity
    table['pressure_hpa'] = np.full(n, pressure_hpa)
    table['temperature_k'] = np.full(n, temperature_k)
    table['humidity_clipped'] = humidity_raw != humidity
    for name, unit in [('midpoint_jd_utc','d'), ('mjd_utc','d'), ('exposure_seconds','s'),
                       ('berv_kms','km / s'), ('altitude_deg','deg'),
                       ('zenith_angle_deg','deg'), ('humidity_percent','percent'),
                       ('pressure_hpa','hPa'), ('temperature_k','K')]:
        table[name].unit = unit
    table.meta.update({'schema': 'exoatmosim.observing_conditions.v1',
                       'time_scale': 'utc', 'wavelength_medium': 'air',
                       'humidity_model': 'stationary_exponential_gaussian_clipped',
                       'humidity_mean_percent': float(humidity_percent),
                       'humidity_sigma_percentage_points': float(humidity_sigma_percent),
                       'humidity_timescale_minutes': float(humidity_timescale_minutes),
                       'atmosphere_seed': seed,
                       'humidity_clipped_count': int(np.count_nonzero(table['humidity_clipped']))})
    return table


def write_observing_conditions(path, table):
    """Write an Astropy ECSV whose physical units and provenance survive reads."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write(path, format='ascii.ecsv', overwrite=True)
    return path


def read_observing_conditions(path, schedule=None):
    """Read and optionally verify an ECSV against the original schedule."""
    from astropy.table import Table
    table = Table.read(path, format='ascii.ecsv')
    required = ('exposure_id', 'midpoint_jd_utc', 'airmass', 'berv_kms',
                'zenith_angle_deg', 'humidity_percent', 'pressure_hpa', 'temperature_k')
    if any(c not in table.colnames for c in required):
        raise ValueError('Observing conditions table is missing required columns')
    ids = np.asarray(table['exposure_id'], int)
    if not np.array_equal(ids, np.arange(len(table))):
        raise ValueError('Exposure IDs must be consecutive and zero-based')
    for name in ('midpoint_jd_utc','airmass','berv_kms','zenith_angle_deg',
                 'humidity_percent','pressure_hpa','temperature_k'):
        if not np.all(np.isfinite(np.asarray(table[name], float))):
            raise ValueError(f'Non-finite observing conditions: {name}')
    rh = np.asarray(table['humidity_percent'], float)
    if np.any((rh < 0) | (rh > 100)):
        raise ValueError('Relative humidity must be between 0 and 100%')
    if schedule is not None:
        if len(table) != len(schedule.midpoint_jd_utc) or not np.allclose(
                np.asarray(table['midpoint_jd_utc'], float), schedule.midpoint_jd_utc,
                atol=1e-10, rtol=0):
            raise ValueError('Conditions timestamps do not match schedule')
    return table


def plot_observing_conditions(table, path, *, snr_fits=None, tellurics_enabled=True):
    """Plot PWV (or humidity), airmass, and optional per-exposure VIOLA S/N.

    The VIOLA native backend writes a wavelength-dependent S/N per resolution
    element in snr_res_0.fits. Plot the median over finite positive wavelength
    pixels for each exposure, avoiding a misleading single-pixel S/N.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from exoatmosim.plotting import publication_style
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    t = (np.asarray(table['midpoint_jd_utc'], float) -
         float(table['midpoint_jd_utc'][0])) * 24
    snr = None
    if snr_fits is not None:
        from astropy.io import fits
        snr_grid = np.asarray(fits.getdata(snr_fits), dtype=float)
        if snr_grid.ndim != 2 or snr_grid.shape[0] != len(t):
            raise ValueError('VIOLA S/N FITS must have shape (exposures, pixels)')
        valid = np.isfinite(snr_grid) & (snr_grid > 0)
        if not np.all(valid.any(axis=1)):
            raise ValueError('Each exposure must contain valid positive S/N pixels')
        snr = np.array([np.median(row[mask]) for row, mask in zip(snr_grid, valid)])
    with publication_style():
        fig, axes = plt.subplots(3 if snr is not None else 2, 1,
                                 figsize=(10, 8 if snr is not None else 6),
                                 sharex=True, constrained_layout=True)
        if not tellurics_enabled:
            axes[0].text(0.5, 0.5, 'No tellurics considered',
                         transform=axes[0].transAxes, ha='center', va='center')
            axes[0].set_ylabel('PWV [mm]')
            axes[0].set_yticks([])
        elif 'pwv_mm' in table.colnames:
            axes[0].plot(t, table['pwv_mm'])
            axes[0].set_ylabel('PWV [mm]')
        else:
            axes[0].plot(t, table['humidity_percent'])
            axes[0].set_ylabel('Relative humidity [%]')
        axes[1].plot(t, table['airmass'])
        axes[1].set_ylabel('Airmass')
        if snr is not None:
            axes[2].plot(t, snr, marker='o', markersize=3)
            axes[2].set_ylabel('Median S/N / res. element')
            axes[2].ticklabel_format(axis='y', style='plain', useOffset=False)
        axes[-1].set_xlabel('Hours since first exposure')
        fig.align_ylabels(axes)
        fig.savefig(path, dpi=180)
        plt.close(fig)
    return path


def sample_pwv(schedule, *, mean_mm=2., sigma_mm=0.1, timescale_minutes=30., seed=42):
    """Stationary, seeded correlated PWV values in mm, one per exposure."""
    mids=np.asarray(schedule.midpoint_jd_utc,float)
    if (mids.ndim!=1 or len(mids)==0 or np.any(np.diff(mids)<=0) or
        not np.isfinite([mean_mm,sigma_mm,timescale_minutes]).all() or
        mean_mm<=0 or sigma_mm<0 or timescale_minutes<=0):
        raise ValueError('Invalid PWV sampling parameters or timestamps')
    rng=np.random.default_rng(seed)
    z=np.zeros(len(mids))
    if sigma_mm: z[0]=rng.normal()
    for i in range(1,len(mids)):
        rho=np.exp(-(mids[i]-mids[i-1])*1440/timescale_minutes)
        z[i]=rho*z[i-1]+np.sqrt(max(0.,1.-rho*rho))*rng.normal() if sigma_mm else 0.
    # strict positive PWV; grid-range checks happen in grid_transmission
    return np.maximum(0.001,mean_mm+sigma_mm*z)
