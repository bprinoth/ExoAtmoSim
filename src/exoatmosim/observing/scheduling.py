"""Independent exposure scheduler for ExoAtmoSim (no external mock imports).

Times are exposure MIDPOINTS. ``Tobs`` is an observing-window midpoint as
Julian Date UTC. Without it, find the nearest upper meridian passage (maximum
geometric altitude) to ``Tc``; ``Tc`` is used only as an approximate search
anchor, not as the observation epoch in the output time scale.

This is intentionally distinct from the legacy mock's phase-based timing.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math
import numpy as np

from exoatmosim.config import Config
from exoatmosim.plotting import publication_style


@dataclass(frozen=True)
class ObservationSchedule:
    center_jd_utc: float
    start_jd_utc: np.ndarray
    midpoint_jd_utc: np.ndarray
    end_jd_utc: np.ndarray
    duration_seconds: float
    readout_seconds: float
    window_hours: float
    altitude_deg: np.ndarray
    airmass: np.ndarray
    mode: str
    berv_kms: np.ndarray | None = None

    @property
    def n_exposures(self) -> int:
        return len(self.midpoint_jd_utc)


def _as_config(config: Config | str | Path) -> Config:
    return config if isinstance(config, Config) else Config.from_file(config)


def _exposure_seconds(cfg: Config) -> float:
    """Read the single, required per-exposure integration time (seconds)."""
    if 'exposure_time' not in cfg:
        raise ValueError('Missing exposure_time in config (seconds)')
    value = cfg.get_float('exposure_time')
    if not np.isfinite(value) or value <= 0:
        raise ValueError('exposure_time must be finite and positive (seconds)')
    return value


def exposure_midpoints(center_jd: float, hours: float, texp: float, readout: float):
    """Return complete exposures centered inside a window, no extra overhead.

    The first exposure starts at the window start, with no leading readout.
    Trailing unused seconds are permitted. All timestamps are in JD.
    """
    vals = [center_jd, hours, texp, readout]
    if not all(np.isfinite(vals)) or hours <= 0 or texp <= 0 or readout < 0:
        raise ValueError('center, hours, texp and readout must be finite, with positive window/exposure')
    available_s = hours * 3600
    n = int(np.floor((available_s + readout + 1e-9) / (texp + readout)))
    if n < 1:
        raise ValueError('Observing window is shorter than one exposure')
    start = center_jd - hours / 48.0 + np.arange(n, dtype=float) * (texp + readout) / 86400.0
    middle = start + texp / (2.0 * 86400.0)
    end = start + texp / 86400.0
    if end[-1] > center_jd + hours / 48 + 1e-9:
        raise AssertionError('An exposure exceeds the observing window')
    return start, middle, end


def _sky_location(cfg: Config):
    try:
        from astropy.coordinates import SkyCoord, EarthLocation
        import astropy.units as u
    except ImportError as exc:
        raise ImportError('Astropy is required for target altitude and BERV. Install exoatmosim[orbits].') from exc
    # legacy config convention: RA hourangle, DEC degrees, lat/long degrees.
    target = SkyCoord(ra=cfg['RA'], dec=cfg['DEC'], unit=(u.hourangle, u.deg), frame='icrs')
    if str(cfg.get('observatory', '')).strip().lower() in ('armazones', 'cerro armazones'):
        # ESO ELT site: 24°36' S, 70°11'30" W, ~3046 m.
        # The NASA template's legacy lat/long/elev describe Paranal;
        # do not silently use those for ANDES.
        location = EarthLocation.from_geodetic(
            lon=(-70 - 11/60 - 30/3600) * u.deg,
            lat=-24.6 * u.deg,
            height=3046 * u.m)
    else:
        location = EarthLocation.from_geodetic(
            lon=cfg.get_float('long') * u.deg,
            lat=cfg.get_float('lat') * u.deg,
            height=cfg.get_float('elev') * u.m,
        )
    return target, location


def _altitude(target, location, jd_utc):
    from astropy.coordinates import AltAz
    from astropy.time import Time
    t = Time(jd_utc, format='jd', scale='utc')
    return np.atleast_1d(target.transform_to(AltAz(obstime=t, location=location, pressure=0)).alt.deg)


def nearest_maximum_altitude_jd(anchor_jd_utc: float, target, location) -> float:
    """Find geometric upper culmination nearest an approximate UTC JD anchor.

    We maximize altitude on a sampled sidereal-day interval; quadratic
    interpolation gives sub-minute placement. Not a weather/visibility optimizer.
    """
    sidereal_day = 0.9972695663
    trial = anchor_jd_utc + np.linspace(-0.5, 0.5, 1441) * sidereal_day
    alt = _altitude(target, location, trial)
    i = int(np.argmax(alt))
    if 0 < i < len(trial) - 1:
        x = (alt[i - 1] - alt[i + 1]) / (2.0 * (alt[i - 1] - 2 * alt[i] + alt[i + 1])) if alt[i-1] - 2*alt[i] + alt[i+1] != 0 else 0.0
        return float(trial[i] + np.clip(x, -1, 1) * (trial[1] - trial[0]))
    return float(trial[i])


def observation_window_hours(cfg: Config) -> float:
    """Use transit duration plus one hour of baseline on each side.

    ``T14`` comes from NASA Exoplanet Archive, in HOURS. An explicit
    ``observation_window_hours`` takes precedence. Without that setting,
    use ``T14 + 2`` hours to provide one hour of baseline on either side.
    Historical ``t`` is not an observing-window setting.
    """
    if 'observation_window_hours' in cfg:
        hours = cfg.get_float('observation_window_hours')
    elif 'T14' in cfg and math.isfinite(cfg.get_float('T14')) and cfg.get_float('T14') > 0:
        hours = cfg.get_float('T14') + 2.0
    else:
        raise ValueError('Set observation_window_hours (total hours), or provide a positive T14 (transit duration in hours)')
    if not np.isfinite(hours) or hours <= 0:
        raise ValueError('Observation-window hours must be positive and finite')
    return hours


def make_observation_schedule(config: Config | str | Path, *, compute_berv: bool = True) -> ObservationSchedule:
    """Schedule exposures by Tobs or nearest upper culmination, then evaluate geometry.

    ``Tobs``: Julian Date UTC of observing-window *center*.
    ``Tc``: approximate JD of transit center used only to select the nearby
    maximum-altitude passage when Tobs is missing. Tc may be BJD_TDB; this
    introduces a <~minute search-anchor offset, NOT a timestamp conversion.
    """
    cfg = _as_config(config)
    target, location = _sky_location(cfg)
    duration = _exposure_seconds(cfg)
    readout = cfg.get_float('readout', 0.0)
    hours = observation_window_hours(cfg)
    if 'Tobs' in cfg:
        center = cfg.get_float('Tobs')
        mode = 'Tobs'
    else:
        center = nearest_maximum_altitude_jd(cfg.get_float('Tc'), target, location)
        mode = 'max_altitude'
    start, mids, end = exposure_midpoints(center, hours, duration, readout)
    altitude = _altitude(target, location, mids)
    # True geometric airmass: do not clip a below-horizon target to 1.0.
    z = 90.0 - altitude
    airmass = np.full_like(altitude, np.nan)
    visible = altitude > 0
    airmass[visible] = 1.0 / (np.cos(np.deg2rad(z[visible])) + 0.50572 * (96.07995 - z[visible])**(-1.6364))
    berv = None
    if compute_berv:
        import astropy.units as u
        from astropy.time import Time
        time = Time(mids, format='jd', scale='utc', location=location)
        berv = np.asarray(target.radial_velocity_correction(kind='barycentric', obstime=time).to_value(u.km / u.s))
    return ObservationSchedule(center, start, mids, end, duration, readout, hours, altitude, airmass, mode, berv)


def save_schedule(schedule: ObservationSchedule, directory: str | Path, *, plots: bool = True) -> list[Path]:
    """Write a reproducible tabular schedule and observation-only diagnostics."""
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    output = path / 'observation_schedule.csv'
    with output.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['start_jd_utc', 'mid_jd_utc', 'end_jd_utc', 'altitude_deg', 'airmass', 'berv_kms'])
        for i in range(schedule.n_exposures):
            writer.writerow([schedule.start_jd_utc[i], schedule.midpoint_jd_utc[i], schedule.end_jd_utc[i], schedule.altitude_deg[i], schedule.airmass[i], schedule.berv_kms[i] if schedule.berv_kms is not None else ''])
    outputs = [output]
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        x = (schedule.midpoint_jd_utc - schedule.center_jd_utc) * 24
        png = path / 'observation_geometry.png'
        with publication_style():
            fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
            axes[0].plot(x, schedule.altitude_deg, 'o-', ms=3)
            axes[0].set_ylabel('Altitude (deg)')
            axes[1].plot(x, schedule.airmass, 'o-', ms=3)
            axes[1].set_ylabel('Airmass')
            axes[1].set_xlabel('Hours relative to observing-window center')
            if schedule.mode == 'max_altitude':
                axes[0].axvline(0, color='0.45', ls='--', lw=1.5, label='Idealized mid-transit / culmination')
                axes[0].legend()
            fig.tight_layout()
            fig.savefig(png, dpi=180)
            plt.close(fig)
        outputs.append(png)
    return outputs
