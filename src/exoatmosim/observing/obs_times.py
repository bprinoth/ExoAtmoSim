"""Parse tayph-compatible obs_times files without guessing Astropy column names."""
from pathlib import Path
import numpy as np


def read_obs_times_mjd(path: str | Path) -> np.ndarray:
    """Return the MJD in the first column of a headerless or headered table.

    A header is allowed only before any numeric rows. Reject empty, malformed,
    non-finite and non-increasing sequences rather than reverting to a fake orbit.
    """
    path = Path(path)
    timestamps = []
    with path.open(encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            first = line.split()[0]
            try:
                mjd = float(first)
            except ValueError as exc:
                if not timestamps and first.lower() in {'mjd', 'time_mjd', 'obstime', 'time'}:
                    continue
                raise ValueError(f'{path}:{line_no}: invalid MJD {first!r}') from exc
            if not np.isfinite(mjd):
                raise ValueError(f'{path}:{line_no}: non-finite MJD')
            timestamps.append(mjd)
    if not timestamps:
        raise ValueError(f'{path}: no observation timestamps found')
    values = np.array(timestamps, dtype=float)
    if np.any(np.diff(values) <= 0):
        raise ValueError(f'{path}: observation timestamps must be strictly increasing')
    return values


def write_tayph_obs_times(path: str | Path, schedule) -> Path:
    """Write tayph's five-column observation table from the physical schedule.

    Columns: MJD(UTC), ISO(UTC), exposure seconds, geometric airmass,
    barycentric radial-velocity correction in km/s. All quantities correspond
    to the exposure midpoint. Does not replace airmass with exposure time.
    """
    from astropy.time import Time

    path = Path(path)
    mids = np.asarray(schedule.midpoint_jd_utc, dtype=float)
    airmass = np.asarray(schedule.airmass, dtype=float)
    if schedule.berv_kms is None:
        raise ValueError('Cannot write complete obs_times: scheduler has no BERV; use compute_berv=True')
    berv = np.asarray(schedule.berv_kms, dtype=float)
    if mids.ndim != 1 or any(a.shape != mids.shape for a in (airmass, berv)):
        raise ValueError('MJD, airmass and BERV must be matching one-dimensional exposure arrays')
    if not len(mids) or np.any(~np.isfinite(mids)) or np.any(np.diff(mids) <= 0):
        raise ValueError('Observation midpoints must be finite, nonempty and increasing')
    if np.any(~np.isfinite(airmass)) or np.any(airmass < 1):
        raise ValueError('Observation contains undefined/invalid airmass (check target is above horizon)')
    if np.any(~np.isfinite(berv)):
        raise ValueError('Observation contains non-finite barycentric velocities')
    exposure = float(schedule.duration_seconds)
    if not np.isfinite(exposure) or exposure <= 0:
        raise ValueError('Exposure duration must be positive and finite')
    times = Time(mids, format='jd', scale='utc')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as handle:
        for mjd, iso, am, bc in zip(times.mjd, times.isot, airmass, berv):
            handle.write(f'{mjd:.11f}\t{iso}\t{exposure:.4f}\t{am:.8f}\t{bc:.8f}\n')
    return path
