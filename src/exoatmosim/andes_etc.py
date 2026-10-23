"""Read manually exported ANDES ETC v5.6 curves and interpolate in airmass.

Input wavelengths are ETC vacuum nm; outputs retain the original wavelength
convention. S/N is per resolution element. No implicit magnitude, seeing or
exposure-time scaling: all nine CSVs must share the same ETC settings.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import re
import numpy as np

AIRMASSES = (1.0, 1.5, 3.0)
ARMS = ('UBV', 'RIZ', 'YJH')
FIELDS = {
    'wavelength': 'wavelength [nm]',
    'snr': 'snr [per res. elem.]',
    'source': 'object signal [e-/res. elem.]',
    'noise': 'total noise [e-/res. elem.]',
    'background': 'background [e-/s/res. elem.]',
    'atmosphere': 'Efficiencies [Atm]',
}


def _resolve_csv(folder: Path, arm: str) -> Path:
    exact = folder / f'{arm}.csv'
    if exact.is_file():
        return exact
    # Preserve the official timestamp-based export filenames.
    candidates = [p for p in folder.glob('*.csv')
                  if re.search(r'(^|[_-])' + arm + r'(?=\.csv$)', p.name, re.I)]
    if len(candidates) != 1:
        raise FileNotFoundError(
            f'Expected exactly one {arm} ETC export in {folder} '
            f'(or a named {arm}.csv); found {len(candidates)}.')
    return candidates[0]


def _read_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline='') as fh:
        header = next(csv.reader(fh))
    missing = set(FIELDS.values()) - set(header)
    if missing:
        raise ValueError(f'{path}: missing ETC columns {sorted(missing)}')
    columns = np.loadtxt(path, delimiter=',', skiprows=1,
                         usecols=[header.index(FIELDS[name]) for name in FIELDS],
                         unpack=True, ndmin=2)
    if columns.ndim != 2 or columns.shape[0] != len(FIELDS):
        raise ValueError(f'{path}: invalid ETC numeric array')
    data = dict(zip(FIELDS, columns))
    wl = data['wavelength']
    # ETC rows are concatenated order/segment samples. Exported wavelength
    # precision can yield equal adjacent values at order boundaries. Preserve
    # those rows (and their distinct noise values) rather than deduplicating.
    if wl.size < 3 or not np.isfinite(wl).all() or np.any(np.diff(wl) < 0):
        raise ValueError(f'{path}: wavelengths must be finite and nondecreasing')
    for name, values in data.items():
        if not np.isfinite(values).all() or (name != 'wavelength' and np.any(values < 0)):
            raise ValueError(f'{path}: invalid {name} values')
    if not np.any(data['snr'] > 0):
        raise ValueError(f'{path}: S/N is zero everywhere')
    return data


@dataclass(frozen=True)
class ETCArm:
    wavelength_nm: np.ndarray
    snr_resel: np.ndarray
    source_e_resel: np.ndarray
    noise_e_resel: np.ndarray
    background_e_s_resel: np.ndarray
    atmosphere_efficiency: np.ndarray


class AndesETCGrid:
    """Validated 3-airmass by 3-arm local ETC grid.

    Interpolates source electrons and noise **variance** linearly in airmass,
    at each wavelength. S/N is then source/sqrt(variance), rather than a
    direct linear interpolation of S/N. Wavelength values use linear
    interpolation on the ETC grid. No extrapolation is permitted.
    """

    def __init__(self, directory: str | Path):
        self.directory = Path(directory)
        self.files = {}
        self.data = {}
        for arm in ARMS:
            rows = []
            for am in AIRMASSES:
                folder = self.directory / f'airmass_{am:.1f}'
                path = _resolve_csv(folder, arm)
                row = _read_csv(path)
                rows.append(row)
                self.files[(am, arm)] = path
            ref = rows[0]['wavelength']
            if any(len(row['wavelength']) != len(ref) or
                   not np.allclose(row['wavelength'], ref, rtol=0, atol=1e-7)
                   for row in rows[1:]):
                raise ValueError(f'{arm}: ETC wavelengths differ between airmass files; '
                                 'regenerate with consistent settings')
            self.data[arm] = {
                'wavelength': ref,
                **{name: np.stack([row[name] for row in rows])
                   for name in FIELDS if name != 'wavelength'},
            }

    def evaluate(self, arm: str, wavelength_nm, airmass) -> ETCArm:
        arm = str(arm).upper()
        if arm not in self.data:
            raise ValueError(f'Unknown ETC spectrograph {arm!r}')
        wl = np.asarray(wavelength_nm, dtype=float)
        am = np.asarray(airmass, dtype=float)
        if wl.ndim != 1 or wl.size < 1 or am.ndim != 1 or am.size < 1:
            raise ValueError('Wavelength and exposure airmass must be nonempty 1D arrays')
        if not np.isfinite(wl).all() or not np.isfinite(am).all():
            raise ValueError('ETC interpolation coordinates must be finite')
        x = self.data[arm]['wavelength']
        if wl.min() < x[0] or wl.max() > x[-1]:
            raise ValueError(f'{arm}: requested wavelengths outside ETC range '
                             f'[{x[0]:g}, {x[-1]:g}] nm')
        if am.min() < AIRMASSES[0] or am.max() > AIRMASSES[-1]:
            raise ValueError('ETC airmass outside calibrated range [1.0, 3.0]')

        # np.interp requires strictly increasing support. At equal exported
        # wavelengths use the mean of all rows at that rounded wavelength.
        # This preserves every row in the stored ETC grid and produces a
        # deterministic value at the otherwise ambiguous boundary.
        unique_x, starts, counts = np.unique(x, return_index=True,
                                               return_counts=True)
        if len(unique_x) != len(x):
            x_interp = unique_x
        else:
            x_interp = x

        index = np.clip(np.searchsorted(AIRMASSES, am, side='right') - 1, 0, 1)
        t = (am - np.asarray(AIRMASSES)[index]) / (
            np.asarray(AIRMASSES)[index+1] - np.asarray(AIRMASSES)[index])

        def interpolate(name: str, *, squared=False):
            y = self.data[arm][name]
            if len(x_interp) != len(x):
                y = np.stack([np.add.reduceat(row, starts) / counts for row in y])
            low = np.stack([np.interp(wl, x_interp, row) for row in y])
            if squared:
                low = low**2
            interp = low[index] * (1-t[:, None]) + low[index+1] * t[:, None]
            return np.sqrt(np.maximum(interp, 0)) if squared else interp

        source = interpolate('source')
        noise = interpolate('noise', squared=True)
        background = interpolate('background')
        atmospheric = interpolate('atmosphere')
        snr = np.divide(source, noise, out=np.zeros_like(source), where=noise > 0)
        return ETCArm(wl.copy(), snr, source, noise, background, atmospheric)
