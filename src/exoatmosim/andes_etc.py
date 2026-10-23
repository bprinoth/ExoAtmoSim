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


def spectrum_snr(grid: AndesETCGrid, wavelength_nm, airmass, resolving_power):
    """Return per-sample SNR for ANDES concatenated optical+IR segments.

    The ETC provides SNR per resolution element. Each simulated wavelength
    sample is assigned to one ETC arm; never average across arms. Convert
    to independent-sample SNR using local pixels/resolution-element.
    Assumes the model sampling is sufficiently finer than a resolution
    element; the observation noise model remains an approximation.
    """
    wavelength_nm = np.asarray(wavelength_nm, float)
    airmass = np.asarray(airmass, float)
    if wavelength_nm.ndim != 1 or len(wavelength_nm) < 2 or np.any(np.diff(wavelength_nm) <= 0):
        raise ValueError('ANDES wavelength samples must be strictly increasing')
    if not np.isfinite(resolving_power) or resolving_power <= 0:
        raise ValueError('Invalid ANDES resolving power')
    out = np.empty((len(airmass), len(wavelength_nm)))
    # At arm overlaps prefer UBV below 621.7, RIZ up to 950,
    # and YJH above 950; boundaries follow actual loaded data.
    bounds = {arm: (grid.data[arm]['wavelength'][0],
                    grid.data[arm]['wavelength'][-1]) for arm in ARMS}
    # Use nearest arm to the wavelength, favoring YJH at its lower edge.
    choices = np.full(len(wavelength_nm), -1, dtype=int)
    for k, arm in enumerate(ARMS):
        lo, hi = bounds[arm]
        mask = (wavelength_nm >= lo) & (wavelength_nm <= hi)
        choices[mask] = k
    if np.any(choices < 0):
        missing = wavelength_nm[choices < 0]
        raise ValueError(f'ETC does not cover ANDES wavelengths: {missing.min():.4f}-{missing.max():.4f} nm')
    for k, arm in enumerate(ARMS):
        selected = choices == k
        if np.any(selected):
            out[:, selected] = grid.evaluate(arm, wavelength_nm[selected], airmass).snr_resel

    # Local sample width in log-lambda, appropriate even for irregular grids.
    dln = np.gradient(np.log(wavelength_nm))
    pixels_per_resel = 1. / (float(resolving_power) * dln)
    if not np.isfinite(pixels_per_resel).all() or np.any(pixels_per_resel < 1):
        raise ValueError('ANDES order sampling is not at least one sample per resolution element')
    snr_pixel = out / np.sqrt(pixels_per_resel[None, :])
    if not np.isfinite(snr_pixel).all() or np.any(snr_pixel <= 0):
        raise ValueError('ETC produced zero or non-finite S/N; inspect ETC coverage/noise')
    return snr_pixel
