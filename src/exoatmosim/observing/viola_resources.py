"""Export real VIOLA wavelength samples to bundled ExoAtmoSim order resources.

The VIOLA ETC wavelength output, not tayph.mock, defines these sample locations.
No synthetic order boundaries are invented; one contiguous segment is used until
VIOLA supplies documented detector/order segmentation.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

from .viola import read_viola_csv

BANDS = ('K', 'H', 'Y')


def export_viola_order_grid(csv_path, band, output_dir, *, wavelength_medium='air', overwrite=False):
    """Save a (1, Npixel) measured air-nm wavelength resource from VIOLA CSV.

    These are order *sampling grids*; a single row is NOT a claim that the
    instrument has just one physical echelle order.
    """
    band = str(band).upper()
    if band not in BANDS:
        raise ValueError(f'VIOLA band must be one of {BANDS}')
    wave, _, _ = read_viola_csv(csv_path, wavelength_medium=wavelength_medium)
    if wave.size < 2:
        raise ValueError('VIOLA requires at least two wavelength samples')
    output = Path(output_dir) / f'VIOLA_{band}_orders_nm.npy'
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not overwrite:
        existing = np.load(output, allow_pickle=False)
        if existing.shape == (1, wave.size) and np.array_equal(existing[0], wave):
            return output
        raise FileExistsError(f'{output} already exists with a different grid; pass overwrite=True deliberately')
    np.save(output, wave[None, :])
    return output
