"""Instrument resource discovery, without copying order files per planet.

The actual tayph order_def file schema is backend-specific; validation of an
instrument order file for the lightweight simulator is separate from tayph.
"""
from pathlib import Path
from importlib.resources import files
from .config import Config


def resolve_instrument_file(config: Config, config_path: str | Path, *, key: str = 'instrument_order_file') -> Path:
    """Find an explicit file: absolute, planet-relative, or shared-library-relative.

    Shared root is selected with ``instrument_library_path`` in the same config.
    Explicit resolution rather than implicit CWD searches prevents accidental
    reuse of another planet's instrument setup.
    """
    if key not in config:
        if 'instrument' not in config:
            raise ValueError(f'Missing {key} or instrument: select a bundled instrument or explicit order table')
        return instrument_resource(config['instrument'])
    raw = Path(config[key]).expanduser()
    config_dir = Path(config_path).resolve().parent
    if raw.is_absolute():
        path = raw
    elif 'instrument_library_path' in config:
        root = Path(config['instrument_library_path']).expanduser()
        root = root if root.is_absolute() else config_dir / root
        path = root / raw
    else:
        path = config_dir / raw
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f'Instrument definition not found: {path}')
    return path


# The bundled inputs are measured wavelength grids, not (min,max) tables.
# Missing instrument order files are intentionally NOT replaced by guesses.
_BUNDLED = {
    'CRIRES_K': ('CRIRES_K_orders_nm.npy', None),
    'ANDES_K': ('ANDES_K_orders_nm.npy', None),
    'VIOLA_K': ('VIOLA_K_orders_nm.npy', None),
    'VIOLA_H': ('VIOLA_H_orders_nm.npy', None),
    'VIOLA_Y': ('VIOLA_Y_orders_nm.npy', None),
    'IGRINS_HK': ('IGRINS_HK_orders_nm.npy', None),
    'IGRINS2_HK': ('IGRINS2_HK_orders_nm.npy', None),
    'SPIROU': ('SPIROU_orders_nm.npy', 'SPIROU_blaze.npy'),
    'HARPSN': (None, 'HARPSN_blaze.npy'),
    'HARPS': (None, None),
    'ESPRESSO': ('ESPRESSO_orders_nm.npy', 'ESPRESSO_blaze.npy'),
}
_ALIASES = {'HARPS-N': 'HARPSN', 'IGRINS-2': 'IGRINS2_HK', 'CRIRES+': 'CRIRES_K'}


def canonical_instrument(instrument: str) -> str:
    """Normalize supported instrument aliases for filenames and tayph mode."""
    name = instrument.upper().replace(' ', '_')
    name = _ALIASES.get(name, name)
    if name not in _BUNDLED:
        raise ValueError(f'Unknown instrument {instrument!r}; available: {", ".join(_BUNDLED)}')
    return name


def available_instruments() -> dict[str, dict[str, bool]]:
    """Report supplied (not extrapolated) instrument resources."""
    return {name: {'orders': order is not None, 'blaze': blaze is not None}
            for name, (order, blaze) in _BUNDLED.items()}


def instrument_resource(instrument: str, kind: str = 'orders') -> Path:
    """Return path to a packaged resource; fail transparently for missing data.

    For zip-imported distributions callers should use importlib.resources.as_file;
    ExoAtmoSim's normal wheel installation installs these resources as files.
    """
    name = canonical_instrument(instrument)
    if kind not in ('orders', 'blaze'):
        raise ValueError('kind must be orders or blaze')
    filename = _BUNDLED[name][0 if kind == 'orders' else 1]
    if filename is None:
        raise FileNotFoundError(f'{name}: no valid bundled {kind} file. '
                                'Provide a real calibrated resource; the uploaded order file was empty or missing.')
    resource = files('exoatmosim.resources.instruments').joinpath(filename)
    if not resource.is_file():
        raise FileNotFoundError(f'Missing packaged resource: {filename}. '
                                'For VIOLA, export actual ETC wavelength samples using '
                                'examples/export_viola_orders.py; no synthetic grid is provided.')
    return Path(str(resource))


def instrument_wavelength_bounds_um(config: Config, config_path: str | Path) -> tuple[float, float]:
    """Full measured order coverage plus an explicit margin, in microns.

    Values are measured from the bundled per-pixel order definition; no
    hard-coded wavelength limits are used for any instrument.
    ``instrument_margin_nm`` defaults to 5 nm per side. Explicit ``wl_min``
    and ``wl_max`` can override either calculated boundary independently.
    """
    import numpy as np
    order_path = resolve_instrument_file(config, config_path)
    data = np.asarray(np.load(order_path, allow_pickle=False) if order_path.suffix == '.npy'
                      else np.loadtxt(order_path), dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f'Invalid instrument order table shape: {data.shape}')
    valid = data[np.isfinite(data) & (data > 0)]
    if valid.size < 2:
        raise ValueError('Instrument order table contains no valid wavelength coverage')
    margin = config.get_float('instrument_margin_nm', 5.0)
    if not np.isfinite(margin) or margin < 0:
        raise ValueError('instrument_margin_nm must be finite and nonnegative')
    low = config.get_float('wl_min', (float(valid.min()) - margin) / 1000.0)
    high = config.get_float('wl_max', (float(valid.max()) + margin) / 1000.0)
    if not np.isfinite(low) or not np.isfinite(high) or low <= 0 or low >= high:
        raise ValueError('Invalid wavelength bounds from instrument orders and overrides')
    if low * 1000 > valid.min() or high * 1000 < valid.max():
        raise ValueError('Explicit wavelength bounds do not cover all instrument orders')
    return low, high
