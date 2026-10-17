"""Lightweight, explicit synthetic spectroscopic observations.

This is *not* a replacement for previous external mock implementation: no blaze, detector model,
stellar line spectrum or tellurics are included. All inputs are testable.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from .config import Config
from .instruments import resolve_instrument_file


@dataclass(frozen=True)
class SpectralOrders:
    """Each row lists wavelength endpoints [nm] for one detector order.

    For .npy inputs the accepted format is strictly shape (n_orders, 2).
    Other legacy order definition encodings must be converted explicitly, rather
    than silently interpreting pixel/indices as wavelengths.
    """
    ranges_nm: np.ndarray

    def __post_init__(self):
        a = np.asarray(self.ranges_nm, dtype=float)
        if a.ndim != 2 or a.shape[1] != 2 or len(a) == 0:
            raise ValueError('Order limits must have shape (n_orders, 2), in nm')
        if not np.isfinite(a).all() or np.any(a[:, 0] <= 0) or np.any(a[:, 1] <= a[:, 0]):
            raise ValueError('Invalid order wavelength limits in nm')
        object.__setattr__(self, 'ranges_nm', a.copy())

    @classmethod
    def load(cls, path: str | Path) -> 'SpectralOrders':
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f'Instrument order definition missing: {path}')
        if path.suffix == '.npy':
            array = np.load(path, allow_pickle=False)
        elif path.suffix in {'.txt', '.dat', '.csv'}:
            array = np.loadtxt(path, comments='#', delimiter=',' if path.suffix == '.csv' else None, ndmin=2)
        else:
            raise ValueError('Supported order formats: .npy, .dat, .txt, .csv')
        # Real previous mock framework order definitions are per-pixel wavelength grids.
        # Reduce each row to its covered window for the approximate simulator.
        # The original unmodified grids remain bundled for external integration.
        if array.ndim == 2 and array.shape[1] > 2:
            if not np.isfinite(array).all() or np.any(array <= 0):
                raise ValueError('Wavelength grid contains invalid values')
            array = np.column_stack((array.min(axis=1), array.max(axis=1)))
        return cls(array)


@dataclass(frozen=True)
class ObservationSettings:
    order_file: Path
    snr: float
    resolving_power: float
    pixels_per_resolution_element: float = 3.0
    random_seed: int = 42

    @classmethod
    def from_config(cls, config_path: str | Path) -> 'ObservationSettings':
        path = Path(config_path).resolve()
        config = Config.from_file(path)
        order_file = resolve_instrument_file(config, path)
        snr = config.get_float('mock_snr')
        resolving_power = config.get_float('mock_resolution', 140000.)
        sampling = config.get_float('mock_pixels_per_resolution_element', 3.)
        seed = config.get_int('mock_seed', 42)
        if not np.isfinite(snr) or snr <= 0 or not np.isfinite(resolving_power) or resolving_power <= 0:
            raise ValueError('mock_snr and mock_resolution must be positive finite numbers')
        if not np.isfinite(sampling) or sampling < 2:
            raise ValueError('mock_pixels_per_resolution_element must be at least 2')
        return cls(order_file, snr, resolving_power, sampling, seed)


@dataclass(frozen=True)
class SyntheticObservation:
    wavelength_nm: tuple[np.ndarray, ...]
    flux: tuple[np.ndarray, ...]
    uncertainty: tuple[np.ndarray, ...]
    rv_planet_kms: np.ndarray
    phase: np.ndarray

    def save(self, path: str | Path) -> Path:
        """Explicit opt-in save; ragged orders are separate numbered NPZ fields."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        values = {'rv_planet_kms': self.rv_planet_kms, 'phase': self.phase,
                  'n_orders': np.array(len(self.wavelength_nm))}
        for i, (wave, flux, err) in enumerate(zip(self.wavelength_nm, self.flux, self.uncertainty)):
            values[f'wavelength_nm_{i}'] = wave
            values[f'flux_{i}'] = flux
            values[f'uncertainty_{i}'] = err
        np.savez_compressed(path, **values)
        return path

    def plot(self, path: str | Path, order_index: int = 0) -> Path:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        data = np.asarray(self.flux[order_index], dtype=float)
        if data.ndim != 2 or data.shape[0] != len(self.phase) or data.shape[1] != len(self.wavelength_nm[order_index]):
            raise ValueError('Mock flux shape disagrees with exposure/wavelength axes')
        if not np.isfinite(data).all():
            bad = int(np.count_nonzero(~np.isfinite(data)))
            raise ValueError(f'Cannot plot mock order {order_index}: {bad}/{data.size} non-finite flux pixels')
        fig, ax = plt.subplots(figsize=(10, 5))
        img = ax.pcolormesh(self.wavelength_nm[order_index], np.arange(len(self.phase)),
                             data, shading='auto', cmap='viridis')
        fig.colorbar(img, ax=ax, label='Normalised stellar flux')
        ax.set(xlabel='Wavelength [nm]', ylabel='Exposure index', title='Synthetic order (flat stellar continuum)')
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path


def simulate_transmission(transmission, orders: SpectralOrders, settings: ObservationSettings,
                          phase: np.ndarray, rv_planet_kms: np.ndarray,
                          transit_weight: np.ndarray | None = None) -> SyntheticObservation:
    """Shift differential transit signal onto log-lambda pixels; add Gaussian noise.

    Input ``transmission`` is the v0.7.1 result named transmission = 1 -
    POSEIDON output. We subtract its median baseline; only the differential
    spectral signature is injected, *not* the broadband transit light curve.
    This approximate simulation is intended for integration tests, not SNR/yield
    forecasts. No exposure smearing, limb darkening, star, blaze or tellurics.
    """
    phase = np.asarray(phase, float)
    rv = np.asarray(rv_planet_kms, float)
    if phase.ndim != 1 or len(phase) == 0 or rv.shape != phase.shape:
        raise ValueError('phase and rv_planet_kms must be matching 1-D arrays')
    weight = np.ones_like(phase) if transit_weight is None else np.asarray(transit_weight, float)
    if weight.shape != phase.shape or np.any((weight < 0) | (weight > 1)):
        raise ValueError('Transit weights must be within [0,1]')
    wave = np.asarray(transmission.wavelength_um, float) * 1000
    signal = np.asarray(transmission.transmission, float)
    if wave.ndim != 1 or signal.shape != wave.shape or np.any(np.diff(wave) <= 0):
        raise ValueError('Input transmission wavelengths must increase')
    if not np.isfinite(signal).all():
        raise ValueError('POSEIDON transmission contains non-finite values before mock injection')
    differential = signal - np.median(signal)
    rng = np.random.default_rng(settings.random_seed)
    wout, fout, eout = [], [], []
    for lo, hi in orders.ranges_nm:
        if lo < wave[0] or hi > wave[-1]:
            raise ValueError(f'Order {lo:g}-{hi:g} nm outside POSEIDON wavelength interval')
        dlog = 1 / (settings.resolving_power * settings.pixels_per_resolution_element)
        n = int(np.floor(np.log(hi / lo) / dlog)) + 1
        if n < 2 or n > 2_000_000:
            raise ValueError('Order yields invalid or excessive number of pixels')
        w = lo * np.exp(np.arange(n) * dlog)
        shifted = w[None, :] / (1 + rv[:, None] / 299792.458)
        model = 1 - weight[:, None] * np.interp(shifted, wave, differential, left=0, right=0)
        sigma = np.full(model.shape, 1 / settings.snr)
        noisy = model + rng.normal(0, 1 / settings.snr, size=model.shape)
        if not np.isfinite(noisy).all():
            raise ValueError(f'Non-finite simulated flux for order {lo:g}-{hi:g} nm; inspect injection and SNR')
        wout.append(w)
        fout.append(noisy)
        eout.append(sigma)
    return SyntheticObservation(tuple(wout), tuple(fout), tuple(eout), rv.copy(), phase.copy())
