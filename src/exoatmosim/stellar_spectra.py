"""Stellar spectra: generic data object and tayph/PHOENIX provider.

The provider imports tayph only when called, allowing lightweight tests and
other stellar-spectrum sources without a tayph dependency at import time.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np

from .config import Config
from .objects import Star
from .backend_paths import configure_backend_path


@dataclass(frozen=True)
class StellarSpectrum:
    """Wavelength in nm and flux density in erg/s/cm^2/nm."""
    wavelength_nm: np.ndarray
    flux_erg_s_cm2_nm: np.ndarray
    wavelength_medium: str = 'air'
    source: str = 'phoenix'

    def __post_init__(self):
        w = np.array(self.wavelength_nm, dtype=float, copy=True)
        f = np.array(self.flux_erg_s_cm2_nm, dtype=float, copy=True)
        if w.ndim != 1 or len(w) < 2 or f.shape != w.shape:
            raise ValueError('Stellar wavelength and flux must be matching 1-D arrays (at least two samples)')
        if not np.all(np.isfinite(w)) or np.any(w <= 0) or np.any(np.diff(w) <= 0):
            raise ValueError('Stellar wavelengths must be finite, positive and increasing')
        if not np.all(np.isfinite(f)):
            raise ValueError('Stellar flux must be finite')
        if self.wavelength_medium not in ('air', 'vacuum'):
            raise ValueError('wavelength_medium must be air or vacuum')
        w.setflags(write=False)
        f.setflags(write=False)
        object.__setattr__(self, 'wavelength_nm', w)
        object.__setattr__(self, 'flux_erg_s_cm2_nm', f)

    def plot(self, path: str | Path) -> Path:
        """Write an optional diagnostic; nothing is saved during retrieval."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(11, 4))
        ax.plot(self.wavelength_nm, self.flux_erg_s_cm2_nm, lw=0.6)
        ax.set_xlabel(f'Wavelength (nm, {self.wavelength_medium})')
        ax.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ nm$^{-1}$)')
        ax.set_title(f'Stellar spectrum — {self.source}')
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path


@dataclass(frozen=True)
class PhoenixSpectrumProvider:
    """Retrieve PHOENIX model spectra from tayph, optionally crop to a range."""
    wavelength_min_um: float | None = None
    wavelength_max_um: float | None = None
    vacuum: bool = False
    tayph_path: str | None = None

    @classmethod
    def from_config(cls, config: Config, *, config_path: str | Path = 'config') -> 'PhoenixSpectrumProvider':
        source = config['stellar_spectrum_source'].strip().lower() if 'stellar_spectrum_source' in config else 'phoenix'
        if source != 'phoenix':
            raise ValueError(f'Stellar spectrum source {source!r} not handled by PhoenixSpectrumProvider')
        # tayph config convention: air=True requests air wavelengths.
        air = config.get_bool('air', True)
        if 'instrument' in config or 'instrument_order_file' in config:
            from .instruments import instrument_wavelength_bounds_um
            lower, upper = instrument_wavelength_bounds_um(
                config, config_path)
        else:
            lower, upper = config.get_float('wl_min', 0.6), config.get_float('wl_max', 2.5)
        lo = config.get_float('stellar_wl_min', lower)
        hi = config.get_float('stellar_wl_max', upper)
        return cls(lo, hi, vacuum=not air,
                   tayph_path=config["tayph_path"] if "tayph_path" in config else None)

    def get_spectrum(self, star: Star) -> StellarSpectrum:
        if star.teff_k is None or star.logg_cgs is None:
            raise ValueError('PHOENIX requires stellar Teff and logg in config')
        if (self.wavelength_min_um is None) != (self.wavelength_max_um is None):
            raise ValueError('Specify both wavelength bounds or neither')
        if self.wavelength_min_um is not None and not (0 < self.wavelength_min_um < self.wavelength_max_um):
            raise ValueError('Invalid stellar wavelength bounds')
        if self.tayph_path is not None:
            configure_backend_path(Config.from_mapping({"tayph_path": self.tayph_path}),
                                   key="tayph_path", package="tayph")
        import astropy.units as u
        import tayph.phoenix as pho
        flux = (pho.get_phoenix_model_spectrum(star.teff_k, star.logg_cgs)
                * u.erg / u.s / u.cm**3).to_value(u.erg / u.s / u.cm**2 / u.nm)
        wave = (pho.get_phoenix_wavelengths(vacuum=self.vacuum) * u.Angstrom).to_value(u.nm)
        wave = np.asarray(wave)
        flux = np.asarray(flux)
        if self.wavelength_min_um is not None:
            select = (wave >= self.wavelength_min_um * 1000) & (wave <= self.wavelength_max_um * 1000)
            if np.count_nonzero(select) < 2:
                raise ValueError('Requested wavelength interval has fewer than two PHOENIX samples')
            wave, flux = wave[select], flux[select]
        return StellarSpectrum(wave, flux, 'vacuum' if self.vacuum else 'air', 'phoenix')
