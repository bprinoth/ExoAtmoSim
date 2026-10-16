"""Small, dependency-free physical objects constructed from legacy configs.

Units are explicit in attribute names. These hold data; no chemistry, file I/O,
or radiative-transfer backend is imported here.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .config import Config

# IAU nominal solar and Jovian radii; Jupiter mass and CODATA G in SI.
R_SUN_M = 6.957e8
R_JUP_M = 7.1492e7
M_JUP_KG = 1.89813e27
G_SI = 6.67430e-11


def _optional_float(config: Config, key: str) -> float | None:
    if key not in config:
        return None
    value = config.get_float(key)
    return value if math.isfinite(value) else None


def _positive(name: str, value: float | None) -> None:
    if value is not None and (not math.isfinite(value) or value <= 0):
        raise ValueError(f'{name} must be positive and finite')


@dataclass(frozen=True)
class Star:
    radius_solar: float | None = None
    mass_solar: float | None = None
    teff_k: float | None = None
    logg_cgs: float | None = None
    spectral_type: str | None = None

    def __post_init__(self):
        for name in ('radius_solar', 'mass_solar', 'teff_k'):
            _positive(name, getattr(self, name))

    @classmethod
    def from_config(cls, config: Config) -> Star:
        return cls(radius_solar=_optional_float(config, 'Rs'),
                   mass_solar=_optional_float(config, 'Ms'),
                   teff_k=_optional_float(config, 'Teff'),
                   logg_cgs=_optional_float(config, 'logg'),
                   spectral_type=config['stellar_type'] if 'stellar_type' in config else None)


@dataclass(frozen=True)
class Planet:
    radius_jupiter: float | None = None
    mass_jupiter: float | None = None
    logg_cgs: float | None = None
    equilibrium_temperature_k: float | None = None
    reference_pressure_bar: float | None = None

    def __post_init__(self):
        for name in ('radius_jupiter', 'mass_jupiter', 'equilibrium_temperature_k',
                     'reference_pressure_bar'):
            _positive(name, getattr(self, name))
        if self.logg_cgs is not None and not math.isfinite(self.logg_cgs):
            raise ValueError('logg_cgs must be finite')

    @property
    def calculated_logg_cgs(self) -> float | None:
        """Gravity calculated independently from mass and radius, if available."""
        if self.mass_jupiter is None or self.radius_jupiter is None:
            return None
        g_m_s2 = G_SI * self.mass_jupiter * M_JUP_KG / (self.radius_jupiter * R_JUP_M)**2
        return math.log10(g_m_s2 * 100.)

    @property
    def gravity_discrepancy_dex(self) -> float | None:
        """Configured logg minus gravity inferred from mass/radius; never silently overwritten."""
        calculated = self.calculated_logg_cgs
        if calculated is None or self.logg_cgs is None:
            return None
        return self.logg_cgs - calculated

    @classmethod
    def from_config(cls, config: Config) -> Planet:
        return cls(radius_jupiter=_optional_float(config, 'Rp'),
                   mass_jupiter=_optional_float(config, 'Mp'),
                   logg_cgs=_optional_float(config, 'loggp'),
                   equilibrium_temperature_k=_optional_float(config, 'Teq'),
                   reference_pressure_bar=_optional_float(config, 'p0'))
