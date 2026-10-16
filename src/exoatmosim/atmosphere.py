"""Backend-independent temperature/pressure structure constructed from legacy Config.

Pressure is in bar; temperature in K; both arrays sorted increasing pressure.
The Guillot branch delegates to petitRADTRANS, matching the legacy approach.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import numpy as np

from .config import Config
from .objects import Planet


@dataclass(frozen=True)
class Atmosphere:
    pressure_bar: np.ndarray
    temperature_k: np.ndarray
    altitude_km: np.ndarray | None = None
    tp_source: str = 'unknown'

    def __post_init__(self) -> None:
        p = np.array(self.pressure_bar, dtype=float, copy=True)
        t = np.array(self.temperature_k, dtype=float, copy=True)
        if p.ndim != 1 or t.ndim != 1 or len(p) != len(t) or len(p) < 2:
            raise ValueError('pressure and temperature must be 1-D arrays with matching length >= 2')
        if not np.all(np.isfinite(p)) or np.any(p <= 0) or np.any(np.diff(p) <= 0):
            raise ValueError('pressures must be positive, finite and strictly increasing in bar')
        if not np.all(np.isfinite(t)) or np.any(t <= 0):
            raise ValueError('temperatures must be positive and finite in K')
        if self.altitude_km is not None:
            z = np.array(self.altitude_km, dtype=float, copy=True)
            if z.shape != p.shape or not np.all(np.isfinite(z)):
                raise ValueError('altitude must match pressure shape and be finite')
            z.setflags(write=False)
            object.__setattr__(self, 'altitude_km', z)
        p.setflags(write=False)
        t.setflags(write=False)
        object.__setattr__(self, 'pressure_bar', p)
        object.__setattr__(self, 'temperature_k', t)

    @classmethod
    def from_config(cls, config: Config, *, config_dir: str | Path = '.', planet: Planet | None = None) -> 'Atmosphere':
        choice = config['tp_profile'].strip() if 'tp_profile' in config else 'isothermal'
        mode = choice.lower()
        if mode in ('isothermal', 'guillot'):
            layers = config.get_int('layers', 150)
            if layers < 2:
                raise ValueError('layers must be >= 2')
            # Existing scripts used 10^-9 to 10^1 bar by default.
            lo = config.get_float('log_p_min', -9.0)
            hi = config.get_float('log_p_max', 1.0)
            if 'pressure_min' in config:
                pmin = config.get_float('pressure_min')
                if not math.isfinite(pmin) or pmin <= 0:
                    raise ValueError('pressure_min must be positive and finite [bar]')
                lo = math.log10(pmin)
            if 'pressure_max' in config:
                pmax = config.get_float('pressure_max')
                if not math.isfinite(pmax) or pmax <= 0:
                    raise ValueError('pressure_max must be positive and finite [bar]')
                hi = math.log10(pmax)
            if 'n_pressure_layers' in config:
                layers = config.get_int('n_pressure_layers')
                if layers < 2:
                    raise ValueError('n_pressure_layers must be >= 2')
            if not all(map(math.isfinite, (lo, hi))) or lo >= hi:
                raise ValueError('log_p_min must be smaller than log_p_max')
            pressure = np.logspace(lo, hi, layers)
            teq = config.get_float('Teq')
            if not math.isfinite(teq) or teq <= 0:
                raise ValueError('Teq must be positive and finite')
            if mode == 'isothermal':
                temperature = np.full(layers, config.get_float('T_iso', teq))
            else:
                if planet is None:
                    planet = Planet.from_config(config)
                logg = planet.logg_cgs if planet.logg_cgs is not None else planet.calculated_logg_cgs
                if logg is None:
                    raise ValueError('Guillot profile requires loggp or both Mp and Rp')
                try:
                    from petitRADTRANS.physics import guillot_global
                except ImportError as exc:
                    raise ImportError('Guillot profiles require petitRADTRANS; install a compatible version') from exc
                temperature = np.asarray(guillot_global(
                    pressure, config.get_float('kappa_IR', 0.01),
                    config.get_float('gamma', 0.4), 10 ** logg,
                    config.get_float('T_int', 200.0), teq), dtype=float)
            return cls(pressure, temperature, tp_source=mode)

        path = Path(choice).expanduser()
        if not path.is_absolute():
            path = Path(config_dir) / path
        # Normal file format from the legacy Planet class: T[K], P[bar].
        # Kaltenegger profiles have distinct syntax and require an explicit reader later.
        raw = np.loadtxt(path, comments='#', ndmin=2)
        if raw.shape[1] != 2:
            raise ValueError('T-P file must contain exactly two columns: temperature [K], pressure [bar]')
        order = np.argsort(raw[:, 1])
        return cls(raw[order, 1], raw[order, 0], tp_source=str(path))

    def plot_tp(self, output_path: str | Path, *, title: str | None = None) -> Path:
        """Save a T-P diagnostic; plotting is opt-in and imports matplotlib lazily."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError('T-P plots require matplotlib: pip install matplotlib') from exc
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(6, 7))
        try:
            ax.plot(self.temperature_k, self.pressure_bar, linewidth=1.8)
            ax.set_yscale('log')
            ax.invert_yaxis()
            ax.set_xlabel('Temperature (K)')
            ax.set_ylabel('Pressure (bar)')
            ax.set_title(title or f'T–P profile ({self.tp_source})')
            ax.grid(True, which='both', alpha=0.25)
            fig.tight_layout()
            fig.savefig(output_path, dpi=180)
        finally:
            plt.close(fig)
        return output_path
