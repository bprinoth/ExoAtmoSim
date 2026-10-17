"""Backend-independent abundance tables: pressure [bar], species linear VMR."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from .atmosphere import Atmosphere


@dataclass(frozen=True)
class AbundanceProfile:
    pressure_bar: np.ndarray
    mixing_ratios: dict[str, np.ndarray]
    bulk_species: tuple[str, ...]

    def __post_init__(self):
        p = np.array(self.pressure_bar, dtype=float, copy=True)
        if p.ndim != 1 or len(p) < 2 or not np.all(np.isfinite(p)) or np.any(p <= 0) or np.any(np.diff(p) <= 0):
            raise ValueError('pressure must be positive, finite and strictly increasing (bar)')
        names = tuple(self.mixing_ratios)
        if not names or len(names) != len(set(names)):
            raise ValueError('At least one unique species required')
        bulk = tuple(self.bulk_species)
        if not bulk or len(bulk) != len(set(bulk)) or any(s not in names for s in bulk):
            raise ValueError('bulk_species must explicitly name unique species present in abundance data')
        checked = {}
        for name, values in self.mixing_ratios.items():
            a = np.array(values, dtype=float, copy=True)
            if a.shape != p.shape or not np.all(np.isfinite(a)) or np.any(a < 0) or np.any(a > 1):
                raise ValueError(f'invalid linear VMR profile for {name}')
            a.setflags(write=False)
            checked[name] = a
        total = np.sum(list(checked.values()), axis=0)
        if np.any(total > 1 + 1e-8):
            raise ValueError('sum of volume mixing ratios exceeds one')
        p.setflags(write=False)
        object.__setattr__(self, 'pressure_bar', p)
        object.__setattr__(self, 'mixing_ratios', checked)
        object.__setattr__(self, 'bulk_species', bulk)

    @classmethod
    def from_file(cls, path: str | Path, *, bulk_species: str | tuple[str, ...] | list[str]):
        """Header: pressure H2 He H2O ...; comments start with # or !."""
        lines = []
        for line in Path(path).read_text(encoding='utf-8').splitlines():
            line = line.split('#', 1)[0].split('!', 1)[0].strip()
            if line:
                lines.append(line)
        if len(lines) < 3:
            raise ValueError('abundance table needs a header and at least two data rows')
        names = lines[0].split()
        if names[0].lower() != 'pressure' or len(names) < 2 or len(set(names)) != len(names):
            raise ValueError('header must be: pressure speciesA speciesB ... (unique names)')
        try:
            data = np.asarray([[float(x) for x in line.split()] for line in lines[1:]], dtype=float)
        except ValueError as exc:
            raise ValueError('abundance rows must contain numeric linear VMRs') from exc
        if data.ndim != 2 or data.shape[1] != len(names):
            raise ValueError('row column count does not match abundance header')
        order = np.argsort(data[:, 0])
        b = tuple(bulk_species.split()) if isinstance(bulk_species, str) else tuple(bulk_species)
        return cls(data[order, 0], {name: data[order, idx] for idx, name in enumerate(names[1:], 1)}, b)

    def to_file(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        names = list(self.mixing_ratios)
        matrix = np.column_stack([self.pressure_bar] + [self.mixing_ratios[n] for n in names])
        with path.open('w', encoding='utf-8') as f:
            f.write('# Pressure [bar], abundances: linear volume mixing ratios\n')
            f.write('pressure ' + ' '.join(names) + '\n')
            np.savetxt(f, matrix, fmt='%.16e')
        return path

    def on_atmosphere(self, atmosphere: Atmosphere, *, extrapolate: bool = False) -> 'AbundanceProfile':
        """Linear VMR interpolation versus log10(P), no silent extrapolation."""
        target = atmosphere.pressure_bar
        p = self.pressure_bar
        if not extrapolate and (target[0] < p[0] * (1 - 1e-12) or target[-1] > p[-1] * (1 + 1e-12)):
            raise ValueError('T-P pressure grid extends beyond abundance grid; explicit extrapolate=True required')
        x = np.log10(p)
        out = {s: np.interp(np.log10(target), x, vmr, left=vmr[0], right=vmr[-1])
               for s, vmr in self.mixing_ratios.items()}
        return AbundanceProfile(target, out, self.bulk_species)

    def as_poseidon_composition(self):
        """Bridge to existing experimental adapter; species ordering belongs to backend."""
        from .poseidon import Composition
        return Composition(self.mixing_ratios)

    def plot(self, output_path: str | Path, *, title: str = 'Atmospheric abundances') -> Path:
        """Opt-in log-VMR/pressure diagnostic; zero values are not plotted."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 7))
        try:
            for species, abundance in self.mixing_ratios.items():
                valid = abundance > 0
                if np.any(valid):
                    ax.plot(abundance[valid], self.pressure_bar[valid], label=species)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.invert_yaxis()
            ax.set_xlabel('Linear volume mixing ratio (VMR)')
            ax.set_ylabel('Pressure (bar)')
            ax.set_title(title)
            ax.grid(True, which='both', alpha=0.25)
            ax.legend()
            fig.tight_layout()
            fig.savefig(path, dpi=180)
        finally:
            plt.close(fig)
        return path
