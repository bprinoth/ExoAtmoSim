"""Fetukines chemical-network output converter (log10(VMR) columns)."""
from __future__ import annotations
from pathlib import Path
import numpy as np
from ..atmosphere import Atmosphere
from ..profiles import AbundanceProfile


class FetukinesConverter:
    """Read a column table with temperature[K], pressure[bar] and log10(VMR) species.

    Accepts a comment-prefixed header (e.g. '! altitude temperature[K] pressure[bar] ...')
    or an uncommented header. Non-species columns can be named in ignored_columns.
    """
    def convert(self, path: str | Path, *, bulk_species: str | tuple[str, ...],
                ignored_columns: tuple[str, ...] = ('altitude', 'altitude[km]', 'height', 'height[km]', 'MMW[g]')):
        lines = Path(path).read_text(encoding='utf-8').splitlines()
        header = None
        numeric_lines = []
        for raw in lines:
            stripped = raw.strip()
            if not stripped:
                continue
            if stripped.startswith(('!', '#')):
                candidate = stripped.lstrip('!#').strip().split()
                if 'temperature[K]' in candidate and 'pressure[bar]' in candidate:
                    header = candidate
                continue
            tokens = stripped.split()
            if header is None:
                if 'temperature[K]' in tokens and 'pressure[bar]' in tokens:
                    header = tokens
                    continue
                raise ValueError('Fetukines header must name temperature[K] and pressure[bar]')
            numeric_lines.append(tokens)
        if header is None or not numeric_lines:
            raise ValueError('Missing Fetukines header or data')
        if len(set(header)) != len(header):
            raise ValueError('Duplicate Fetukines header columns')
        if any(len(row) != len(header) for row in numeric_lines):
            raise ValueError('Fetukines data columns do not match header')
        try:
            data = np.asarray(numeric_lines, dtype=float)
        except ValueError as exc:
            raise ValueError('Non-numeric Fetukines data') from exc
        p = data[:, header.index('pressure[bar]')]
        t = data[:, header.index('temperature[K]')]
        order = np.argsort(p)
        atmosphere = Atmosphere(p[order], t[order], tp_source=str(path))
        ignored = set(ignored_columns) | {'temperature[K]', 'pressure[bar]'}
        species = [name for name in header if name not in ignored]
        if not species:
            raise ValueError('No Fetukines species columns')
        values = {name: np.power(10.0, data[order, header.index(name)]) for name in species}
        b = tuple(bulk_species.split()) if isinstance(bulk_species, str) else tuple(bulk_species)
        chemistry = AbundanceProfile(atmosphere.pressure_bar, values, b)
        return atmosphere, chemistry

    def export(self, path, *, bulk_species, output_dir, **kwargs):
        atmosphere, abundances = self.convert(path, bulk_species=bulk_species, **kwargs)
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        tp = output / 'tp_profile.dat'
        with tp.open('w', encoding='utf-8') as f:
            f.write('# Temperature [K], pressure [bar]\n')
            np.savetxt(f, np.column_stack((atmosphere.temperature_k, atmosphere.pressure_bar)), fmt='%.16e')
        abundances.to_file(output / 'abundances.dat')
        return tp, output / 'abundances.dat'
