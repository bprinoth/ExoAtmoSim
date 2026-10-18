"""Kaltenegger/Rugheimer Earth-like atmosphere: N2 residual, opacity-filtered."""
from pathlib import Path
import warnings
import numpy as np
from ..profiles import AbundanceProfile

class KalteneggerConverter:
    def convert(self, tp_path, chemistry_path, *, bulk_species):
        from ..poseidon import read_kaltenegger_pt, read_kaltenegger_chemistry
        atmosphere = read_kaltenegger_pt(tp_path)
        chemistry = read_kaltenegger_chemistry(chemistry_path)
        values = {s: np.interp(atmosphere.altitude_km, z, x, left=x[0], right=x[-1])
                  for s, (z, x) in chemistry.items()}
        bulk = tuple(bulk_species.split()) if isinstance(bulk_species, str) else tuple(bulk_species)
        return atmosphere, AbundanceProfile(atmosphere.pressure_bar, values, bulk)

    def load_earth_like(self, tp_path, chemistry_path):
        """Fill N2 using complete inventory, THEN omit species absent from POSEIDON."""
        from ..poseidon import read_kaltenegger_pt, read_kaltenegger_chemistry, SUPPORTED_SPECIES
        atmosphere = read_kaltenegger_pt(tp_path)
        chemistry = read_kaltenegger_chemistry(chemistry_path)
        if 'N2' in chemistry:
            raise ValueError('Kaltenegger chemistry already contains N2')
        values = {s: np.interp(atmosphere.altitude_km, z, x, left=x[0], right=x[-1])
                  for s, (z, x) in chemistry.items()}
        if not values or any(not np.isfinite(x).all() or np.any(x < 0) for x in values.values()):
            raise ValueError('Invalid Kaltenegger species mixing ratios')
        summed = np.sum(list(values.values()), axis=0)
        if np.any(summed > 1 + 1e-10):
            raise ValueError('Kaltenegger volume mixing ratios exceed 1: cannot fill N2')
        values['N2'] = np.maximum(0, 1 - summed)
        omitted = {s: x for s, x in values.items() if s not in SUPPORTED_SPECIES and s != 'N2'}
        kept = {s: x for s, x in values.items() if s in SUPPORTED_SPECIES or s == 'N2'}
        omitted_sum = (np.sum(list(omitted.values()), axis=0) if omitted
                       else np.zeros_like(summed))
        self.last_inventory = {s: x.copy() for s, x in values.items()}
        self.last_omitted_opacities = tuple(sorted(omitted))
        self.last_omitted_max_vmr = float(np.max(omitted_sum))
        if omitted:
            warnings.warn('Kaltenegger POSEIDON calculation omits opacities for: ' +
                          ', '.join(sorted(omitted)) +
                          f'; maximum combined omitted VMR={self.last_omitted_max_vmr:.3e}. '
                          'N2 and supported species were not renormalized.',
                          UserWarning, stacklevel=2)
        return atmosphere, AbundanceProfile(atmosphere.pressure_bar, kept, ('N2',))

    def export(self, tp_path, chemistry_path, *, bulk_species, output_dir):
        atm, abundances = self.convert(tp_path, chemistry_path, bulk_species=bulk_species)
        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        tp_file = directory / 'tp_profile.dat'
        with tp_file.open('w', encoding='utf8') as h:
            h.write('# Temperature [K], pressure [bar]\n')
            np.savetxt(h, np.column_stack((atm.temperature_k, atm.pressure_bar)), fmt='%.16e')
        return tp_file, abundances.to_file(directory / 'abundances.dat')
