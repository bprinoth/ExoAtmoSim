"""Read original Kaltenegger altitude profiles and translate to pressure VMR table."""
from __future__ import annotations
from pathlib import Path
import numpy as np
from ..atmosphere import Atmosphere
from ..profiles import AbundanceProfile


class KalteneggerConverter:
    def convert(self, tp_path: str | Path, chemistry_path: str | Path, *, bulk_species: str | tuple[str, ...]):
        from ..poseidon import read_kaltenegger_pt, read_kaltenegger_chemistry
        atmosphere = read_kaltenegger_pt(tp_path)
        chemistry = read_kaltenegger_chemistry(chemistry_path)
        if atmosphere.altitude_km is None:
            raise ValueError('Kaltenegger profiles require altitude grid')
        values = {name: np.interp(atmosphere.altitude_km, altitude, vmr,
                                  left=vmr[0], right=vmr[-1])
                  for name, (altitude, vmr) in chemistry.items()}
        return atmosphere, AbundanceProfile(atmosphere.pressure_bar, values,
                                            tuple(bulk_species.split()) if isinstance(bulk_species, str) else tuple(bulk_species))

    def export(self, tp_path, chemistry_path, *, bulk_species, output_dir):
        atm, chem = self.convert(tp_path, chemistry_path, bulk_species=bulk_species)
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        tp = output / 'tp_profile.dat'
        with tp.open('w', encoding='utf-8') as f:
            f.write('# Temperature [K], pressure [bar]\n')
            np.savetxt(f, np.column_stack((atm.temperature_k, atm.pressure_bar)), fmt='%.16e')
        chem.to_file(output / 'abundances.dat')
        return tp, output / 'abundances.dat'
