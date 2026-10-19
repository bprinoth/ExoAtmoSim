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

    def load_earth_like(self, tp_path, chemistry_path):
        """Read an Earth analogue with unlisted N2 as the residual bulk gas.

        The complete file composition determines residual N2 first. Only
        after that are species lacking POSEIDON opacity removed from the
        radiative-transfer composition, without renormalizing other VMRs.
        Omitted species and abundances are recorded and warned about.
        """
        from ..poseidon import read_kaltenegger_pt, read_kaltenegger_chemistry, SUPPORTED_SPECIES
        atmosphere = read_kaltenegger_pt(tp_path)
        raw = read_kaltenegger_chemistry(chemistry_path)
        if 'N2' in raw:
            raise ValueError('N2 is already supplied in the chemistry file; do not add residual nitrogen twice')
        if atmosphere.altitude_km is None:
            raise ValueError('Kaltenegger T-P profile requires altitude levels')
        altitude = atmosphere.altitude_km
        values = {name: np.interp(altitude, z, vmr, left=vmr[0], right=vmr[-1])
                  for name, (z, vmr) in raw.items()}
        if not values or any(not np.all(np.isfinite(x)) or np.any(x < 0) for x in values.values()):
            raise ValueError('Invalid Kaltenegger mixing ratios')
        total = np.sum(list(values.values()), axis=0)
        if np.any(total > 1.0 + 1e-10):
            raise ValueError('Kaltenegger mixing ratios exceed unity; cannot construct residual N2')
        values['N2'] = 1.0 - total
        # The complete chemical inventory must not be re-assigned to N2.
        # Only POSEIDON-supported species are given line opacities. The
        # unmodelled species remain in a recorded inventory with their VMRs.
        unsupported = {name: v for name, v in values.items()
                       if name not in SUPPORTED_SPECIES and name != 'N2'}
        neglected = (np.sum(list(unsupported.values()), axis=0)
                     if unsupported else np.zeros_like(total))
        supported = {name: v for name, v in values.items()
                     if name in SUPPORTED_SPECIES or name == 'N2'}
        # This is NOT an opaque N2 renormalization. The RT abundance sum is
        # intentionally < 1 wherever unsupported chemistry is present.
        import warnings
        if unsupported:
            max_missing = float(np.max(neglected))
            warnings.warn(
                'Kaltenegger radiative transfer omits opacities for ' +
                ', '.join(sorted(unsupported)) +
                f'; their combined maximum VMR is {max_missing:.3e}. '
                'Their abundances are NOT transferred to N2. '
                'Review physical significance before publishing spectra.',
                UserWarning, stacklevel=2)
        self.last_inventory = {name: x.copy() for name, x in values.items()}
        self.last_omitted_opacities = tuple(sorted(unsupported))
        self.last_omitted_max_vmr = float(np.max(neglected))
        return atmosphere, AbundanceProfile(atmosphere.pressure_bar, supported, ('N2',))

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
