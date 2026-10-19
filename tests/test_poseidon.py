import tempfile
import unittest
from pathlib import Path
import numpy as np

from exoatmosim import Atmosphere, Config, Star, Planet
from exoatmosim.poseidon import (Composition, read_kaltenegger_pt,
    read_kaltenegger_chemistry, composition_on_atmosphere,
    select_legacy_profiles, PoseidonTransmissionModel, validate_hminus_continuum)


class TestPoseidonBoundary(unittest.TestCase):
    def test_kaltenegger_pt_and_chemistry(self):
        with tempfile.TemporaryDirectory() as root:
            pt = Path(root) / 'pt.txt'
            chem = Path(root) / 'chem.txt'
            pt.write_text('! header\nZ,T,P 0 300 1000\nZ,T,P 10 250 100\n')
            chem.write_text('Molecule N2\nN2_MR 0.8, ^\n 0.7\nN2_ALT 0, 10\nMolecule O2\nO2_MR 0.2, 0.3\nO2_ALT 0, 10\n')
            a = read_kaltenegger_pt(pt)
            self.assertTrue(np.allclose(a.pressure_bar, [0.1, 1.0]))
            self.assertTrue(np.allclose(a.altitude_km, [10, 0]))
            c = composition_on_atmosphere(read_kaltenegger_chemistry(chem), a)
            np.testing.assert_allclose(c.mixing_ratios['N2'], [0.7, 0.8])

    def test_profile_selection(self):
        with tempfile.TemporaryDirectory() as root:
            r = Path(root)
            for directory, name in [('PT_profiles', 'AltPT_M3A_0.0Ga.txt'),
                                    ('chem_profiles', 'Chem_M3A_0.0Ga.txt')]:
                d = r / directory
                d.mkdir()
                (d / name).write_text('test')
            config = Config.from_mapping({'stellar_type': 'M2.5V', 'Gyr_age': '4.6'})
            pt, chem = select_legacy_profiles(config, r)
            self.assertIn('M3A', pt.name)
            self.assertIn('M3A', chem.name)

    def test_backend_validation_without_poseidon(self):
        config = Config.from_mapping({'Rs': 1, 'Teff': 5700, 'logg': 4.4,
            'Rp': 1.0, 'Mp': 1.0, 'Teq': 1200, 'resolution': 300000})
        atmosphere = Atmosphere(np.array([0.1, 1.0]), np.array([1000, 1100]))
        model = PoseidonTransmissionModel.from_config(config)
        good = Composition({'O2': [0.2, 0.3], 'H2O': [0.01, 0.02]})
        model.validate(atmosphere, good, 'O2', Planet.from_config(config), Star.from_config(config))
        with self.assertRaisesRegex(ValueError, 'Species'):
            model.validate(atmosphere, Composition({'O2': [0.2, 0.3], 'FAKE': [0, 0]}),
                           'O2', Planet.from_config(config), Star.from_config(config))

    def test_hydrogen_helium_and_hminus_fastchem_mix(self):
        config = Config.from_mapping({'Rs': 1, 'Teff': 8000, 'logg': 4.1,
            'Rp': 1.5, 'Mp': 2.0, 'Teq': 2500})
        atmosphere = Atmosphere(np.array([1e-15, 1e-6, 1.0]),
                                np.array([3000, 2600, 2200]))
        chemistry = Composition({
            'H2': [1e-7, 0.1, 0.7], 'He': [0.1, 0.1, 0.1],
            'H': [0.8, 0.7, 0.1], 'H-': [1e-12] * 3,
            'e-': [1e-5] * 3, 'Fe': [1e-6] * 3,
            'Fe+': [1e-7] * 3, 'Ti': [1e-8] * 3,
            'Ti+': [1e-9] * 3,
        })
        PoseidonTransmissionModel().validate(
            atmosphere, chemistry, ('H2', 'He'),
            Planet.from_config(config), Star.from_config(config))

    def test_hminus_continuum_registration(self):
        chem = {'H2': [0.7], 'He': [0.1], 'H': [0.1],
                'H-': [1e-12], 'e-': [1e-5], 'Fe': [1e-5]}
        validate_hminus_continuum({'bf_species': ['H-bf'], 'ff_pairs': ['H-ff']}, chem)
        with self.assertRaisesRegex(RuntimeError, 'H-bf'):
            validate_hminus_continuum({'bf_species': [], 'ff_pairs': ['H-ff']}, chem)
        with self.assertRaisesRegex(RuntimeError, 'H-ff'):
            validate_hminus_continuum({'bf_species': ['H-bf'], 'ff_pairs': []}, chem)
        with self.assertRaisesRegex(ValueError, 'requires H'):
            validate_hminus_continuum({'bf_species': ['H-bf'], 'ff_pairs': []}, {'H-': [1e-12]})

    def test_wrong_layer_count(self):
        config = Config.from_mapping({'Rs': 1, 'Teff': 5700, 'logg': 4.4,
            'Rp': 1, 'Mp': 1, 'Teq': 1000})
        with self.assertRaisesRegex(ValueError, 'layer count'):
            PoseidonTransmissionModel().validate(
                Atmosphere([0.1, 1], [900, 1000]), Composition({'O2': [0.2, 0.3, 0.4]}),
                'O2', Planet.from_config(config), Star.from_config(config))


if __name__ == '__main__':
    unittest.main()
