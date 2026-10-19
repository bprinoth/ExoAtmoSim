"""Regression checks for Earth-like residual nitrogen and NASA type provenance."""
import unittest
import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import numpy as np

from exoatmosim.config import Config
from exoatmosim.nasa import parameters_from_row
from exoatmosim.converters.kaltenegger import KalteneggerConverter
from exoatmosim.poseidon import select_legacy_profiles, TransmissionRun


class KalteneggerNitrogenTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("astropy"), "astropy required")
    def test_nasa_spectral_type_is_not_mapped_during_fetch(self):
        self.assertEqual(parameters_from_row({'st_spectype': 'M2.5V'})['stellar_type'], 'M2.5V')

    def test_mapping_occurs_only_when_profiles_selected(self):
        with TemporaryDirectory() as d:
            root = Path(d)
            (root / 'PT_profiles').mkdir()
            (root / 'chem_profiles').mkdir()
            (root / 'PT_profiles' / 'AltPT_M3A_0.0Ga.txt').touch()
            (root / 'chem_profiles' / 'Chem_M3A_0.0Ga.txt').touch()
            c = Config.from_mapping({'stellar_type':'M2.5V', 'Gyr_age':'4.6'})
            t,chem=select_legacy_profiles(c,root)
            self.assertIn('M3A',t.name)
            self.assertEqual(c['stellar_type'],'M2.5V')

    def test_residual_n2_with_uploaded_profile(self):
        root = Path(__file__).parent/'fixtures'/'kaltenegger'
        with self.assertWarns(UserWarning):
            atmosphere, chemistry=KalteneggerConverter().load_earth_like(
                root/'AltPT_F0V_0.0Ga.txt', root/'Chem_F0V_0.0Ga.txt')
        self.assertEqual(chemistry.bulk_species, ('N2',))
        self.assertEqual(len(atmosphere.pressure_bar), 100)
        self.assertTrue(np.allclose(sum(chemistry.mixing_ratios.values()),1.,atol=1e-10))
        self.assertTrue(np.all((chemistry.mixing_ratios['N2'] > 0.78) &
                               (chemistry.mixing_ratios['N2'] < 0.80)))
        self.assertTrue(np.allclose(chemistry.mixing_ratios['O2'], .21))

    def test_rejects_overfilled_chemistry(self):
        root = Path(__file__).parent/'fixtures'/'kaltenegger'
        with TemporaryDirectory() as d:
            src=(root/'Chem_F0V_0.0Ga.txt').read_text()
            changed=src.replace('0.21E+00','1.10E+00')
            target=Path(d)/'chem.txt';target.write_text(changed)
            with self.assertRaisesRegex(ValueError,'exceed unity'):
                KalteneggerConverter().load_earth_like(root/'AltPT_F0V_0.0Ga.txt',target)

    def test_transmission_run_consumes_kaltenegger_without_fastchem(self):
        root = Path(__file__).parent/'fixtures'/'kaltenegger'
        with TemporaryDirectory() as d:
            dest=Path(d)
            config=Config.from_mapping({
                'stellar_type':'F0V', 'Gyr_age':'4.6', 'chemistry':'kaltenegger',
                'kaltenegger_profiles_dir':str(dest),
                'Rp':'1', 'Mp':'1', 'Teq':'280', 'Rs':'1', 'Teff':'7200',
                'logg':'4.3', 'instrument':'CRIRES_K'
            })
            (dest/'PT_profiles').mkdir();(dest/'chem_profiles').mkdir()
            import shutil
            shutil.copy(root/'AltPT_F0V_0.0Ga.txt', dest/'PT_profiles')
            shutil.copy(root/'Chem_F0V_0.0Ga.txt', dest/'chem_profiles')
            cfg_path=dest/'config';cfg_path.write_text(config.to_text())
            with self.assertWarns(UserWarning):
                out=TransmissionRun.from_config(cfg_path)
            self.assertEqual(out.abundances.bulk_species,('N2',))
            self.assertEqual(out.config['stellar_type'],'F0V')


class KalteneggerAutomaticSelectionTests(unittest.TestCase):
    def test_all_available_ages_are_matched_and_nasa_values_preserved(self):
        from exoatmosim.poseidon import select_kaltenegger_profiles
        with TemporaryDirectory() as d:
            root = Path(d)
            (root/'PT_profiles').mkdir()
            (root/'chem_profiles').mkdir()
            for age in ('0.0', '0.8', '2.0', '3.9'):
                (root/'PT_profiles'/f'AltPT_M3A_{age}Ga.txt').touch()
                (root/'chem_profiles'/f'Chem_M3A_{age}Ga.txt').touch()
            for stellar_age, expected in [(4.6,'0.0'),(3.8,'0.8'),(2.6,'2.0'),(0.7,'3.9')]:
                config = Config.from_mapping({'stellar_type':'M2.5V','Gyr_age':str(stellar_age)})
                pt, chem = select_kaltenegger_profiles(config, root)
                self.assertEqual(pt.name, f'AltPT_M3A_{expected}Ga.txt')
                self.assertEqual(chem.name, f'Chem_M3A_{expected}Ga.txt')
                self.assertEqual(config['stellar_type'],'M2.5V')
                self.assertEqual(config['Gyr_age'],str(stellar_age))

    def test_missing_age_is_not_assumed_present_day(self):
        from exoatmosim.poseidon import select_kaltenegger_profiles
        with TemporaryDirectory() as root:
            with self.assertRaisesRegex(ValueError,'Gyr_age'):
                select_kaltenegger_profiles(Config.from_mapping({'stellar_type':'G2V'}),root)

    def test_solar_class_maps_to_sun(self):
        from exoatmosim.nasa import closest_spectral_type
        self.assertEqual(closest_spectral_type('G2V'), 'Sun')
