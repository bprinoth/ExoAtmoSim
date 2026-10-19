"""Regression coverage for external Kaltenegger profiles without third-party data."""
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
from exoatmosim.config import Config
from exoatmosim.poseidon import select_legacy_profiles, TransmissionRun
from exoatmosim.converters.kaltenegger import KalteneggerConverter

class KalteneggerIntegrationTests(unittest.TestCase):
    def test_missing_local_data_guidance(self):
        with tempfile.TemporaryDirectory() as d, patch('exoatmosim.backend_paths.discover_site_config', return_value=None):
            cfg = Path(d) / 'config'
            cfg.write_text('chemistry kaltenegger\n', encoding='utf8')
            with self.assertRaisesRegex(FileNotFoundError, 'zenodo.org'):
                TransmissionRun.from_config(cfg)

    def test_preserved_nasa_and_epoch_choice(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for sub in ('PT_profiles', 'chem_profiles'):
                (root/sub).mkdir()
            for age in ('0.0','0.8','2.0','3.9'):
                (root/'PT_profiles'/f'AltPT_Sun_{age}Ga.txt').write_text('x')
                (root/'chem_profiles'/f'Chem_Sun_{age}Ga.txt').write_text('x')
            for stellar_age, wanted in ((4.6,'0.0'),(3.8,'0.8'),(2.6,'2.0'),(.7,'3.9')):
                cfg = Config.from_text(f'stellar_type G2V\nGyr_age {stellar_age}\n')
                pt, chem = select_legacy_profiles(cfg,root)
                self.assertEqual(pt.name,f'AltPT_Sun_{wanted}Ga.txt')
                self.assertEqual(cfg['stellar_type'],'G2V')

    def test_n2_filled_before_opacity_filtering(self):
        with tempfile.TemporaryDirectory() as d:
            pt=Path(d)/'pt.txt'
            chem=Path(d)/'chem.txt'
            pt.write_text('Z,T,P 0 300 1000\nZ,T,P 1 290 900\n')
            chem.write_text('Molecule O2\nO2_MR .21,.21\nO2_ALT 0,1\n'
                            'Molecule H2O2\nH2O2_MR .00002,.00003\nH2O2_ALT 0,1\n')
            converter=KalteneggerConverter()
            with self.assertWarns(UserWarning):
                atmosphere, abundances=converter.load_earth_like(pt,chem)
            np.testing.assert_allclose(np.sort(abundances.mixing_ratios['N2']),[.78997,.78998])
            self.assertNotIn('H2O2',abundances.mixing_ratios)
            self.assertIn('H2O2',converter.last_inventory)

    def test_local_site_path_is_allowed_but_not_copied_to_planet_config(self):
        from exoatmosim.nasa import load_site_config, format_science_config
        with tempfile.TemporaryDirectory() as d:
            site_file = Path(d) / 'site_paths.conf'
            site_file.write_text(
                'kaltenegger_profiles_dir /tmp/downloaded_kaltenegger\n',
                encoding='utf8')
            site = load_site_config(site_file)
            self.assertEqual(site['kaltenegger_profiles_dir'],
                             '/tmp/downloaded_kaltenegger')
            result = format_science_config({'stellar_type': 'G2V', 'Gyr_age': 4.6},
                                           site_values=site, include_workflow=False)
            self.assertNotIn('kaltenegger_profiles_dir', result)
            self.assertIn('stellar_type', result)
