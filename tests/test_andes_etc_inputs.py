"""Manual ANDES ETC input-folder regression tests (no external ETC access)."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from exoatmosim.andes_etc_inputs import scaffold_andes_etc
from exoatmosim import workflow


class AndesETCInputTests(unittest.TestCase):
    def test_scaffold_creates_nine_slots_without_overwriting_files(self):
        with tempfile.TemporaryDirectory() as td:
            planet = Path(td) / 'TRAPPIST-1_b'
            root = scaffold_andes_etc(planet)
            expected = {root / f'airmass_{am}' / f'{arm}.csv'
                        for am in ('1.0', '1.5', '3.0')
                        for arm in ('UBV', 'RIZ', 'YJH')}
            self.assertEqual(len(expected), 9)
            self.assertTrue(all(path.parent.is_dir() for path in expected))
            self.assertTrue(all(not path.exists() for path in expected))
            self.assertIn('exposure_time', (root / 'README.md').read_text())
            example = root / 'airmass_1.0' / 'YJH.csv'
            example.write_text('example data')
            scaffold_andes_etc(planet)
            self.assertEqual(example.read_text(), 'example data')

    def test_armazones_location_uses_config_coordinates(self):
        from exoatmosim.config import Config
        from exoatmosim.observing.scheduling import _sky_location
        cfg = Config.from_text(
            'RA 23:06:29.28\\nDEC -05:02:29\\n'
            'observatory armazones\\n'
            'lat -24.5892\\nlong -70.1916\\nelev 3046\\n')
        _, location = _sky_location(cfg)
        self.assertAlmostEqual(location.lon.deg, -70.1916, places=5)
        self.assertAlmostEqual(location.lat.deg, -24.5892, places=5)
        self.assertAlmostEqual(location.height.to_value('m'), 3046, delta=2)

    def test_andes_site_defaults_are_applied_to_archive_rows(self):
        from unittest.mock import MagicMock, patch
        from exoatmosim.nasa import find_parameters
        from exoatmosim.config import Config
        # Do not query NASA; exercise the real config-writing function.
        import sys
        from types import ModuleType
        package = ModuleType('astroquery')
        ipac = ModuleType('astroquery.ipac')
        nexsci = ModuleType('astroquery.ipac.nexsci')
        archive = ModuleType('astroquery.ipac.nexsci.nasa_exoplanet_archive')
        archive.NasaExoplanetArchive = MagicMock()
        archive.NasaExoplanetArchive.query_criteria.return_value = [object()]
        with tempfile.TemporaryDirectory() as td:
            with patch.dict(sys.modules, {'astroquery': package,
                                          'astroquery.ipac': ipac,
                                          'astroquery.ipac.nexsci': nexsci,
                                          'astroquery.ipac.nexsci.nasa_exoplanet_archive': archive}), \\
                 patch('exoatmosim.nasa.parameters_from_row',
                       return_value={'RA':'23:06:29', 'DEC':'-05:02:29',
                                     'lat':-24.6276,'long':-70.4051,'elev':2635}), \\
                 patch('exoatmosim.instruments.instrument_resource'):
                find_parameters('TRAPPIST-1 b', data_dir=td, instrument='ANDES',
                                print_params=False)
            cfg = Config.from_file(Path(td) / 'TRAPPIST-1_b' / 'config')
            self.assertEqual(cfg['observatory'], 'armazones')
            self.assertAlmostEqual(cfg.get_float('lat'), -24.5892)
            self.assertAlmostEqual(cfg.get_float('long'), -70.1916)
            self.assertAlmostEqual(cfg.get_float('elev'), 3046)

    def test_make_config_scaffolds_only_andes(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'planet' / 'config'
            path.parent.mkdir()
            with patch('exoatmosim.workflow.resolve_config', return_value=path):
                path.write_text('instrument ANDES\nexposure_time 60\n')
                with redirect_stdout(StringIO()):
                    self.assertEqual(workflow.make_config('planet', instrument='ANDES'), path)
                self.assertTrue((path.parent / 'andes_etc' / 'airmass_3.0').is_dir())
                path.write_text('instrument CRIRES_K\nexposure_time 60\n')
                import shutil
                shutil.rmtree(path.parent / 'andes_etc')
                with redirect_stdout(StringIO()):
                    workflow.make_config('planet', instrument='CRIRES_K')
                self.assertFalse((path.parent / 'andes_etc').exists())


if __name__ == '__main__':
    unittest.main()
