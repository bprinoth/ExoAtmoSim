import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from exoatmosim.backend_paths import discover_site_config
from exoatmosim.observing.telluric_grid import resolve_grid_file


class TestTelluricSitePaths(unittest.TestCase):
    def test_discovery_from_planet_checkout_without_cwd(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'source'
            planet = root / 'data' / 'WASP-189_b' / 'config'
            planet.parent.mkdir(parents=True)
            planet.write_text('instrument ESPRESSO\n')
            settings = root / 'config' / 'site_paths.conf'
            settings.parent.mkdir()
            grids = root / 'private_grids'
            grids.mkdir()
            expected = grids / 'telluric_paranal.h5'
            expected.touch()
            settings.write_text('telluric_grid_dir ../private_grids\n')
            with patch.dict(os.environ, {'EXOATMOSIM_SITE_CONFIG': ''}):
                with patch('exoatmosim.backend_paths.Path.cwd', return_value=Path(tmp)):
                    self.assertEqual(discover_site_config(config_path=planet), settings)
                    self.assertEqual(resolve_grid_file('paranal', config_path=planet), expected)

    def test_env_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'grid').mkdir()
            target = root/'grid'/'telluric_mauna_kea.h5'
            target.touch()
            conf = root/'site_paths.conf'
            conf.write_text('telluric_grid_dir grid\n')
            with patch.dict(os.environ, {'EXOATMOSIM_SITE_CONFIG': str(conf)}):
                self.assertEqual(resolve_grid_file('mauna_kea'), target)

    def test_missing_grid_is_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            conf = Path(tmp)/'site_paths.conf'
            conf.write_text('telluric_grid_dir missing\n')
            with patch.dict(os.environ, {'EXOATMOSIM_SITE_CONFIG': str(conf)}):
                with self.assertRaises(FileNotFoundError):
                    resolve_grid_file('paranal')


if __name__ == '__main__':
    unittest.main()
