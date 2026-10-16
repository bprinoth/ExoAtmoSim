"""Fast pipeline orchestration tests without scientific backends."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples'))
from run_full_crires_pipeline import resolve_config, run
from exoatmosim.config import Config


class FullCriresTests(unittest.TestCase):
    def test_existing_config_no_nasa_call(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / 'config'
            config.write_text('instrument CRIRES_K\nchemistry fastchem\ntellurics False\n')
            self.assertEqual(resolve_config(config, planet=None, site_config=None), config)

    def test_missing_config_and_conflicting_arguments(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                resolve_config(Path(tmp) / 'missing', planet=None, site_config=None)
            with self.assertRaises(ValueError):
                resolve_config(Path(tmp) / 'config', planet='TOI-1842 b', site_config=None)

    def test_invalid_mode_and_existing_outputs_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            planet = Path(tmp) / 'target'
            planet.mkdir()
            config = planet / 'config'
            config.write_text('instrument ESPRESSO\nchemistry fastchem\ntellurics False\n')
            with self.assertRaisesRegex(ValueError, 'CRIRES_K'):
                run(config)
            config.write_text('instrument CRIRES_K\nchemistry fastchem\ntellurics False\n')
            (planet / 'target.fits').write_text('placeholder')
            with self.assertRaises(FileExistsError):
                run(config)

    def test_tellurics_enabled_is_not_rejected_by_orchestrator(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / 'config'
            config.write_text('instrument CRIRES_K\nchemistry fastchem\ntellurics True\n')
            with self.assertRaisesRegex(ValueError, 'bulk_species'):
                run(config)


if __name__ == '__main__':
    unittest.main()
