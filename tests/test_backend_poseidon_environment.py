import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from exoatmosim.config import Config
from exoatmosim.backend_paths import configure_poseidon_environment


class TestPoseidonEnvironment(unittest.TestCase):
    def test_configured_directories_set_both_variables(self):
        with tempfile.TemporaryDirectory() as d:
            stellar = Path(d) / 'stellar_grids'
            stellar.mkdir()
            config = Config.from_mapping({'poseidon_input_data': d, 'pysyn_cdbs': stellar})
            with patch.dict(os.environ, {}, clear=True):
                result = configure_poseidon_environment(config)
                self.assertEqual(result['POSEIDON_input_data'], str(Path(d).resolve()))
                self.assertEqual(os.environ['PYSYN_CDBS'], str(stellar.resolve()))

    def test_missing_config_preserves_existing_environment(self):
        with patch.dict(os.environ, {'POSEIDON_input_data': '/existing'}, clear=True):
            self.assertEqual(configure_poseidon_environment(Config.from_mapping({'x': 'y'})), {})
            self.assertEqual(os.environ['POSEIDON_input_data'], '/existing')

    def test_missing_directory_fails_without_side_effects(self):
        with tempfile.TemporaryDirectory() as d:
            config = Config.from_mapping({'poseidon_input_data': d,
                                          'pysyn_cdbs': str(Path(d) / 'missing')})
            with patch.dict(os.environ, {}, clear=True):
                with self.assertRaisesRegex(FileNotFoundError, 'pysyn_cdbs'):
                    configure_poseidon_environment(config)
                self.assertNotIn('POSEIDON_input_data', os.environ)


if __name__ == '__main__':
    unittest.main()
