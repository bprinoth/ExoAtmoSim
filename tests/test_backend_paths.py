import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from exoatmosim import Config
from exoatmosim.backend_paths import configure_backend_path


class TestBackendPaths(unittest.TestCase):
    def test_absent_key_is_noop(self):
        with patch.object(sys, 'path', list(sys.path)) as paths:
            self.assertIsNone(configure_backend_path(Config.from_mapping({'air': True}),
                                                     key='tayph_path', package='tayph'))
            self.assertEqual(paths, sys.path)

    def test_explicit_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pkg = root / 'tayph'
            pkg.mkdir()
            (pkg / '__init__.py').write_text('')
            with patch.object(sys, 'path', list(sys.path)):
                selected = configure_backend_path(Config.from_mapping({'tayph_path': str(root)}),
                                                  key='tayph_path', package='tayph')
                self.assertEqual(selected, root)
                self.assertEqual(sys.path[0], str(root))

    def test_relative_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'source'
            (root / 'POSEIDON').mkdir(parents=True)
            (root / 'POSEIDON' / '__init__.py').write_text('')
            with patch.object(sys, 'path', list(sys.path)):
                selected = configure_backend_path(Config.from_mapping({'poseidon_path': 'source'}),
                                                  key='poseidon_path', package='POSEIDON', config_dir=directory)
                self.assertEqual(selected, root)

    def test_missing_checkout_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                configure_backend_path(Config.from_mapping({'tayph_path': directory}),
                                       key='tayph_path', package='tayph')


if __name__ == '__main__':
    unittest.main()
