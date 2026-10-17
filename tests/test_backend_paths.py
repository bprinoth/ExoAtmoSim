import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from exoatmosim import Config
from exoatmosim.backend_paths import configure_backend_path


class TestBackendPaths(unittest.TestCase):
    def test_relative_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'source'
            (root / 'POSEIDON').mkdir(parents=True)
            (root / 'POSEIDON' / '__init__.py').write_text('')
            with patch.object(sys, 'path', list(sys.path)):
                selected = configure_backend_path(Config.from_mapping({'poseidon_path': 'source'}),
                                                  key='poseidon_path', package='POSEIDON', config_dir=directory)
                self.assertEqual(selected, root)

if __name__ == '__main__':
    unittest.main()
