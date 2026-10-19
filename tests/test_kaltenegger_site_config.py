"""Data are external; missing site configuration must explain how to obtain them."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from exoatmosim.poseidon import TransmissionRun

class KalteneggerSiteConfigTests(unittest.TestCase):
    def test_missing_site_config_guides_download(self):
        with TemporaryDirectory() as td, patch('exoatmosim.backend_paths.discover_site_config', return_value=None):
            cfg = Path(td) / 'config'
            cfg.write_text('chemistry kaltenegger\n', encoding='utf8')
            with self.assertRaisesRegex(FileNotFoundError, 'zenodo.org/records/14829587'):
                TransmissionRun.from_config(cfg)
