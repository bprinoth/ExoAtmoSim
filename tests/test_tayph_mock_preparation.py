import tempfile
import unittest
from pathlib import Path

from exoatmosim.tayph_mock import prepare_tayph_mock


class TestTayphPreparation(unittest.TestCase):
    def test_wasp189b_staging_idempotent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            planet_dir = root / 'data' / 'WASP-189_b'
            planet_dir.mkdir(parents=True)
            cfg = planet_dir / 'config'
            cfg.write_text('instrument\tESPRESSO\nP\t2.724033\n')
            prepared = prepare_tayph_mock(cfg, project_root=root)
            self.assertEqual(prepared.runfile.read_text(),
                             'datapath\tdata/WASP-189_b/\n'
                             'mockdatapath\tdata/WASP-189_b/\n'
                             'mockdata\tWASP-189_b.fits\n'
                             'transit\tTrue\n')
            self.assertTrue(prepared.order_definition.is_file())
            self.assertTrue(prepared.blaze.is_file())
            self.assertFalse(prepared.expected_mock_fits.exists())
            prepare_tayph_mock(cfg, project_root=root)
            self.assertEqual(cfg.read_text(), 'instrument\tESPRESSO\nP\t2.724033\n')

    def test_does_not_overwrite_user_resource(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            planet_dir = root / 'data' / 'WASP-189_b'
            planet_dir.mkdir(parents=True)
            cfg = planet_dir / 'config'
            cfg.write_text('instrument\tESPRESSO\n')
            (planet_dir / 'blaze_ESPRESSO.npy').write_bytes(b'user version')
            with self.assertRaises(FileExistsError):
                prepare_tayph_mock(cfg, project_root=root)
