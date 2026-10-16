import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from exoatmosim.tayph_mock import prepare_tayph_mock, run_tayph_mock
from exoatmosim.instruments import available_instruments


class MultiInstrumentStaging(unittest.TestCase):
    def config(self, root, instrument):
        planet = Path(root) / 'data' / 'TOI-1842_b'
        planet.mkdir(parents=True)
        cfg = planet / 'config'
        cfg.write_text(f'instrument {instrument}\n')
        return planet, cfg

    def test_bundled_resource_status(self):
        resources = available_instruments()
        self.assertTrue(resources['CRIRES_K']['orders'])
        self.assertFalse(resources['CRIRES_K']['blaze'])
        self.assertTrue(resources['SPIROU']['blaze'])

    def test_stages_all_available_order_grids(self):
        for inst in ('ESPRESSO', 'SPIROU', 'CRIRES_K', 'IGRINS_HK', 'IGRINS2_HK', 'ANDES_K'):
            with self.subTest(inst=inst), tempfile.TemporaryDirectory() as td:
                d, cfg = self.config(td, inst)
                result = prepare_tayph_mock(cfg, project_root=td)
                self.assertEqual(result.order_definition.name, f'order_def_{inst}.npy')
                self.assertTrue(result.order_definition.is_file())
                if inst in ('ESPRESSO', 'SPIROU'):
                    self.assertTrue(result.blaze.is_file())
                    self.assertEqual(np.load(result.blaze).shape, np.load(result.order_definition).shape)
                else:
                    self.assertIsNone(result.blaze)
                    self.assertFalse((d / f'blaze_{inst}.npy').exists())
                self.assertEqual(result.runfile.parent, d)
                prepare_tayph_mock(cfg, project_root=td)

    def test_missing_blaze_is_left_to_tayph(self):
        with tempfile.TemporaryDirectory() as td:
            d,cfg = self.config(td, 'CRIRES_K')
            (d / 'TOI-1842_b.fits').write_bytes(b'placeholder')
            result,kw = run_tayph_mock(cfg, project_root=td, dry_run=True)
            self.assertIsNone(result.blaze)
            self.assertEqual(kw['mode'], 'CRIRES_K')

    def test_explicit_crires_blaze_and_shape_check(self):
        with tempfile.TemporaryDirectory() as td:
            d,cfg = self.config(td, 'CRIRES_K')
            source = prepare_tayph_mock(cfg, project_root=td).order_definition
            shape = np.load(source).shape
            result_blaze = d / 'blaze_CRIRES_K.npy'
            custom = d / 'real_crires_blaze.npy'
            np.save(custom, np.ones((2,3)))
            cfg.write_text('instrument CRIRES_K\ninstrument_blaze_file real_crires_blaze.npy\n')
            with self.assertRaisesRegex(ValueError, 'blaze shape'):
                prepare_tayph_mock(cfg, project_root=td)
            np.save(custom, np.ones(shape))
            result = prepare_tayph_mock(cfg, project_root=td)
            self.assertTrue(result.blaze.is_file())
            self.assertEqual(result.blaze.name, 'blaze_CRIRES_K.npy')

    def test_user_resources_are_preserved(self):
        with tempfile.TemporaryDirectory() as td:
            d,cfg = self.config(td, 'SPIROU')
            (d/'order_def_SPIROU.npy').write_bytes(b'custom')
            with self.assertRaises(FileExistsError):
                prepare_tayph_mock(cfg, project_root=td)

if __name__ == '__main__':
    unittest.main()
