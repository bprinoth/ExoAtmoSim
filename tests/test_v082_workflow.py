import unittest
import numpy as np
from exoatmosim.config import Config
from exoatmosim.atmosphere import Atmosphere
from exoatmosim.nasa import format_science_config
from exoatmosim.instruments import instrument_resource
from exoatmosim.observations import SpectralOrders

class TestWorkflowV082(unittest.TestCase):
    def test_wasp_pressure_range(self):
        cfg = Config.from_text('Teq 2500\ntp_profile isothermal\npressure_min 1e-15\npressure_max 10\nn_pressure_layers 200\n')
        atm = Atmosphere.from_config(cfg)
        self.assertEqual(len(atm.pressure_bar), 200)
        self.assertAlmostEqual(atm.pressure_bar[0], 1e-15, delta=1e-27)
        self.assertAlmostEqual(atm.pressure_bar[-1], 10)
    def test_comment_config_keeps_original_keys(self):
        mapping = {'P': 1.2, 'Teff': 8000, 'observatory': 'paranal', 'new_future_key': 'test'}
        cfg = Config.from_text(format_science_config(mapping))
        self.assertTrue(set(mapping).issubset(set(cfg)))
        self.assertEqual(cfg['chemistry'], 'fastchem')
        self.assertNotIn('mock_n_exposures', cfg)
        self.assertIn('# --- Host star', cfg.to_text())
        self.assertIn('# --- Observatory', cfg.to_text())
    def test_espresso_grid(self):
        grid = np.load(instrument_resource('ESPRESSO'))
        blaze = np.load(instrument_resource('ESPRESSO', 'blaze'))
        self.assertEqual(grid.shape, blaze.shape)
        self.assertEqual(len(SpectralOrders.load(instrument_resource('ESPRESSO')).ranges_nm), grid.shape[0])

if __name__ == '__main__':
    unittest.main()
