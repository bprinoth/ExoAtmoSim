import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from exoatmosim.config import Config
from exoatmosim.instruments import instrument_wavelength_bounds_um
from exoatmosim.poseidon import PoseidonTransmissionModel

class TestInstrumentCoverage(unittest.TestCase):
    def test_espresso_orders_define_model_range(self):
        cfg = Config.from_mapping({'instrument': 'ESPRESSO', 'instrument_margin_nm': 5})
        a, b = instrument_wavelength_bounds_um(cfg, 'config')
        self.assertTrue(0.36 < a < 0.39)
        self.assertTrue(0.78 < b < 0.81)
        model = PoseidonTransmissionModel.from_config(cfg)
        self.assertAlmostEqual(model.wavelength_min_um, a)
        self.assertAlmostEqual(model.wavelength_max_um, b)
    def test_overrides_must_cover_orders(self):
        cfg = Config.from_mapping({'instrument': 'ESPRESSO', 'wl_min': 0.6})
        with self.assertRaises(ValueError):
            instrument_wavelength_bounds_um(cfg, 'config')
    def test_custom_orders_are_relative_to_config(self):
        with TemporaryDirectory() as temp:
            folder = Path(temp)
            np.save(folder/'orders.npy', [[510., 520.], [530., 540.]])
            cfg = Config.from_mapping({'instrument_order_file': 'orders.npy'})
            a, b = instrument_wavelength_bounds_um(cfg, folder/'config')
            self.assertAlmostEqual(a, .505)
            self.assertAlmostEqual(b, .545)
