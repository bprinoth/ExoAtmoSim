"""Regression tests for real bundled instrumental arrays and future placeholders."""
import unittest
import numpy as np
from exoatmosim.instruments import available_instruments, instrument_resource
from exoatmosim.observations import SpectralOrders


class TestPackagedInstruments(unittest.TestCase):
    def test_instrument_names_stay_supported(self):
        inventory = available_instruments()
        for name in ('ESPRESSO', 'HARPS', 'HARPSN', 'SPIROU', 'CRIRES_K', 'IGRINS_HK', 'IGRINS2_HK'):
            self.assertIn(name, inventory)
        self.assertTrue(inventory['ESPRESSO']['orders'])
        self.assertTrue(inventory['ESPRESSO']['blaze'])
        self.assertFalse(inventory['HARPS']['orders'])
        self.assertFalse(inventory['HARPSN']['orders'])
        self.assertTrue(inventory['HARPSN']['blaze'])

    def test_valid_order_grids_bundled_and_loadable(self):
        for name in ('ESPRESSO', 'SPIROU', 'CRIRES_K', 'IGRINS_HK', 'IGRINS2_HK'):
            with self.subTest(name=name):
                source = instrument_resource(name)
                raw = np.load(source, allow_pickle=False)
                self.assertGreater(raw.shape[1], 2)
                reduced = SpectralOrders.load(source).ranges_nm
                self.assertEqual(len(reduced), len(raw))
                self.assertTrue(np.all(reduced[:, 1] > reduced[:, 0]))

    def test_espresso_170_strictly_increasing_orders(self):
        wave = np.load(instrument_resource('ESPRESSO', 'orders'), allow_pickle=False)
        blaze = np.load(instrument_resource('ESPRESSO', 'blaze'), allow_pickle=False)
        self.assertEqual(wave.shape, (170, 9111))
        self.assertEqual(blaze.shape, wave.shape)
        self.assertTrue(np.isfinite(wave).all())
        self.assertTrue(np.all(np.diff(wave, axis=1) > 0))

    def test_valid_blazes_bundled(self):
        for name in ('ESPRESSO', 'SPIROU', 'HARPSN'):
            with self.subTest(name=name):
                values = np.load(instrument_resource(name, 'blaze'), allow_pickle=False)
                self.assertEqual(values.ndim, 2)
                # Real blaze calibration arrays can contain NaNs in invalid pixels.
                # Preserve masks rather than inventing values in those pixels.
                self.assertGreater(np.isfinite(values).mean(), 0.9)
                self.assertFalse(np.isinf(values).any())

    def test_missing_orders_fail_explicitly(self):
        for name in ('HARPS', 'HARPS-N'):
            with self.subTest(name=name), self.assertRaisesRegex(FileNotFoundError, 'no valid bundled orders'):
                instrument_resource(name)
