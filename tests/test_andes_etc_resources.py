import unittest
import numpy as np
from exoatmosim.instruments import load_instrument_orders
from exoatmosim.observing.order_mock import simulate_orders

class AndesETCTests(unittest.TestCase):
    def test_combined_and_k(self):
        optical = load_instrument_orders('ANDES')
        infrared = load_instrument_orders('ANDES_K')
        self.assertGreater(len(optical), 40)
        self.assertEqual(len(infrared), 19)
        self.assertEqual(sum(len(w) for w in optical), 31350 + 7808 + 162285)
        self.assertTrue(all(2 <= len(w) <= 4096 for w in optical + infrared))
        self.assertLess(min(w[0] for w in optical), 350)
        self.assertGreater(max(w[-1] for w in optical), 1790)
        self.assertGreater(min(w[0] for w in infrared), 1900)

    def test_jagged_segments(self):
        orders = [np.linspace(1000,1001,4), np.linspace(1002,1003,3)]
        wave=np.linspace(999,1004,1001)
        result=simulate_orders(orders,wave,np.ones_like(wave),wave,np.full_like(wave,.99),
            np.array([0.]),np.array([0.]),np.array([1.]),snr=100,add_noise=False)
        self.assertEqual([v[0].shape for v in result],[(1,4),(1,3)])
