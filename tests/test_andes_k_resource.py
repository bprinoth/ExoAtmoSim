"""Experimental ANDES K segments preserve the supplied provisional wavelength samples."""
import unittest
import numpy as np
from exoatmosim.instruments import load_instrument_orders

class AndesKResourceTests(unittest.TestCase):
    def test_provisional_chunks_and_valid_wavelengths(self):
        parts = load_instrument_orders('ANDES_K')
        self.assertEqual(len(parts), 19)
        self.assertEqual(sum(map(len, parts)), 76328)
        self.assertTrue(all(np.isfinite(w).all() and np.all(np.diff(w)>0) for w in parts))
