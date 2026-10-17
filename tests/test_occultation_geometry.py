"""Independent reference cases for projected circle-disc overlap."""
import unittest
import numpy as np

from exoatmosim.observing.geometry import occulted_fraction

class TransitGeometryTests(unittest.TestCase):
    def test_central_transit(self):
        self.assertAlmostEqual(occulted_fraction(0.0, 0.1, 1.0), 0.01)

    def test_behind_star_is_not_negative_flux(self):
        self.assertEqual(occulted_fraction(0.0, 0.1, -1.0), 0.0)
        self.assertEqual(occulted_fraction(1.0, 0.1, -1.0), 0.0)

    def test_no_overlap_and_external_tangency(self):
        self.assertEqual(occulted_fraction(1.1, 0.1, 1.0), 0.0)
        self.assertEqual(occulted_fraction(2.0, 0.1, 1.0), 0.0)

    def test_full_overlap_to_internal_tangency(self):
        self.assertAlmostEqual(occulted_fraction(0.9, 0.1, 1.0), 0.01)

    def test_partial_overlap_bounds(self):
        x = occulted_fraction(1.0, 0.1, 1.0)
        self.assertTrue(0.0 < x < 0.01)
        self.assertAlmostEqual(x, occulted_fraction(1.0, 0.1, 1.0))

    def test_invalid_geometry(self):
        for d, rp, zp in [(np.nan, .1, 1), (-.01, .1, 1), (.1, 0, 1), (.1, .1, np.nan)]:
            with self.assertRaises(ValueError):
                occulted_fraction(d, rp, zp)

    def test_continuity_near_contacts(self):
        for contact in (0.9, 1.1):
            v = [occulted_fraction(x, .1, 1) for x in (contact-1e-6, contact, contact+1e-6)]
            self.assertLess(max(v)-min(v), 1e-5)

if __name__ == '__main__':
    unittest.main()
