"""Regression tests for mixed Astropy/numeric units in RVTraceEstimator."""
import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples'))
from run_independent_crires import velocity_values_kms

class TestRunnerRVUnits(unittest.TestCase):
    def test_numeric_kms(self):
        np.testing.assert_allclose(velocity_values_kms([1.0, -2.0]), [1, -2])

    def test_dimensionless_legacy_planet_rv(self):
        try:
            import astropy.units as u
        except ImportError:
            self.skipTest('Astropy optional')
        legacy = np.array([14., -25.]) * u.dimensionless_unscaled
        np.testing.assert_allclose(velocity_values_kms(legacy), [14., -25.])

    def test_physical_speed_conversion(self):
        try:
            import astropy.units as u
        except ImportError:
            self.skipTest('Astropy optional')
        np.testing.assert_allclose(velocity_values_kms(np.array([1200., -2300.]) * u.m / u.s), [1.2, -2.3])

    def test_incompatible_units_error(self):
        try:
            import astropy.units as u
        except ImportError:
            self.skipTest('Astropy optional')
        with self.assertRaises(u.UnitConversionError):
            velocity_values_kms(np.array([3.]) * u.m)
