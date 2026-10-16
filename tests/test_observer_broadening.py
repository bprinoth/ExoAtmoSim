"""Physical regression tests for the independent observer-frame and LSF path."""
import unittest
import numpy as np
from exoatmosim.observing.spectral_processing import (
    observer_frame_velocities, broaden_spectrum, C_KMS)
from exoatmosim.observing.order_mock import simulate_orders


class TestObserverFrame(unittest.TestCase):
    def test_observer_velocity_components(self):
        star, planet = observer_frame_velocities(
            np.array([1., -2.]), np.array([-120., 130.]), 13.,
            np.array([7., 9.]))
        np.testing.assert_allclose(star, [7., 2.])
        np.testing.assert_allclose(planet, [-114., 134.])

    def test_missing_berv_or_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            observer_frame_velocities([1., 2.], [1., 2.], 0., [1.])
        with self.assertRaises(ValueError):
            observer_frame_velocities([1.], [1.], np.nan, [1.])

    def test_spectral_shift_recovers_expected_velocity(self):
        wave = np.linspace(1928., 1932., 12001)
        line = 1. - 0.5*np.exp(-0.5*((wave-1930.)/0.003)**2)
        orders = wave[4000:8000][None, :]
        observer_vel = 6.
        result = simulate_orders(orders, wave, line, wave, np.ones_like(wave),
                                 np.array([observer_vel]), np.array([0.]),
                                 np.array([0.]), add_noise=False)
        sampled_wavelength, flux, _ = result[0]
        min_wavelength = sampled_wavelength[0, np.argmin(flux[0])]
        np.testing.assert_allclose(min_wavelength, 1930.*(1.+observer_vel/C_KMS), atol=0.0004)

    def test_instrument_broadening_constant_preserved(self):
        wave = np.linspace(1900., 2000., 2000)
        _, flux = broaden_spectrum(wave, np.ones_like(wave)*8.5, resolution=100000.)
        np.testing.assert_allclose(flux, 8.5, rtol=1e-12)

    def test_instrument_broadening_reduces_line_depth(self):
        wave = np.linspace(1929.9, 1930.1, 4001)
        narrow = 1.-0.8*np.exp(-0.5*((wave-1930.)/0.0005)**2)
        w, out = broaden_spectrum(wave, narrow, resolution=100000.)
        self.assertGreater(np.min(out), np.min(narrow))
        self.assertLess(np.min(out), 1.)
        self.assertGreater(w.size, 100)

    def test_rotation_changes_line_shape(self):
        try:
            import PyAstronomy
        except ImportError:
            self.skipTest('PyAstronomy unavailable locally')
        wave = np.linspace(1929.5, 1930.5, 1501)
        line = 1.-0.5*np.exp(-0.5*((wave-1930.)/0.005)**2)
        w, out = broaden_spectrum(wave, line, vsini_kms=8.)
        self.assertGreater(np.min(out), np.min(line))
        self.assertAlmostEqual(w[0], wave[0]); self.assertAlmostEqual(w[-1], wave[-1])


if __name__ == '__main__':
    unittest.main()
