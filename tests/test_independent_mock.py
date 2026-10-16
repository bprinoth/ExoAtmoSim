import tempfile
import unittest
from pathlib import Path
import numpy as np
from exoatmosim.observing.order_mock import simulate_orders
from exoatmosim.observing.scheduling import observation_window_hours
from exoatmosim.config import Config


class WindowTests(unittest.TestCase):
    def test_T14_plus_two_supersedes_legacy_t(self):
        c = Config.from_mapping({'T14':3.2,'t':8.6141658})
        self.assertAlmostEqual(observation_window_hours(c), 5.2)
    def test_override(self):
        c = Config.from_mapping({'T14':3.2,'t':8.6141658,'observation_window_hours':4.0})
        self.assertEqual(observation_window_hours(c), 4.0)
    def test_legacy_fallback(self):
        self.assertEqual(observation_window_hours(Config.from_mapping({'t':2.5})), 2.5)

class OrderTests(unittest.TestCase):
    def setUp(self):
        self.w = np.linspace(1950, 1970, 1001)
        self.o = np.tile(np.linspace(1954, 1964, 100), (2, 1))
        self.star = 1.0 - 0.1 * np.exp(-0.5*((self.w-1958)/0.08)**2)
        self.pl = 1.0 - 0.04 * np.exp(-0.5*((self.w-1960)/0.08)**2)
    def calc(self, *, noise=False, mask=None):
        return simulate_orders(self.o,self.w,self.star,self.w,self.pl,
                               [0., 0.], [0., 0.], mask if mask is not None else [0.,1.],
                               snr=100, seed=7, add_noise=noise)
    def test_deterministic_planet_injection(self):
        result = self.calc()[0]
        wave, flux, noise = result
        np.testing.assert_allclose(flux[0], np.interp(wave[0],self.w,self.star),rtol=0,atol=1e-12)
        expect = np.interp(wave[1],self.w,self.star)*np.interp(wave[1],self.w,self.pl)
        np.testing.assert_allclose(flux[1],expect,rtol=0,atol=1e-12)
    def test_noise_seed_and_shape(self):
        a=self.calc(noise=True);b=self.calc(noise=True)
        self.assertEqual(len(a),2)
        self.assertEqual(a[0][1].shape,(2,100))
        np.testing.assert_array_equal(a[0][1],b[0][1])
        self.assertTrue(np.any(a[0][1]!=self.calc()[0][1]))
    def test_reject_missing_wavelength_coverage(self):
        with self.assertRaisesRegex(ValueError,'not fully covered'):
            simulate_orders(self.o,self.w[:30],self.star[:30],self.w,self.pl,[0],[0],[1])

try:
    from astropy.io import fits
    HAS_ASTROPY=True
except ImportError:
    HAS_ASTROPY=False

@unittest.skipUnless(HAS_ASTROPY, "requires astropy")
class OutputTests(unittest.TestCase):
    def test_roundtrip(self):
        from exoatmosim.observing.order_mock import write_tayph_orders

        with tempfile.TemporaryDirectory() as temp:
            orders = OrderTests()
            orders.setUp()
            simulated = orders.calc()

            products = write_tayph_orders(simulated, temp)
            self.assertEqual(len(products), 4)

            with fits.open(Path(temp) / "order_1.fits") as hdul:
                self.assertEqual(hdul[0].data.shape, (2, 100))

            with self.assertRaises(FileExistsError):
                write_tayph_orders(simulated, temp)