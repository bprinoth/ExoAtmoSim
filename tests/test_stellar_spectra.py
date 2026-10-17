import unittest
from unittest.mock import patch
import sys
from types import ModuleType
import numpy as np
from exoatmosim import Config, Star, PhoenixSpectrumProvider, StellarSpectrum

class TestStellarSpectrum(unittest.TestCase):
    def test_spectrum_validation(self):
        with self.assertRaises(ValueError):
            StellarSpectrum([200, 100], [1, 1])

    def test_provider_config(self):
        c = Config.from_mapping({"Teff": 5100, "logg": 4.5, "air": "True",
                                 "wl_min": 0.5, "wl_max": 0.6})
        p = PhoenixSpectrumProvider.from_config(c)
        self.assertFalse(p.vacuum)
        self.assertEqual((p.wavelength_min_um, p.wavelength_max_um), (0.5, 0.6))

    def test_fake_phoenix_provider(self):
        try:
            import astropy.units
        except ImportError:
            self.skipTest("astropy optional dependency not present")
        from exoatmosim import phoenix as pho
        c = Config.from_mapping({"Teff": 5100, "logg": 4.5, "air": "True",
                                 "wl_min": 0.5, "wl_max": 0.6})
        with patch.object(pho, "get_phoenix_model_spectrum", return_value=np.array([1.,2.,3.,4.])), \
             patch.object(pho, "get_phoenix_wavelengths", return_value=np.array([4000.,5000.,6000.,7000.])):
            result = PhoenixSpectrumProvider.from_config(c).get_spectrum(Star.from_config(c))
        np.testing.assert_allclose(result.wavelength_nm, [500., 600.])
        np.testing.assert_allclose(result.flux_erg_s_cm2_nm, [2e-7, 3e-7])
