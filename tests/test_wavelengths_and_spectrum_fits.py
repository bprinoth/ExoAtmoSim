"""Internal air/vacuum and planetary FITS tests (no external previous mock framework)."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from exoatmosim.wavelengths import airtovac, vactoair
from exoatmosim.spectrum_fits import write_spectrum_fits
try:
    from astropy.io import fits
except ImportError:
    fits = None

class WavelengthTests(unittest.TestCase):
    def test_conversion_and_round_trip(self):
        vacuum = np.array([400., 500., 1000., 2200.])
        air = vactoair(vacuum)
        self.assertTrue(np.all(air < vacuum))
        np.testing.assert_allclose(airtovac(air), vacuum, rtol=3e-7)

    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            vactoair([np.nan])

@unittest.skipIf(fits is None, 'astropy required')
class FITSTests(unittest.TestCase):
    def test_export_schema_and_overwrite_protection(self):
        spec = SimpleNamespace(wavelength_um=np.array([0.4,0.5,0.6]),
                               transmission=np.array([0.99,0.98,0.97]))
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / 'planet.fits'
            write_spectrum_fits(spec, dest)
            with self.assertRaises(FileExistsError):
                write_spectrum_fits(spec, dest)
            with fits.open(dest) as h:
                self.assertEqual(len(h), 3)
                np.testing.assert_allclose(h[1].data, vactoair(spec.wavelength_um*1000))
                np.testing.assert_allclose(h[2].data, spec.transmission)
