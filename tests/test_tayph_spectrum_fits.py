import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace, ModuleType
from unittest.mock import patch
import sys
import numpy as np
from exoatmosim.tayph_mock import write_tayph_spectrum_fits

try:
    from astropy.io import fits
except ImportError:
    fits = None


@unittest.skipIf(fits is None, 'astropy optional dependency not installed')
class TestTayphSpectrumFITS(unittest.TestCase):
    def test_schema_units_and_overwrite_protection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'WASP-189_b.fits'
            spectrum = SimpleNamespace(wavelength_um=np.array([0.37, 0.38, 0.39]),
                                       transmission=np.array([0.99, 0.98, 0.97]))
            tayph = ModuleType('tayph')
            ops = ModuleType('tayph.operations')
            ops.vactoair = lambda wave: np.asarray(wave) - 0.1  # deterministic test stub
            with patch.dict(sys.modules, {'tayph': tayph, 'tayph.operations': ops}):
                write_tayph_spectrum_fits(spectrum, path)
                with self.assertRaises(FileExistsError):
                    write_tayph_spectrum_fits(spectrum, path)
                write_tayph_spectrum_fits(spectrum, path, overwrite=True)
            with fits.open(path) as h:
                self.assertEqual(len(h), 3)
                np.testing.assert_array_equal(h[0].data, [0.])
                np.testing.assert_allclose(h[1].data, [369.9, 379.9, 389.9])
                np.testing.assert_allclose(h[2].data, spectrum.transmission)

    def test_bad_grid_is_rejected(self):
        spectrum = SimpleNamespace(wavelength_um=np.array([0.4, 0.3]),
                                   transmission=np.array([1., 1.]))
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                write_tayph_spectrum_fits(spectrum, Path(tmp) / 'bad.fits')
