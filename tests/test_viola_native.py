"""Fast, deterministic VIOLA ETC noise and photometric reference tests."""
import unittest
import tempfile
from pathlib import Path
import numpy as np
from exoatmosim.config import Config
from exoatmosim.observing.viola_native import (
    calibrate_phoenix_flux, electron_budget, viola_etc_variance)

class ViolaNativeTests(unittest.TestCase):
    def test_variance_matches_etc_runner(self):
        s=np.array([100., 400.]); sky=np.array([5., 20.])
        th=np.array([2., 8.]); dark=np.array([1., 1.])
        expected=s+2*(sky+th+dark)+20.**2*4
        np.testing.assert_allclose(viola_etc_variance(s,sky,th,dark),expected)
        np.testing.assert_allclose(viola_etc_variance(s,sky,th,dark,nod=False),
                                   s+sky+th+dark+1600)

    def test_eight_meter_collecting_area_and_magnitude(self):
        wave=np.linspace(1500,1700,11)
        shape=np.ones_like(wave)
        f10=calibrate_phoenix_flux(wave,shape,10.,'H')
        f11=calibrate_phoenix_flux(wave,shape,11.,'H')
        np.testing.assert_allclose(f11/f10,10**(-0.4))
        c8=Config.from_mapping({'resolution':300000,'viola_telescope_diameter_m':8})
        c2=Config.from_mapping({'resolution':300000,'viola_telescope_diameter_m':2})
        star8, variance8, terms8=electron_budget(wave,f10,np.zeros_like(wave),c8,exposure_s=150)
        star2, _, _=electron_budget(wave,f10,np.zeros_like(wave),c2,exposure_s=150)
        # The VIOLA reference aperture is fixed at 8 m; outdated config overrides are ignored.
        np.testing.assert_allclose(star8,star2)
        self.assertTrue(np.all(star8 > 0))
        self.assertTrue(np.all(variance8 >= star8))
        self.assertTrue(np.all(terms8['thermal_e'] >= 0))

    def test_broadening_returns_matching_wavelength_axis_for_viola(self):
        # VIOLA's internal geometric grid may be resampled by the LSF
        # implementation. Regression: do not interpolate the returned flux
        # against the ORIGINAL high-resolution grid.
        from exoatmosim.observing.spectral_processing import broaden_spectrum
        high = np.geomspace(1900., 1901., 1500)
        flux = 1. - 0.4 * np.exp(-0.5 * ((high - 1900.5) / 0.002)**2)
        blurred_wave, blurred_flux = broaden_spectrum(high, flux, resolution=300000.)
        self.assertEqual(len(blurred_wave), len(blurred_flux))
        self.assertTrue(np.all(np.diff(blurred_wave) > 0))
        sampled = np.interp(np.linspace(1900.1, 1900.9, 100), blurred_wave, blurred_flux)
        self.assertTrue(np.isfinite(sampled).all())
        self.assertLess(sampled.min(), 1.0)
        # The implementation is free to resample: input and output axes
        # must NOT be assumed equal in length.

    def test_oversampling_does_not_duplicate_photons(self):
        cfg=Config.from_mapping({'resolution':300000})
        center=np.array([1600.,1600.001,1600.002])
        flux=np.ones(3)*1e-12
        star,_,_=electron_budget(center,flux,np.zeros(3),cfg,exposure_s=100)
        self.assertLess(np.ptp(star)/np.mean(star),1e-4)
        # Each sample represents the SAME-sized resolution element, not an
        # independent detector bin. Adjacent samples are correlated physically.


    def test_observe_rejects_nan_vsini_before_creating_products(self):
        from exoatmosim.independent_crires import run
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / 'config'
            config.write_text('instrument VIOLA_K\nresolution 300000\nvsini nan\n')
            with self.assertRaisesRegex(ValueError, 'vsini.*finite'):
                run(config, None)
            self.assertFalse((Path(temp) / 'simulations').exists())

    def test_observe_rejects_nonfinite_resolution(self):
        from exoatmosim.independent_crires import run
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / 'config'
            config.write_text('instrument VIOLA_K\nresolution nan\nvsini 4\n')
            with self.assertRaisesRegex(ValueError, 'resolution.*finite'):
                run(config, None)

    def test_rotation_and_lsf_are_separate_operations(self):
        from exoatmosim.observing.spectral_processing import broaden_spectrum
        wave = np.linspace(2100., 2101., 2001)
        line = 1. - .3 * np.exp(-.5 * ((wave - 2100.5) / .003)**2)
        instrumental_wave, instrumental = broaden_spectrum(
            wave, line, resolution=300000.)
        self.assertEqual(len(instrumental_wave), len(instrumental))
        try:
            from PyAstronomy import pyasl
        except ImportError:
            self.skipTest('PyAstronomy unavailable')
        rotated_wave, rotated = broaden_spectrum(
            wave, line, resolution=None, vsini_kms=5.)
        combined_wave, combined = broaden_spectrum(
            rotated_wave, rotated, resolution=300000.)
        # fastRotBroad uses a finite convolution window: its boundary pixels
        # can be depressed by edge padding. Compare the isolated line core,
        # not minima across the full wavelength interval.
        core_instrument = np.abs(instrumental_wave - 2100.5) < 0.06
        core_combined = np.abs(combined_wave - 2100.5) < 0.06
        self.assertGreater(np.min(combined[core_combined]),
                           np.min(instrumental[core_instrument]))
        self.assertTrue(np.isfinite(combined).all())


    def test_viola_exposure_loop_preserves_stellar_rest_frame_grid(self):
        # Regression: the LSF's returned wavelength axis must never replace
        # the stellar-rest-frame axis that subsequent RV shifts interpolate.
        import ast
        import inspect
        import textwrap
        from exoatmosim.observing.viola_native import make_native_viola_observation
        tree = ast.parse(textwrap.dedent(inspect.getsource(make_native_viola_observation)))
        exposure_loop = next(node for node in ast.walk(tree)
                             if isinstance(node, ast.For)
                             and isinstance(node.target, ast.Name)
                             and node.target.id == 'i')
        assignments = [target.id for node in ast.walk(exposure_loop)
                       if isinstance(node, (ast.Assign, ast.AnnAssign))
                       for target in (node.targets if isinstance(node, ast.Assign) else [node.target])
                       for target in (target.elts if isinstance(target, (ast.Tuple, ast.List)) else [target])
                       if isinstance(target, ast.Name)]
        self.assertNotIn('star_wave', assignments)
        self.assertIn('lsf_star_wave', assignments)

if __name__=='__main__': unittest.main()
