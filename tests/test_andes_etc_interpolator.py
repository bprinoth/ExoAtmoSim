"""Fast standalone tests for nine-file ANDES ETC interpolation."""
import tempfile
import unittest
from pathlib import Path
import numpy as np

from exoatmosim.andes_etc import AndesETCGrid


class AndesETCInterpolatorTests(unittest.TestCase):
    def write_grid(self, root):
        header = ('wavelength [nm],snr [per res. elem.],'
                  'object signal [e-/res. elem.],'
                  'total noise [e-/res. elem.],'
                  'background [e-/s/res. elem.],Efficiencies [Atm]')
        for am in (1.0, 1.5, 3.0):
            folder = root / f'airmass_{am:.1f}'
            folder.mkdir()
            for arm in ('UBV', 'RIZ', 'YJH'):
                # ETC-style timestamped names (no renaming required).
                file = folder / f'20261023_ETCv5.6_{arm}.csv'
                with file.open('w') as f:
                    f.write(header + '\n')
                    for wavelength in (500., 501., 502.):
                        src = 100. / am
                        noise = 10. * am
                        f.write(f'{wavelength},{src/noise},{src},{noise},3,0.7\n')

    def test_timestamped_exports_and_variance_interpolation(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_grid(root)
            grid = AndesETCGrid(root)
            result = grid.evaluate('RIZ', [500., 501.5], [1., 1.25, 3.])
            self.assertEqual(result.snr_resel.shape, (3, 2))
            self.assertAlmostEqual(result.snr_resel[0, 0], 10)
            self.assertAlmostEqual(result.snr_resel[-1, 0], (100./3.)/30.)
            self.assertAlmostEqual(result.source_e_resel[1, 0], 100*(1+1/1.5)/2)
            self.assertAlmostEqual(result.noise_e_resel[1, 0]**2,
                                   (10**2 + 15**2)/2)
            np.testing.assert_allclose(result.snr_resel[:, 0], result.snr_resel[:, 1])

    def test_rounded_order_boundary_duplicate(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_grid(root)
            for am in (1.0, 1.5, 3.0):
                path = root / f'airmass_{am:.1f}' / '20261023_ETCv5.6_UBV.csv'
                lines = path.read_text().splitlines()
                lines.insert(2, lines[1])  # same rounded wavelength at boundary
                path.write_text('\n'.join(lines) + '\n')
            result = AndesETCGrid(root).evaluate('UBV', [500., 501.], [1.25])
            self.assertEqual(result.snr_resel.shape, (1, 2))
            self.assertTrue(np.isfinite(result.snr_resel).all())

    def test_per_resel_to_pixel_noise(self):
        from exoatmosim.andes_etc import spectrum_snr
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_grid(root)
            grid = AndesETCGrid(root)
            # Wavelength interval of 1/300000 in log(lambda), so N~3 for R=100000.
            wave = 500. * np.exp(np.arange(3) / 300000.)
            result = spectrum_snr(grid, wave, [1.0], resolving_power=100000.)
            np.testing.assert_allclose(result[0], 10. / np.sqrt(3.), rtol=1e-4)

    def test_etc_noise_callable_preserves_telluric_noise_floor(self):
        from exoatmosim.observing.order_mock import simulate_orders
        wave = np.array([500., 500.002, 500.004])
        kwargs = dict(order_wavelengths_nm=[wave],
                      stellar_wave_nm=np.array([499., 500., 500.002, 500.004, 501.]),
                      stellar_flux=np.ones(5),
                      planet_wave_nm=np.array([499., 500., 500.002, 500.004, 501.]),
                      planet_transmission=np.ones(5),
                      rv_star_kms=[0.], rv_planet_kms=[0.],
                      transit_weight=[0.], add_noise=False,
                      snr=lambda i, wl: np.full((1, len(wl)), 20.))
        clear = simulate_orders(**kwargs)[0]
        telluric = simulate_orders(
            **kwargs, telluric_wavelength_nm=wave,
            telluric_transmission=np.full((1, len(wave)), .5))[0]
        np.testing.assert_allclose(clear[2], .05)
        np.testing.assert_allclose(telluric[1], .5)
        np.testing.assert_allclose(telluric[2], .05)

    def test_sparse_computational_sampling_keeps_resel_snr(self):
        from exoatmosim.andes_etc import spectrum_snr
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_grid(root)
            grid = AndesETCGrid(root)
            # Coarser than a resolution element: these are not detector pixels.
            wave = np.array([500., 500.1, 500.2])
            result = spectrum_snr(grid, wave, [1.0], resolving_power=100000.)
            np.testing.assert_allclose(result, 10.)

    def test_zero_etc_source_remains_zero_not_floored(self):
        from exoatmosim.andes_etc import spectrum_snr
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_grid(root)
            for am in (1.0, 1.5, 3.0):
                path = root / f'airmass_{am:.1f}' / '20261023_ETCv5.6_UBV.csv'
                rows = path.read_text().splitlines()
                # Force a completely dark ETC wavelength bin (noise remains positive).
                fields = rows[2].split(',')
                fields[1] = '0'
                fields[2] = '0'
                rows[2] = ','.join(fields)
                path.write_text('\\n'.join(rows) + '\\n')
            grid = AndesETCGrid(root)
            values = spectrum_snr(grid, np.array([500., 501., 502.]),
                                  [1.0, 1.25, 3.0], 100000.)
            self.assertTrue(np.all(values[:, 1] == 0))
            self.assertTrue(np.all(values[:, [0, 2]] > 0))

    def test_missing_file_and_extrapolation_raise(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.write_grid(root)
            (root / 'airmass_3.0' / '20261023_ETCv5.6_YJH.csv').unlink()
            with self.assertRaises(FileNotFoundError):
                AndesETCGrid(root)
            self.write_grid_missing_replacement(root)
            grid = AndesETCGrid(root)
            with self.assertRaisesRegex(ValueError, 'airmass outside'):
                grid.evaluate('UBV', [501.], [3.1])
            with self.assertRaisesRegex(ValueError, 'outside ETC range'):
                grid.evaluate('UBV', [503.], [1.])

    @staticmethod
    def write_grid_missing_replacement(root):
        file = root / 'airmass_3.0' / '20261023_ETCv5.6_YJH.csv'
        with file.open('w') as fh:
            fh.write('wavelength [nm],snr [per res. elem.],'
                     'object signal [e-/res. elem.],total noise [e-/res. elem.],'
                     'background [e-/s/res. elem.],Efficiencies [Atm]\n')
            for w in (500., 501., 502.):
                fh.write(f'{w},1,1,1,3,0.7\n')


if __name__ == '__main__':
    unittest.main()
