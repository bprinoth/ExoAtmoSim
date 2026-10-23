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
