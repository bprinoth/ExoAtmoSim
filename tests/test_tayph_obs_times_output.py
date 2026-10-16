import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from exoatmosim.observing.obs_times import read_obs_times_mjd, write_tayph_obs_times

try:
    import astropy  # noqa: F401
    ASTROPY = True
except ImportError:
    ASTROPY = False


@unittest.skipUnless(ASTROPY, 'requires astropy')
class TayphObservationTableTests(unittest.TestCase):
    def schedule(self, **changes):
        data = dict(midpoint_jd_utc=np.array([2458926.34, 2458926.35]),
                    airmass=np.array([1.10, 1.27]),
                    berv_kms=np.array([21.5, 21.6]),
                    duration_seconds=150.0)
        data.update(changes)
        return SimpleNamespace(**data)

    def test_five_columns_and_values(self):
        with tempfile.TemporaryDirectory() as td:
            path = write_tayph_obs_times(Path(td) / 'obs_times', self.schedule())
            rows = [line.split() for line in path.read_text().splitlines()]
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(len(row) == 5 for row in rows))
            self.assertAlmostEqual(float(rows[0][2]), 150)
            self.assertAlmostEqual(float(rows[0][3]), 1.10)
            self.assertAlmostEqual(float(rows[1][3]), 1.27)
            self.assertAlmostEqual(float(rows[0][4]), 21.5)
            self.assertIn('T', rows[0][1])
            np.testing.assert_allclose(read_obs_times_mjd(path),
                                       self.schedule().midpoint_jd_utc - 2400000.5, atol=1e-8)

    def test_missing_berv_does_not_generate_placeholder(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'obs_times'
            with self.assertRaisesRegex(ValueError, 'no BERV'):
                write_tayph_obs_times(path, self.schedule(berv_kms=None))
            self.assertFalse(path.exists())

    def test_below_horizon_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(ValueError, 'airmass'):
                write_tayph_obs_times(Path(td) / 'obs_times',
                                      self.schedule(airmass=np.array([1.1, np.nan])))
