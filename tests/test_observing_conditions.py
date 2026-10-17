"""Regression checks for self-describing conditions and per-exposure TelFit calls."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from exoatmosim.observing.conditions import (observing_conditions, write_observing_conditions,
                                             read_observing_conditions)
from exoatmosim.observing.obs_times import write_observation_times
from exoatmosim.observing.tellurics import telfit_transmission, vacuum_to_air_nm


def schedule(n=8):
    t = 2460000. + np.arange(n) * 210. / 86400.
    return SimpleNamespace(midpoint_jd_utc=t, altitude_deg=np.linspace(75, 60, n),
                           airmass=np.linspace(1.035, 1.15, n),
                           berv_kms=np.linspace(-8, -7, n), duration_seconds=150.)


try:
    import astropy.table  # noqa: F401
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False


@unittest.skipUnless(HAS_ASTROPY, 'requires Astropy for ECSV tables')
class ConditionsTests(unittest.TestCase):
    def test_reproducible_correlated_humidity_and_ecsv(self):
        s = schedule(103)
        kw = dict(humidity_percent=40, humidity_sigma_percent=5,
                  humidity_timescale_minutes=30, seed=7)
        a = observing_conditions(s, **kw)
        b = observing_conditions(s, **kw)
        np.testing.assert_array_equal(a['humidity_percent'], b['humidity_percent'])
        self.assertGreater(np.std(a['humidity_percent']), 0)
        self.assertTrue(np.all((a['humidity_percent'] >= 0) & (a['humidity_percent'] <= 100)))
        self.assertEqual(a.meta['atmosphere_seed'], 7)
        with tempfile.TemporaryDirectory() as d:
            write_observing_conditions(Path(d) / 'observing_conditions.ecsv', a)
            loaded = read_observing_conditions(Path(d) / 'observing_conditions.ecsv', s)
            np.testing.assert_array_equal(loaded['humidity_percent'], a['humidity_percent'])
            self.assertTrue(loaded['humidity_percent'].unit.is_equivalent('percent'))
            write_observation_times(Path(d) / 'obs_times', s)
            lines = (Path(d) / 'obs_times').read_text().splitlines()
            self.assertEqual(len(lines), 103)
            self.assertTrue(all(len(line.split()) == 5 for line in lines))

    def test_telfit_receives_distinct_states_and_clips_small_undershoot(self):
        s = schedule(3)
        conditions = observing_conditions(s, humidity_percent=50, humidity_sigma_percent=10, seed=4)
        calls = []
        class Modeler:
            def MakeModel(self, **kwargs):
                calls.append(kwargs)
                wvac = np.linspace(1e7/kwargs['highfreq'], 1e7/kwargs['lowfreq'], 200)
                transmission = np.ones(200) * .9
                # Put negative samples INSIDE the requested grid, not only
                # in TelFit's padded boundary (which interpolation discards).
                transmission[70:130] = -1.8e-7
                return SimpleNamespace(x=vacuum_to_air_nm(wvac), y=transmission)
        with self.assertWarns(RuntimeWarning):
            w, spectra = telfit_transmission(s, np.linspace(1925, 1935, 100),
                                             conditions=conditions, modeler_factory=Modeler)
        self.assertEqual(len(calls), 3)
        np.testing.assert_allclose([c['humidity'] for c in calls], conditions['humidity_percent'])
        np.testing.assert_allclose([c['angle'] for c in calls], conditions['zenith_angle_deg'])
        np.testing.assert_allclose([c['pressure'] for c in calls], conditions['pressure_hpa'])
        self.assertEqual(spectra[0, len(w)//2], 0.)
        self.assertAlmostEqual(w[0], 1925.0)
        self.assertAlmostEqual(w[-1], 1935.0)

    def test_bad_alignment_rejected(self):
        s = schedule(3)
        table = observing_conditions(s)
        with self.assertRaisesRegex(ValueError, 'one row'):
            telfit_transmission(s, [1925, 1935], conditions=table[:2], modeler_factory=lambda: None)
        table['midpoint_jd_utc'][0] += 1
        with tempfile.TemporaryDirectory() as d:
            path = write_observing_conditions(Path(d) / 'conditions.ecsv', table)
            with self.assertRaisesRegex(ValueError, 'timestamps'):
                read_observing_conditions(path, s)

if __name__ == '__main__':
    unittest.main()
