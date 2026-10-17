import unittest
import numpy as np
from exoatmosim.observing import exposure_midpoints
from exoatmosim.observing.scheduling import _exposure_seconds
from exoatmosim.config import Config


class SchedulerCadenceTests(unittest.TestCase):
    def test_requested_crires_cadence(self):
        starts, mids, ends = exposure_midpoints(2460000., 8.6141658, 150, 60)
        self.assertEqual(len(mids), 147)
        np.testing.assert_allclose((ends-starts)*86400, 150, atol=1e-4)
        np.testing.assert_allclose(np.diff(mids)*86400, 210, atol=1e-4)
        self.assertLessEqual(ends[-1], 2460000 + 8.6141658/48)

    def test_no_exposure_overrun(self):
        s, m, e = exposure_midpoints(2450000., 1/12, 100, 60)
        self.assertEqual(len(m), 2)
        self.assertLess(e[0], 2450000. + 1/12/48 + 1e-10)

    def test_invalid_inputs(self):
        for args in [(2450000., 0, 20, 0), (2450000., 2, 0, 0), (2450000., 2, 20, -1), (2450000., 0.001, 20, 0)]:
            with self.assertRaises(ValueError):
                exposure_midpoints(*args)

    def test_exposure_time_required_and_positive(self):
        self.assertEqual(_exposure_seconds(Config.from_mapping({'exposure_time': 150})), 150)
        with self.assertRaises(ValueError):
            _exposure_seconds(Config.from_mapping({}))
        with self.assertRaises(ValueError):
            _exposure_seconds(Config.from_mapping({'exposure_time': -2}))
        with self.assertRaises(ValueError):
            _exposure_seconds(Config.from_mapping({'texp': 150}))


try:
    import astropy  # noqa
except ImportError:
    HAVE_ASTROPY = False
else:
    HAVE_ASTROPY = True


@unittest.skipUnless(HAVE_ASTROPY, 'astropy optional: run on mantis with exoatmosim[orbits]')
class TargetGeometryTests(unittest.TestCase):
    def setUp(self):
        from astropy.utils import iers
        iers.conf.auto_download = False
        self.base = {'RA': '12:00:00', 'DEC': '-20:00:00', 'lat': -24.6, 'long': -70.4, 'elev': 2635,
                     'Tc': 2460000.0, 'observation_window_hours': 2, 'exposure_time': 150, 'readout': 60}

    def test_tobs_is_exact_center(self):
        from exoatmosim.observing import make_observation_schedule
        cfg = Config.from_mapping({**self.base, 'Tobs':2460000.125})
        s = make_observation_schedule(cfg, compute_berv=False)
        self.assertEqual(s.mode, 'Tobs')
        self.assertEqual(s.center_jd_utc, 2460000.125)
        self.assertEqual(s.n_exposures, 34)

    def test_max_altitude_near_anchor(self):
        from exoatmosim.observing import make_observation_schedule
        cfg = Config.from_mapping(self.base)
        s = make_observation_schedule(cfg, compute_berv=False)
        self.assertEqual(s.mode, 'max_altitude')
        self.assertLess(abs(s.center_jd_utc-self.base['Tc']), .5)
        self.assertGreater(s.altitude_deg[len(s.altitude_deg)//2], s.altitude_deg[0])
        self.assertGreater(s.altitude_deg[len(s.altitude_deg)//2], s.altitude_deg[-1])
