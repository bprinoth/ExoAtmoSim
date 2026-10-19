import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from exoatmosim.observing.viola_native import sample_viola_pwv, select_nearest_skycalc

class ViolaSkySamplingTests(unittest.TestCase):
    def setUp(self):
        self.schedule = SimpleNamespace(midpoint_jd_utc=2450000. + np.arange(12)*5/1440.)

    def test_pwv_reproducible_and_variable(self):
        a = sample_viola_pwv(self.schedule, sigma_mm=0.4, seed=99)
        b = sample_viola_pwv(self.schedule, sigma_mm=0.4, seed=99)
        np.testing.assert_array_equal(a,b)
        self.assertTrue(np.all(a > 0))
        self.assertGreater(np.ptp(a),0)
        self.assertFalse(np.array_equal(a,sample_viola_pwv(self.schedule,sigma_mm=0.4,seed=100)))

    def test_constant_pwv(self):
        np.testing.assert_array_equal(sample_viola_pwv(self.schedule,mean_mm=2.5,sigma_mm=0),np.full(12,2.5))

    def test_nearest_2d_and_clamping(self):
        candidates = [(Path(f'{a}_{w}.fits'),a,w) for a in (30.,60.,90.) for w in (0.5,2.5,10.)]
        result = select_nearest_skycalc(candidates,[28.,61.,89.],[0.1,2.4,40.])
        self.assertEqual([(a,w) for _,a,w in result],[(30.,0.5),(60.,2.5),(90.,10.)])

    def test_missing_grid_cell_is_an_error(self):
        candidates = [(Path('a'),30.,0.5),(Path('b'),60.,2.5)]
        with self.assertRaisesRegex(ValueError,'Missing SkyCalc'):
            select_nearest_skycalc(candidates,[30.],[2.5])

if __name__ == '__main__': unittest.main()
