import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock
import numpy as np
from exoatmosim.config import Config
from exoatmosim.orbits import RVTraceEstimator
from exoatmosim.manual_chemistry import manual_profile
from exoatmosim.tayph_mock import run_tayph_mock


class TestV0811(unittest.TestCase):
    def test_circular_nan_omega_is_physical(self):
        config = Config.from_mapping({'P':2.7, 'Tc':2450000., 'ecc':0., 'omega':'nan',
                                      'K':.18, 'Ms':2., 'Mp':2.})
        rv = RVTraceEstimator.from_config(config).evaluate(2450000.+np.linspace(-.03,.03,7))
        self.assertTrue(np.isfinite(rv['true_anomaly_rad']).all())
        self.assertTrue(np.isfinite(rv['rv_planet_kms']).all())
        self.assertAlmostEqual(rv['rv_star_kms'][3],0,places=8)

    def test_eccentric_requires_omega(self):
        with self.assertRaisesRegex(ValueError,'finite omega'):
            RVTraceEstimator(2., 2450000., .3, float('nan'), .1, 1000.)

    def test_manual_is_explicit(self):
        atm=type('FakeAtmosphere',(),{'pressure_bar':np.array([1e-8,1e-5,1.])})()
        cfg=Config.from_mapping({'chemistry':'manual','bulk_species':'H2 He',
                                  'manual_abundances':'H2=0.75 He=0.24 H-=1e-8 e-=1e-6'})
        profile=manual_profile(cfg,atm)
        self.assertEqual(profile.mixing_ratios['H-'][0],1e-8)
        self.assertAlmostEqual(profile.mixing_ratios['H2'][2],.75)
        with self.assertRaises(ValueError):
            manual_profile(Config.from_mapping({'bulk_species':'H2 He','manual_abundances':'H2=0.9 He=0.5'}),atm)

    def make_planet(self,root,ecc,omega='nan'):
        d=root/'data'/'WASP-189_b';d.mkdir(parents=True)
        cfg=d/'config';cfg.write_text(f'instrument\tESPRESSO\nP\t2.724\nTc\t2458926.5\necc\t{ecc}\nomega\t{omega}\nK\t0.182\nMs\t2.03\nMp\t1.99\n')
        (d/'WASP-189_b.fits').write_bytes(b'present')
        return cfg

    def test_circular_dispatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cfg=self.make_planet(root,0)
            prepared,kwargs=run_tayph_mock(cfg,project_root=root,dry_run=True)
            self.assertTrue(prepared.runfile.exists())
            self.assertTrue(kwargs['eccentric_orbit'])
            self.assertNotIn('bjds',kwargs)
            self.assertNotIn('phase_range',kwargs)

    def test_eccentric_dispatch_true_anomaly(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cfg=self.make_planet(root,.4,omega='33')
            bjds=2458926.5+np.linspace(-.03,.03,5)
            prepared,kwargs=run_tayph_mock(cfg,project_root=root,bjds=bjds,dry_run=True)
            self.assertTrue(kwargs['eccentric_orbit'])
            np.testing.assert_array_equal(kwargs['bjds'],bjds)
            self.assertNotIn('phase_range', kwargs)  # Computed by real RVTraceEstimator at execution time
            self.assertEqual(kwargs['mode'],'ESPRESSO')

    def test_eccentric_missing_bjds_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cfg=self.make_planet(root,.4,omega='33')
            prepared, kwargs = run_tayph_mock(cfg,project_root=root,dry_run=True)
            self.assertTrue(kwargs['eccentric_orbit'])
            self.assertNotIn('bjds', kwargs)  # Original estimator may read obs_times


if __name__=='__main__':
    unittest.main()
