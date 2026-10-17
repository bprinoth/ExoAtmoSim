"""TelFit exposure spectra can have different input wavelength samplings."""
import unittest
from types import SimpleNamespace
import numpy as np
from exoatmosim.observing.tellurics import telfit_transmission


class GridsTests(unittest.TestCase):
    def test_resample_jittered_exposure_grids(self):
        class Fake:
            counter = 0
            def MakeModel(self, **kwargs):
                i = Fake.counter
                Fake.counter += 1
                x = np.linspace(999.5 - .0005*i, 1002.5 + .0003*i, 90 + i)
                return SimpleNamespace(x=x, y=np.full(x.shape, .5 + i*.1))
        schedule = SimpleNamespace(altitude_deg=np.array([60., 65., 70.]))
        grid = np.linspace(1000.,1002.,21)
        x,y=telfit_transmission(schedule,grid, modeler_factory=Fake)
        self.assertEqual(len(x),21)
        self.assertEqual(y.shape,(3,21))
        self.assertTrue(np.allclose(y[:,0],[.5,.6,.7]))

    def test_no_extrapolation(self):
        class Fake:
            def MakeModel(self, **kwargs):
                return SimpleNamespace(x=np.linspace(1000.2,1001.8,40), y=np.ones(40))
        schedule = SimpleNamespace(altitude_deg=np.array([60.]))
        with self.assertRaisesRegex(ValueError,'cover'):
            telfit_transmission(schedule,np.linspace(1000,1002,20),modeler_factory=Fake)
