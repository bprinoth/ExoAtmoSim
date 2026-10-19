import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np

from exoatmosim.observing.telluric_grid import grid_transmission
from exoatmosim.observing.conditions import sample_pwv

try:
    import h5py
except ImportError:
    h5py = None


@unittest.skipIf(h5py is None,'h5py optional')
class TestTelluricHDF5(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)/'telluric_paranal.h5'
        with h5py.File(self.path,'w') as f:
            f.attrs['config_json']=json.dumps({'wavelength_medium':'air','allow_approx_pwv':False})
            f['airmass']=[1.,2.]
            f['pwv_mm']=[1.,3.]
            for i,xx in enumerate(([500.,501.],[502.,503.])):
                g=f.create_group(f'chunk_{i:03d}')
                g['wavelength_air_nm']=xx
                for a in range(2):
                    for p in range(2):
                        # exact linear optical-depth truth
                        tau=0.01*np.asarray(xx)+0.1*a+0.2*p
                        g[f'transmission_a{a:02d}_p{p:02d}']=np.exp(-tau).astype('f4')

    def test_bilinear_optical_depth_and_chunk_gap(self):
        schedule=SimpleNamespace(airmass=np.array([1.5,1.]))
        waves=np.array([500.,500.5,501.,501.5,502.,503.])
        x,out=grid_transmission(schedule,waves,pwv_mm=[2.,1.],grid_file=self.path)
        np.testing.assert_array_equal(x,waves)
        np.testing.assert_allclose(out[0],np.exp(-(0.01*waves+0.15)),rtol=1e-6)
        np.testing.assert_allclose(out[1],np.exp(-0.01*waves),rtol=1e-6)

    def test_missing_cell_rejected(self):
        with h5py.File(self.path,'a') as f:
            del f['chunk_001']['transmission_a00_p00']
        with self.assertRaises((ValueError,KeyError)):
            grid_transmission(SimpleNamespace(airmass=np.array([1.])),
                              [500.,503.],pwv_mm=[1.],grid_file=self.path)

    def test_outside_wavelength_rejected(self):
        with self.assertRaisesRegex(ValueError,'coverage'):
            grid_transmission(SimpleNamespace(airmass=np.array([1.])),
                              [499.,500.],pwv_mm=[1.],grid_file=self.path)

    def test_outside_pwv_rejected(self):
        with self.assertRaisesRegex(ValueError,'outside supported'):
            grid_transmission(SimpleNamespace(airmass=np.array([1.])),
                              [500.,501.],pwv_mm=[5.],grid_file=self.path)

    def test_pwv_reproducible(self):
        s=SimpleNamespace(midpoint_jd_utc=np.array([2460000.,2460000.01,2460000.02]))
        a=sample_pwv(s,seed=10)
        np.testing.assert_array_equal(a,sample_pwv(s,seed=10))
        self.assertTrue(np.all(a>0))

if __name__=='__main__':unittest.main()
