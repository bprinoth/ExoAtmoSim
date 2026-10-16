import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from exoatmosim.observing.order_mock import simulate_orders
from exoatmosim.observing.viola import viola_observer_spectra, read_viola_csv

class MockCfg:
    def __init__(self):
        self.d = {'viola_band':'K','viola_wavelength_medium':'air','Kmag':'7.0',
                  'Teff':'5000','texp':'150','viola_pwv_mm':'2.5','viola_seeing_arcsec':'1.0'}
    def get(self, k, default=None): return self.d.get(k, default)
    def get_float(self,k,default=None): return float(self.d.get(k, default))
    def get_int(self,k,default=None): return int(self.d.get(k, default))

class ViolaObserverFrameTests(unittest.TestCase):
    def test_stellar_frame_never_shifted_twice(self):
        wave=np.linspace(1928,1932,1001)
        template=1-0.4*np.exp(-0.5*((wave-1930)/0.01)**2)
        # With a 100-km/s velocity, a *second* RV application would be obvious.
        rows=np.vstack([template,template*0.7])
        kwargs=dict(order_wavelengths_nm=np.array([[1929.5,1930.,1930.5]]),
                    stellar_wave_nm=wave, stellar_flux=rows,
                    planet_wave_nm=wave, planet_transmission=np.ones_like(wave),
                    rv_star_kms=np.array([100.,-100.]), rv_planet_kms=np.zeros(2),
                    transit_weight=np.zeros(2),snr=100,add_noise=False,
                    stellar_observer_frame=True)
        result=simulate_orders(**kwargs)[0][1]
        np.testing.assert_allclose(result[0], np.interp(kwargs['order_wavelengths_nm'][0],wave,template))
        np.testing.assert_allclose(result[1], 0.7*result[0])

    def test_refuses_telluric_double_counting(self):
        wave=np.linspace(1900,2000,50)
        with self.assertRaisesRegex(ValueError,'already contains tellurics'):
            simulate_orders(wave[None,:],wave,np.ones((1,50)),wave,np.ones(50),
                            [0.],[0.],[0.],stellar_observer_frame=True,
                            telluric_wavelength_nm=wave,
                            telluric_transmission=np.ones((1,50)))

    def test_exact_rv_forwarded_per_exposure(self):
        rv_seen=[]
        pwv_seen=[]
        def fake_viola(**kwargs):
            rv_seen.append(kwargs['radial_velocity_km_s'])
            pwv_seen.append(kwargs['pwv_mm'])
            folder=Path(kwargs['output_folder'])
            with (folder/'signal_noise.csv').open('w',newline='') as f:
                w=csv.writer(f)
                w.writerow(['wavelength_um','stellar_e_per_resel','snr_per_resel'])
                w.writerow([1.92,12,10]);w.writerow([1.93,13,11])
            return str(folder)
        with tempfile.TemporaryDirectory() as tmp:
            w,f,s=viola_observer_spectra(SimpleNamespace(altitude_deg=[50.,49.]),
                    [6.,6.2],MockCfg(),tmp,run_viola=fake_viola)
        self.assertEqual(rv_seen,[6.,6.2])
        self.assertEqual(pwv_seen,[2.5,2.5])
        self.assertEqual(f.shape,(2,2))
        np.testing.assert_allclose(w,[1920.,1930.])

    def test_variable_pwv_reproducible_and_separate_per_exposure(self):
        from exoatmosim.observing.viola import viola_observer_spectra
        cfg=MockCfg()
        cfg.d.update(viola_pwv_sigma_mm='0.5', viola_pwv_seed='777')
        pwv_samples=[]
        def fake(**kwargs):
            pwv_samples.append(kwargs['pwv_mm'])
            folder=Path(kwargs['output_folder'])
            with (folder/'signal_noise.csv').open('w',newline='') as fh:
                writer=csv.writer(fh)
                writer.writerow(['wavelength_um','stellar_e_per_resel','snr_per_resel'])
                writer.writerow([1.92,12,10]); writer.writerow([1.93,13,11])
            return str(folder)
        with tempfile.TemporaryDirectory() as d:
            schedule=SimpleNamespace(altitude_deg=[50.,49.,48.],
                                     midpoint_jd_utc=np.array([2460000,2460000.002,2460000.004]))
            for attempt in range(2):
                viola_observer_spectra(schedule,[3.,4.,5.],cfg,Path(d)/str(attempt),run_viola=fake)
            self.assertTrue((Path(d)/'0'/'viola_conditions.csv').exists())
            np.testing.assert_array_equal(pwv_samples[:3],pwv_samples[3:])
            self.assertGreater(np.std(pwv_samples[:3]),0)

    def test_rejects_unknown_medium(self):
        with self.assertRaisesRegex(ValueError, 'not confirmed AIR'):
            read_viola_csv('nonexistent',wavelength_medium='vacuum')

if __name__=='__main__': unittest.main()
