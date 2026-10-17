import csv
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from exoatmosim.observing.viola_detector import (
    read_viola_detector_csv, pixel_variance, inject_viola_planet,
    make_viola_observation)


class Cfg:
    d = {'instrument':'VIOLA_K', 'viola_wavelength_medium':'air',
         'Kmag':'8.0','Hmag':'8.5','Ymag':'9.0','Teff':'5000','exposure_time':'150',
         'viola_pwv_mm':'2','viola_pwv_sigma_mm':'0.5',
         'viola_pwv_seed':'17','resolution':'100000', 'tellurics':'False'}
    def get(self, k, default=None): return self.d.get(k, default)
    def get_float(self,k,default=None): return float(self.get(k,default))
    def get_int(self,k,default=None): return int(self.get(k,default))
    def get_bool(self,k,default=None): return str(self.get(k,default)).lower() == 'true'


def write_csv(folder, counts=25):
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/'signal_noise.csv').open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['wavelength_um','stellar_e_per_pixel','stellar_e_per_resel',
                         'snr_per_resel','total_noise_e_per_resel'])
        for wave in (1.925,1.930,1.935):
            writer.writerow((wave,counts,counts*3,5.,15.))


class DetectorTests(unittest.TestCase):
    def test_electrons_per_pixel_and_no_double_rv(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_csv(Path(tmp))
            rec=read_viola_detector_csv(Path(tmp)/'signal_noise.csv')
        np.testing.assert_array_equal(rec.stellar_e_pixel,[25]*3)
        # 15^2 resel variance; 75 stellar shot variance; N_eff=3
        np.testing.assert_allclose(pixel_variance(rec),[75]*3)
        grid, star, expected, sigma, observed=inject_viola_planet(
            [rec,rec], [1910,1920,1925,1930,1935,1940,1950],
            [1,1,1,0.8,1,1,1], [0.,0.], [1.,0.], add_noise=False)
        np.testing.assert_allclose(star,25.)
        self.assertAlmostEqual(expected[0,1],20.)
        self.assertAlmostEqual(expected[1,1],25.)
        np.testing.assert_allclose(observed,expected)
        self.assertGreater(sigma[1,1],sigma[0,1])

    def test_zero_snr_nan_noise_retained_with_quality_mask(self):
        with tempfile.TemporaryDirectory() as td:
            file = Path(td) / 'signal_noise.csv'
            with file.open('w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['wavelength_um', 'stellar_e_per_resel',
                            'stellar_e_per_pixel', 'snr_per_resel',
                            'total_noise_e_per_resel'])
                w.writerow([1.900, .08026744, .02006686, .00188747, 42.52636])
                w.writerow([1.901, 0, 0, 0, 'nan'])
                w.writerow([1.902, 0, 0, 0, 'nan'])
            etc = read_viola_detector_csv(file)
            np.testing.assert_array_equal(etc.valid_noise, [True, False, False])
            var = pixel_variance(etc)
            self.assertTrue(np.isfinite(var[0]))
            self.assertTrue(np.isnan(var[1:]).all())
            output = inject_viola_planet([etc], [1890, 1900, 1910, 1920],
                                          [1, 1, 1, 1], [0], [0], seed=1)
            grid, star, model, sigma, noisy = output
            self.assertEqual(len(grid), 3)
            self.assertEqual(star[0, 1], 0)
            self.assertTrue(np.isnan(sigma[0, 1]))
            self.assertTrue(np.isnan(noisy[0, 1]))
            # Negative physical counts are not masked or silently clipped.
            content = file.read_text().replace('0,0,0,nan', '-1,0,0,nan', 1)
            file.write_text(content)
            with self.assertRaisesRegex(ValueError, 'invalid physical detector'):
                read_viola_detector_csv(file)

    def test_end_to_end_band_and_temporary_cleanup(self):
        try:
            from astropy.io import fits
        except ImportError:
            self.skipTest('Astropy is needed for FITS integration')
        from copy import deepcopy
        cfg=Cfg()
        schedule=SimpleNamespace(altitude_deg=np.array([52.,50.]),
                     midpoint_jd_utc=np.array([2460000.,2460000.003]))
        with tempfile.TemporaryDirectory() as tmp:
            seen=[]
            def fake_viola(**kw):
                seen.append(kw)
                folder=Path(kw['output_folder'])
                write_csv(folder,counts=25 + len(seen))
                os.chdir(folder)  # Deliberately emulate VIOLA's side effect.
                return folder
            original=Path.cwd()
            for band in ('K','H','Y'):
                cfg.d=dict(Cfg.d, instrument='VIOLA_'+band)
                out=Path(tmp)/band
                make_viola_observation(schedule,[6.,6.2],[-5.,5.],[0.,1.],cfg,
                    [1910,1920,1925,1930,1935,1940,1950], [1.,1.,1.,.9,1.,1.,1.],
                    out,run_viola=fake_viola)
                self.assertEqual(Path.cwd(),original)
                self.assertFalse(any(p.is_dir() for p in out.iterdir()))
                self.assertEqual(fits.getdata(out/'stellar_flux_0.fits').shape,(2,3))
                self.assertEqual(fits.getdata(out/'noise_0.fits').shape,(2,3))
                self.assertTrue((out/f'VIOLA_{band}_orders_nm.npy').is_file())
            self.assertEqual([r['band'] for r in seen],['K','K','H','H','Y','Y'])
            self.assertEqual([r['radial_velocity_km_s'] for r in seen],[6.,6.2]*3)
            self.assertNotEqual(seen[0]['pwv_mm'],seen[1]['pwv_mm'])

    def test_recover_existing_viola_without_etc_and_without_deleting_sources(self):
        try:
            from astropy.io import fits
        except ImportError:
            self.skipTest('Astropy needed for VIOLA FITS recovery')
        cfg = Cfg()
        schedule = SimpleNamespace(altitude_deg=np.array([52., 50.]),
                      midpoint_jd_utc=np.array([2460000., 2460000.003]))
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / 'old' / 'viola_stellar'
            for i in range(2):
                write_csv(root / f'exposure_{i:04d}' / 'K', counts=25+i)
            with (root / 'viola_conditions.csv').open('w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(('exposure_id', 'midpoint_minutes', 'altitude_deg',
                                 'stellar_observer_rv_kms', 'pwv_mm'))
                writer.writerow((0, 0, 52, 6, 2.3))
                writer.writerow((1, 4.32, 50, 6.2, 2.4))
            def must_not_call_etc(**kw):
                raise AssertionError('VIOLA ETC must not be called during recovery')
            out = Path(temp) / 'simulations' / 'VIOLA_K' / 'seed_0042'
            make_viola_observation(schedule, [6., 6.2], [-5., 5.], [0., 1.],
                cfg, [1910,1920,1925,1930,1935,1940,1950],
                [1,1,1,.9,1,1,1], out, resume_from=root, run_viola=must_not_call_etc)
            self.assertEqual(fits.getdata(out/'stellar_flux_0.fits').shape, (2,3))
            self.assertTrue((root/'exposure_0001/K/signal_noise.csv').is_file())
            self.assertIn('2.3', (out/'viola_conditions.csv').read_text())
            self.assertTrue((out/'VIOLA_K_orders_nm.npy').is_file())
            with self.assertRaisesRegex(ValueError, 'mismatch'):
                make_viola_observation(schedule, [6., 7.], [-5., 5.], [0.,1.],
                    cfg, [1910,1920,1925,1930,1935,1940,1950],
                    [1,1,1,.9,1,1,1], Path(temp)/'bad',
                    resume_from=root, run_viola=must_not_call_etc)

    def test_rejects_second_telluric_model(self):
        cfg=Cfg(); cfg.d=dict(Cfg.d,tellurics='True')
        with self.assertRaisesRegex(ValueError,'already contains tellurics'):
            make_viola_observation(SimpleNamespace(),[],[],[],cfg,[],[],'.')

if __name__ == '__main__':unittest.main()
