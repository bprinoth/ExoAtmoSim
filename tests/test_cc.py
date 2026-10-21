"""Fast orbital and CLI smoke tests; do not require POSEIDON installed."""
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from exoatmosim.cc import true_anomaly


class CrossCorrelationTests(unittest.TestCase):
    def test_circular_orbit_anomaly_at_conjunction(self):
        omega = np.deg2rad(90.)
        phi = true_anomaly(np.array([2450000., 2450000.25]),
                           2450000., 1., 0., omega)
        np.testing.assert_allclose(phi, [0., np.pi/2], atol=1e-9)

    def test_eccentric_true_anomaly_not_linear_in_time(self):
        jd = 2450000. + np.array([0., .1, .2])
        phi = true_anomaly(jd, 2450000., 1., .5, np.pi/2)
        self.assertTrue(np.isfinite(phi).all())
        self.assertGreater(abs(phi[1] - phi[0]), abs(phi[2] - phi[1]))

    def test_rejects_invalid_eccentricity(self):
        with self.assertRaises(ValueError):
            true_anomaly(np.array([2450000.]), 2450000., 1., 1., 0.)

    def test_expected_kp_from_keplerian_orbit(self):
        import tempfile
        from exoatmosim.cc import expected_planet_kp
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'config'
            config.write_text(
                'Ms 1.0\nMp 0.001\nP 365.25\ninclination 90\necc 0.0\n')
            kp = expected_planet_kp(config)
            self.assertAlmostEqual(kp, 29.78, delta=0.2)
            config.write_text(
                'Ms 1.0\nMp 0.001\nP 365.25\ninclination 90\necc 0.6\n')
            self.assertAlmostEqual(expected_planet_kp(config) / kp, 1.25, places=8)

    def test_read_observation_berv_uses_fifth_column(self):
        import tempfile
        from exoatmosim.cc import read_observation_berv
        with tempfile.TemporaryDirectory() as d:
            obs = Path(d) / 'obs_times'
            obs.write_text(
                '58926.20 2020-03-18T04:51:17 150.0 1.56 21.69844816\n'
                '58926.21 2020-03-18T05:01:17 150.0 1.49 21.68244106\n')
            np.testing.assert_allclose(
                read_observation_berv(d), [21.69844816, 21.68244106])

    def test_template_wavelengths_convert_to_vacuum(self):
        import tempfile
        from astropy.io import fits
        from exoatmosim.cc import _template
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'O2.fits'
            air_nm = np.linspace(1500., 1501., 100)
            fits.HDUList([fits.PrimaryHDU(np.zeros(1)),
                          fits.ImageHDU(air_nm),
                          fits.ImageHDU(np.linspace(.01, .02, 100))]).writeto(path)
            vacuum_um, _ = _template(path)
            self.assertTrue(np.all(vacuum_um > air_nm / 1000.))
            self.assertTrue(np.all(np.diff(vacuum_um) > 0))

    def test_cli_routes_stage(self):
        from exoatmosim.cli import main
        with patch('exoatmosim.cc.cross_correlate_run', return_value=Path('cc')) as fn:
            self.assertEqual(main(['cc','data/WASP-189_b/config','--seed','0',
                                   '--species','O2','O3','H2O','--sysrem','5',
                                   '--cc-only']), 0)
            self.assertEqual(fn.call_args.kwargs['seed'], 0)
            self.assertTrue(fn.call_args.kwargs['cc_only'])
            self.assertEqual(fn.call_args.kwargs['species'], ('O2','O3','H2O'))


if __name__ == '__main__':
    unittest.main()
