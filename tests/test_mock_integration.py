import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from exoatmosim.config import Config
from exoatmosim.instruments import resolve_instrument_file
from exoatmosim.orbits import RVTraceEstimator
from exoatmosim.observations import ObservationSettings, SpectralOrders, simulate_transmission
from exoatmosim.poseidon import TransmissionSpectrum


class TestMockIntegration(unittest.TestCase):
    def test_shared_instrument_library(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'library').mkdir()
            np.savetxt(root / 'library' / 'orders.dat', [[750, 751]])
            (root / 'planet').mkdir()
            path = root / 'planet' / 'config'
            path.write_text('instrument_library_path\t../library\ninstrument_order_file\torders.dat\nmock_snr\t80\n')
            settings = ObservationSettings.from_config(path)
            self.assertEqual(settings.order_file, (root/'library'/'orders.dat').resolve())
            self.assertEqual(len(SpectralOrders.load(settings.order_file).ranges_nm), 1)

    def test_both_eccentric_and_circular_approximate_obs(self):
        with TemporaryDirectory() as tmp:
            cfg_path = Path(tmp) / 'config'
            table = Path(tmp) / 'orders.dat'
            np.savetxt(table, [[750., 750.2]])
            for ecc in (0., 0.6):
                config = Config.from_mapping({'P': 2, 'Tc': 2460000, 'ecc': ecc,
                    'omega': 30, 'K': 0.1, 'Ms': 1., 'Mp': 1.,
                    'instrument_order_file': 'orders.dat', 'mock_snr': 100,
                    'mock_resolution': 100000, 'mock_seed': 5})
                config.write(cfg_path)
                times = 2460000 + np.linspace(-0.04, 0.04, 8)
                rv = RVTraceEstimator.from_config(config).evaluate(times)
                self.assertTrue(np.isfinite(rv['rv_planet_kms']).all())
                wave = np.linspace(0.7, 0.9, 10000)
                absorption = 0.01 + 0.001 * np.exp(-((wave-0.7501)/0.00005)**2)
                mock = simulate_transmission(
                    TransmissionSpectrum(wave, absorption),
                    SpectralOrders.load(table), ObservationSettings.from_config(cfg_path),
                    rv['phase'], rv['rv_planet_kms'])
                self.assertEqual(mock.flux[0].shape[0], len(times))
                self.assertTrue(np.isfinite(mock.flux[0]).all())
                self.assertTrue(np.allclose(mock.uncertainty[0], .01))

    def test_plot_rejects_nonfinite_mock(self):
        from exoatmosim.observations import SyntheticObservation
        with TemporaryDirectory() as tmp:
            obs = SyntheticObservation((np.array([750., 751.]),),
                (np.array([[np.nan, 1.]]),), (np.array([[.01, .01]]),),
                np.array([0.]), np.array([0.]))
            with self.assertRaisesRegex(ValueError, 'non-finite flux'):
                obs.plot(Path(tmp) / 'bad.png')

    def test_order_rejects_placeholder(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / 'invalid.npy'
            path.write_bytes(b'hi')
            with self.assertRaises((ValueError, OSError, EOFError)):
                SpectralOrders.load(path)
