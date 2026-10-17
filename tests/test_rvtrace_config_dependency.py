"""Guard against reintroducing the undefined legacy previous mock framework paramget helper."""
import unittest
import tempfile
from pathlib import Path

from exoatmosim.config import Config


class TestRVConfigMigration(unittest.TestCase):
    def test_all_rvtrace_parameters_use_local_reader(self):
        source = (Path(__file__).resolve().parents[1] / 'src' / 'exoatmosim' / 'rvtrace.py').read_text()
        self.assertNotIn('paramget(', source)
        self.assertIn('Config.from_file', source)

    def test_typed_config_and_sexagesimal_coordinates(self):
        config = Config.from_text('Tc 2459000.0\nP 2.5\nRA 12:13:14\nDEC -45:00:00\n')
        self.assertEqual(config.get_float('Tc'), 2459000.)
        self.assertEqual(config.get_float('P'), 2.5)
        self.assertEqual(config['RA'], '12:13:14')

    def test_estimator_gets_numeric_and_coordinates_if_dependencies_available(self):
        try:
            from exoatmosim.rvtrace import RVTraceEstimator
        except ImportError as exc:
            self.skipTest(f'RV orbit dependencies unavailable: {exc}')
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, 'config').write_text('Tc 2459000\nP 2.5\nRA 12:13:14\nDEC -45:00:00\n')
            est = RVTraceEstimator.__new__(RVTraceEstimator)
            est.config = Config.from_file(Path(tmp, 'config'))
            self.assertEqual(est._parameter('Tc'), 2459000.)
            self.assertEqual(est._parameter('P'), 2.5)
            self.assertEqual(est._parameter('RA'), '12:13:14')


if __name__ == '__main__':
    unittest.main()
