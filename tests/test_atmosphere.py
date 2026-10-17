import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from exoatmosim import Config, Atmosphere


class TestAtmosphere(unittest.TestCase):
    def test_isothermal_defaults(self):
        a = Atmosphere.from_config(Config.from_mapping({'Teq': 1200, 'layers': 8}))
        self.assertEqual(len(a.pressure_bar), 8)
        np.testing.assert_allclose(a.temperature_k, 1200)
        self.assertAlmostEqual(a.pressure_bar[0], 1e-9)
        self.assertAlmostEqual(a.pressure_bar[-1], 10)

    def test_explicit_isothermal_temperature(self):
        a = Atmosphere.from_config(Config.from_mapping({'Teq': 1200, 'T_iso': 900, 'tp_profile': 'isothermal'}))
        np.testing.assert_allclose(a.temperature_k, 900)

    def test_file_sorted_and_relative_to_config(self):
        with tempfile.TemporaryDirectory() as d:
            Path(d, 'my_tp.txt').write_text('500 1\n800 0.01\n700 0.1\n')
            a = Atmosphere.from_config(Config.from_mapping({'tp_profile': 'my_tp.txt'}), config_dir=d)
            np.testing.assert_allclose(a.pressure_bar, [0.01, 0.1, 1])
            np.testing.assert_allclose(a.temperature_k, [800, 700, 500])

    def test_reject_bad_profile(self):
        with self.assertRaises(ValueError):
            Atmosphere([1, 1], [500, 600])
        with self.assertRaises(ValueError):
            Atmosphere.from_config(Config.from_mapping({'Teq': 500, 'layers': 1}))

    def test_plot_writes_diagnostic_to_requested_path(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest('matplotlib not installed')
        with tempfile.TemporaryDirectory() as d:
            a = Atmosphere([1e-3, 1e-2, 1], [400, 500, 600])
            path = Path(d) / 'diagnostics' / 'Test_b' / 'tp_profile.png'
            self.assertEqual(a.plot_tp(path), path)
            self.assertTrue(path.is_file())
            self.assertGreater(path.stat().st_size, 100)

    def test_guillot_calls_prt_with_cgs_gravity(self):
        fake = types.ModuleType('petitRADTRANS.physics')
        calls = []
        def guillot(p, kappa, gamma, g, ti, teq):
            calls.append((kappa, gamma, g, ti, teq))
            return np.full_like(p, 1010.)
        fake.guillot_global = guillot
        with patch.dict(sys.modules, {'petitRADTRANS': types.ModuleType('petitRADTRANS'), 'petitRADTRANS.physics': fake}):
            a = Atmosphere.from_config(Config.from_mapping({'tp_profile': 'guillot', 'Teq': 1200, 'loggp': 3.0, 'layers': 3}))
        self.assertEqual(len(calls), 1)
        self.assertAlmostEqual(calls[0][2], 1000)
        np.testing.assert_allclose(a.temperature_k, 1010)


if __name__ == '__main__':
    unittest.main()
