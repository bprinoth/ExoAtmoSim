import unittest
from types import SimpleNamespace
from exoatmosim.config import Config


class EffectiveTransitCenterTests(unittest.TestCase):
    def test_max_altitude_overrides_only_staged_config(self):
        from examples.run_independent_crires import staged_orbital_config
        source = Config.from_mapping({'Tc': 2459000.0, 'P': 3.2})
        stage = staged_orbital_config(source, SimpleNamespace(mode='max_altitude', center_jd_utc=2459012.34))
        self.assertAlmostEqual(stage.get_float('Tc'), 2459012.34)
        self.assertAlmostEqual(source.get_float('Tc'), 2459000)

    def test_tobs_does_not_override_ephemeris(self):
        from examples.run_independent_crires import staged_orbital_config
        source = Config.from_mapping({'Tc': 2459000.0})
        stage = staged_orbital_config(source, SimpleNamespace(mode='Tobs', center_jd_utc=2459012.34))
        self.assertEqual(stage.get_float('Tc'), source.get_float('Tc'))


class PlotStyleTests(unittest.TestCase):
    def test_context_restores_original_style(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest('matplotlib not installed')
        from exoatmosim.plotting import publication_style
        previous = matplotlib.rcParams['axes.linewidth']
        with publication_style():
            self.assertEqual(matplotlib.rcParams['axes.linewidth'], 3)
            self.assertTrue(matplotlib.rcParams['xtick.top'])
        self.assertEqual(matplotlib.rcParams['axes.linewidth'], previous)
