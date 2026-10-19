import tempfile
import unittest
from pathlib import Path
import numpy as np
from exoatmosim.observing.viola_plots import plot_viola_diagnostics


class ViolaDiagnosticPlotTests(unittest.TestCase):
    def test_creates_plots_without_radiative_transfer(self):
        with tempfile.TemporaryDirectory() as temp:
            w = np.linspace(2000, 2001, 40)
            star = np.tile(1200. - 10. * np.sin(w), (3, 1))
            planet = star.copy()
            planet[1] *= .99
            one, two = plot_viola_diagnostics(temp, w, star, planet,
                    np.linspace(1999, 2002, 100), np.ones(100)*.99,
                    np.array([0, 1, 0]))
            self.assertTrue(one.is_file())
            self.assertTrue(two.is_file())
            self.assertGreater(one.stat().st_size, 1000)
            self.assertGreater(two.stat().st_size, 1000)

    def test_mismatched_axes_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(ValueError):
                plot_viola_diagnostics(temp, np.arange(10), np.ones((2, 9)),
                       np.ones((2, 9)), np.arange(10), np.ones(10), np.ones(2))
