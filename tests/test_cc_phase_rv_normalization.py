"""Display-only phase-RV normalization regression tests."""
import unittest
import numpy as np
from exoatmosim.cc_plots import normalize_phase_rv_rows


class PhaseRVNormalizationTests(unittest.TestCase):
    def test_offsets_and_scales_do_not_create_stripes(self):
        x = np.linspace(-2, 2, 101)
        profile = np.exp(-x*x / 0.12)
        raw = np.vstack([profile, 50 + 4*profile, -30 + .2*profile])
        normalized = normalize_phase_rv_rows(raw)
        np.testing.assert_allclose(np.mean(normalized, axis=1), 0., atol=1e-12)
        np.testing.assert_allclose(np.std(normalized, axis=1, ddof=1), 1., atol=1e-12)
        np.testing.assert_allclose(normalized[0], normalized[1], atol=1e-12)
        np.testing.assert_allclose(normalized[0], normalized[2], atol=1e-12)
        np.testing.assert_allclose(raw[1], 50 + 4*profile)

    def test_flat_rows_masked(self):
        x = np.vstack([np.arange(20., dtype=float), np.ones(20)])
        result = normalize_phase_rv_rows(x)
        self.assertTrue(np.isfinite(result[0]).all())
        self.assertTrue(np.isnan(result[1]).all())
