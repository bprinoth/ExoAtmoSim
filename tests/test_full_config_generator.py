"""Regression tests for the consolidated NASA config generator."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from exoatmosim.config import Config
from exoatmosim.nasa import format_science_config, load_site_config


class GeneratorTests(unittest.TestCase):
    def test_generates_workflow_without_deprecated_mock(self):
        cfg = Config.from_text(format_science_config({'P': 2.7, 'ecc': 0, 'omega': 90}))
        self.assertEqual(cfg['tp_profile'], 'isothermal')
        self.assertEqual(cfg['chemistry'], 'fastchem')
        self.assertEqual(cfg['fastchem_condensation'], 'False')
        self.assertEqual(cfg['rv_plot_extent'], '200')
        self.assertEqual(cfg['tellurics'], 'False')
        self.assertEqual(cfg['pwv'], '2.5')
        self.assertEqual(cfg['pwv_sigma'], '0.5')
        self.assertFalse(any(k.startswith('telfit_') for k in cfg))
        self.assertIn('Observing atmosphere and tellurics',
                      format_science_config({'P': 2.7}))
        self.assertLess(format_science_config({'P': 2.7}).index('# --- Orbital velocities'),
                        format_science_config({'P': 2.7}).index('# --- Observing atmosphere'))
        self.assertEqual(cfg['P'], '2.7')
        self.assertNotIn('mock_snr', cfg)
        self.assertNotIn('mock_n_exposures', cfg)
        self.assertNotIn('rvtrace_path', cfg)
        self.assertNotIn('poseidon_path', cfg)

    def test_site_keys(self):
        with TemporaryDirectory() as d:
            path = Path(d) / 'site.conf'
            path.write_text('poseidon_input_data /data/POSEIDON_inputs/inputs\n'
                            'pysyn_cdbs /data/POSEIDON_inputs/inputs/stellar_grids\n')
            cfg = Config.from_text(format_science_config({'P': 1}, site_values=load_site_config(path)))
            self.assertEqual(cfg['poseidon_input_data'], '/data/POSEIDON_inputs/inputs')
            path.write_text('mock_snr 10\n')
            with self.assertRaises(ValueError):
                load_site_config(path)

    def test_instrument_override_preserves_all_science(self):
        cfg = Config.from_text(format_science_config({'P': 3, 'ecc': 0.3, 'omega': 210}))
        cfg.set('instrument', 'CRIRES_K')
        self.assertEqual(cfg['ecc'], '0.3')
        self.assertEqual(cfg['omega'], '210')
        self.assertEqual(cfg['instrument'], 'CRIRES_K')

if __name__ == '__main__':
    unittest.main()
