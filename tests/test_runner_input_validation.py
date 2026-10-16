import unittest
import numpy as np
from exoatmosim.config import Config

class InputDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        # The runner is an example script, import by adding examples to sys.path.
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples'))
        from run_independent_crires import validate_observation_inputs
        self.validate = validate_observation_inputs
        self.cfg = Config.from_mapping({'Mp': 'nan', 'Ms': 1, 'K': 0.1})
    def test_nonfinite_planet_velocity_identified(self):
        with self.assertRaisesRegex(ValueError, 'Non-finite planet RV.*Mp=nan'):
            self.validate(np.array([0., 0.]), np.array([0., np.nan]), [1., 0.99], self.cfg)
    def test_nonfinite_star_velocity_identified(self):
        with self.assertRaisesRegex(ValueError, 'Non-finite stellar RV.*K=0.1'):
            self.validate(np.array([np.nan]), np.array([0.]), [1.], self.cfg)
    def test_no_transit_is_not_invalid_data(self):
        self.validate(np.array([0.]), np.array([0.]), [1.], self.cfg)
    def test_invalid_transit_geometry(self):
        with self.assertRaisesRegex(ValueError, 'Transit light curve outside'):
            self.validate(np.array([0.]), np.array([0.]), [1.5], self.cfg)
