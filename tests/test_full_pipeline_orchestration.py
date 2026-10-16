
"""Regression tests for the FastChem → POSEIDON → CRIRES+ pipeline."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))

from run_full_crires_pipeline import resolve_config


class FullFlowTests(unittest.TestCase):

    def test_existing_config_not_rewritten(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "config"
            original = "instrument\tCRIRES_K\nchemistry\tfastchem\n"
            p.write_text(original)

            self.assertEqual(
                resolve_config(p, planet=None, site_config=None),
                p.resolve(),
            )
            self.assertEqual(p.read_text(), original)

    def test_conflicting_config_arguments_rejected(self):
        with self.assertRaisesRegex(ValueError, "OR"):
            resolve_config(
                Path("config"),
                planet="TOI-1842 b",
                site_config=None,
            )

    def test_new_planet_config_uses_crires_species(self):
        with tempfile.TemporaryDirectory() as d:
            import os

            previous = os.getcwd()
            try:
                os.chdir(d)

                def fake_nasa(planet, **kwargs):
                    p = Path("data") / "Test_b" / "config"
                    p.parent.mkdir(parents=True)
                    p.write_text(
                        "instrument\tCRIRES_K\n"
                        "chemistry\tmanual\n"
                    )

                with patch(
                    "exoatmosim.nasa.find_parameters",
                    side_effect=fake_nasa,
                ):
                    output = resolve_config(
                        None,
                        planet="Test b",
                        site_config=None,
                    )

                text = output.read_text()
                self.assertIn("chemistry\tfastchem", text)
                self.assertIn("fastchem_species\tH2 He H2O", text)

            finally:
                os.chdir(previous)


if __name__ == "__main__":
    unittest.main()
