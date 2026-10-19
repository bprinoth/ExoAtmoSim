"""Kaltenegger data path belongs to machine-local config, not NASA planet files."""
import tempfile
import unittest
from pathlib import Path
from exoatmosim.nasa import load_site_config, format_science_config


class KalteneggerSiteKeyRegression(unittest.TestCase):
    def test_load_site_path_without_copying_to_planet(self):
        with tempfile.TemporaryDirectory() as tmp:
            site_file = Path(tmp) / "site_paths.conf"
            site_file.write_text(
                "kaltenegger_profiles_dir /some/external/Kaltenegger_profiles\n",
                encoding="utf-8",
            )
            paths = load_site_config(site_file)
            self.assertEqual(paths["kaltenegger_profiles_dir"],
                             "/some/external/Kaltenegger_profiles")
            generated = format_science_config(
                {"stellar_type": "G2V", "Gyr_age": 4.6},
                site_values=paths,
                include_workflow=False,
            )
            self.assertIn("stellar_type", generated)
            self.assertNotIn("kaltenegger_profiles_dir", generated)


if __name__ == "__main__":
    unittest.main()
