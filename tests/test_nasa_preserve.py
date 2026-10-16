import unittest
from exoatmosim.config import Config

class TestNASAPreservation(unittest.TestCase):
    def test_existing_config_updated_without_losing_custom_settings(self):
        original = Config.from_text(
            '# Legacy header\nP  1.0\n# Custom FastChem settings\n'
            'fastchem_elements_file\t/path/abundances.dat\n'
            'bulk_species\tH2 He\n')
        original.set('P', 2.0)
        self.assertIn('# Custom FastChem settings', original.to_text())
        self.assertEqual(original['fastchem_elements_file'], '/path/abundances.dat')
        self.assertEqual(original['bulk_species'], 'H2 He')
