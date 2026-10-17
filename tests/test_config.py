import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from exoatmosim import Config


EXISTING = ("# legacy-compatible config\n"
            "P\t1.2749255\nRp  1.8\n"
            "RA\t07:12:34.56\nDEC\t+01:02:03.40\n"
            "air\tTrue\nexposure_time\t150\n"
            "stellar_type\tM3A\n"
            "custom_passthrough_key\tkeep_me\n\n")


class TestConfig(unittest.TestCase):
    def test_round_trip_preserves_all_lines(self):
        config = Config.from_text(EXISTING)
        self.assertEqual(config.to_text(), EXISTING)
        self.assertEqual(config['custom_passthrough_key'], 'keep_me')

    def test_typed_access(self):
        config = Config.from_text(EXISTING)
        self.assertAlmostEqual(config.get_float('P'), 1.2749255)
        self.assertEqual(config.get_int('exposure_time'), 150)
        self.assertTrue(config.get_bool('air'))
        self.assertEqual(config['RA'], '07:12:34.56')
        self.assertEqual(config.get_float('missing', 2.0), 2.0)

    def test_edits_preserve_unrelated_lines(self):
        config = Config.from_text(EXISTING)
        config.set('exposure_time', 180)
        config.set('chemistry', 'FastChem')
        self.assertIn('Rp  1.8\n', config.to_text())
        self.assertIn('custom_passthrough_key\tkeep_me\n', config.to_text())
        self.assertEqual(config.get_int('exposure_time'), 180)
        self.assertEqual(config['chemistry'], 'FastChem')

    def test_file_round_trip(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'data' / 'WASP-121_b' / 'config'
            config = Config.from_text(EXISTING)
            config.write(path)
            self.assertEqual(path.read_text(), EXISTING)
            self.assertEqual(Config.from_file(path).to_text(), EXISTING)

    def test_reject_ambiguous_or_invalid_input(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            Config.from_text('P 1\nP 2\n')
        with self.assertRaisesRegex(ValueError, 'Invalid'):
            Config.from_text('P\n')
        with self.assertRaisesRegex(ValueError, 'boolean'):
            Config.from_text('air maybe\n').get_bool('air')

    def test_mapping_serialization(self):
        self.assertEqual(Config.from_mapping({'air': True, 'P': 1.2}).to_text(),
                         'air\tTrue\nP\t1.2\n')


if __name__ == '__main__':
    unittest.main()
