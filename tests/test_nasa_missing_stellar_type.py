import unittest
from exoatmosim.nasa import format_science_config
from exoatmosim.config import Config


class TestMissingStellarType(unittest.TestCase):
    def test_blank_stellar_type_from_archive_does_not_break_config(self):
        # The NASA parser should map empty st_spectype to nan before rendering.
        from exoatmosim.nasa import _value
        self.assertIn('stellar_type\tnan', format_science_config({'stellar_type': 'nan'}, include_workflow=False))
