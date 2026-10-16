import tempfile
from pathlib import Path
import unittest
import numpy as np
from exoatmosim.observing.obs_times import read_obs_times_mjd


class TestObsTimesReader(unittest.TestCase):
    def check(self, text):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'obs_times'
            path.write_text(text)
            return read_obs_times_mjd(path)

    def test_headerless_scheduler_format(self):
        np.testing.assert_array_equal(self.check('59000.5\n59000.6\n'), [59000.5, 59000.6])

    def test_tayph_multicolumn(self):
        np.testing.assert_array_equal(self.check('59000.5 2020-01-01 150\n59000.6 2020-01-02 150\n'), [59000.5, 59000.6])

    def test_named_header(self):
        np.testing.assert_array_equal(self.check('MJD texp\n59000.5 150\n59000.6 150\n'), [59000.5, 59000.6])

    def test_bad_content(self):
        for content in ('', 'header only\n', '59000.5\nnan\n', '59000.5\n59000.5\n', '59000.6\n59000.5\n', '59000.5\nbad\n'):
            with self.subTest(content=content), self.assertRaises(ValueError):
                self.check(content)

if __name__ == '__main__': unittest.main()
