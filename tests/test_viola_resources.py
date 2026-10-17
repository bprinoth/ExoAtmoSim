import csv
import tempfile
import unittest
from pathlib import Path
import numpy as np
from exoatmosim.instruments import canonical_instrument, instrument_resource
from exoatmosim.observing.viola_resources import export_viola_order_grid

class ViolaResourceTests(unittest.TestCase):
    def test_bands_are_registered_but_unmeasured_data_is_not_invented(self):
        for band in ('K', 'H', 'Y'):
            self.assertEqual(canonical_instrument('VIOLA_'+band), 'VIOLA_'+band)

    def test_real_csv_becomes_single_segment_grid_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = Path(tmp)/'signal_noise.csv'
            with csv_path.open('w', newline='') as fh:
                writer = csv.writer(fh)
                writer.writerow(['wavelength_um', 'stellar_e_per_resel', 'snr_per_resel'])
                writer.writerows([(1.02, 2, 10), (1.0201, 3, 11), (1.0202, 4, 12)])
            out = export_viola_order_grid(csv_path, 'Y', tmp)
            np.testing.assert_allclose(np.load(out), [[1020., 1020.1, 1020.2]])
            self.assertEqual(out, export_viola_order_grid(csv_path, 'Y', tmp))
            with self.assertRaises(FileExistsError):
                csv_path.write_text('wavelength_um,stellar_e_per_resel,snr_per_resel\n1.1,2,3\n1.2,3,4\n')
                export_viola_order_grid(csv_path, 'Y', tmp)

if __name__ == '__main__': unittest.main()
