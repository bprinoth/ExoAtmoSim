import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from exoatmosim import TransmissionRun, TransmissionSpectrum


class TestStandardPoseidonWorkflow(unittest.TestCase):
    def build(self, directory, *, bulk='CO2 O2', abundance_max_pressure=1):
        d = Path(directory)
        (d / 'tp.dat').write_text('300 0.01\n310 0.1\n320 1.0\n')
        (d / 'abundances.dat').write_text(
            'pressure CO2 O2 H2O\n'
            f'0.01 0.78 0.21 0.01\n'
            f'0.1 0.78 0.21 0.01\n'
            f'{abundance_max_pressure} 0.78 0.21 0.01\n')
        (d / 'config').write_text(
            'Rp\t1\nMp\t1\nTeq\t300\nRs\t1\nTeff\t5700\nlogg\t4.4\n'
            'metallicity\t0\n'
            'tp_profile\ttp.dat\nchemistry\tfile\nabundances_file\tabundances.dat\n'
            f'bulk_species\t{bulk}\nresolution\t100000\nwl_min\t0.7\nwl_max\t1.0\n')
        return d / 'config'

    def test_dry_run_preserves_original_config_and_two_bulk_gases(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.build(temp)
            before = path.read_bytes()
            run = TransmissionRun.from_config(path)
            self.assertEqual(run.abundances.bulk_species, ('CO2', 'O2'))
            self.assertEqual(run.model.wavelength_min_um, 0.7)
            self.assertEqual(len(run.atmosphere.pressure_bar), 3)
            self.assertEqual(before, path.read_bytes())

    def test_missing_bulk_species_is_error(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.build(temp, bulk='H2')
            with self.assertRaisesRegex(ValueError, 'bulk_species'):
                TransmissionRun.from_config(path)

    def test_out_of_range_abundance_grid_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.build(temp, abundance_max_pressure=0.5)
            with self.assertRaisesRegex(ValueError, 'extrapolate'):
                TransmissionRun.from_config(path)

    def test_compute_can_be_stubbed_and_saved_without_poseidon(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.build(temp)
            run = TransmissionRun.from_config(path)
            fake = TransmissionSpectrum(np.array([0.7, 0.8]), np.array([0.01, 0.02]))
            with patch.object(type(run.model), 'compute_transmission', return_value=fake) as compute:
                result = run.compute()
            self.assertIs(result, fake)
            self.assertEqual(compute.call_args.kwargs['bulk_species'], ('CO2', 'O2'))
            a, b = run.save(result, directory=Path(temp) / 'output')
            np.testing.assert_allclose(np.load(a), [0.7, 0.8])
            np.testing.assert_allclose(np.load(b), [0.01, 0.02])


if __name__ == '__main__':
    unittest.main()
