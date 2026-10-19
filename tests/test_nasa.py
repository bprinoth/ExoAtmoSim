import unittest
import importlib.util
from unittest.mock import patch
from pathlib import Path
from tempfile import TemporaryDirectory

from exoatmosim import Config
from exoatmosim.nasa import parameters_from_row, find_parameters

EXPECTED_KEYS = (
    'P a aRstar Rp Mp K vsini RpRstar vsys inclination Tc Tp T14 RA DEC '
    'Teq Teff ecc impact omega Rs Ms logg metallicity stellar_type Hmag Kmag '
    'Jmag Gyr_age loggp lampoo resolution observatory air transit_flag lat '
    'long elev observation_window_hours readout exposure_time'
).split()
ATMOSPHERE_KEYS = ('tp_profile chemistry condensation p0 FastChem_input '
                   'network_path layers kappa_IR gamma T_int C/O Ti/H').split()


@unittest.skipUnless(importlib.util.find_spec('astropy') and importlib.util.find_spec('astroquery'),
                     'Install optional NASA extras: pip install -e .[nasa]')
class TestArchiveCompatibility(unittest.TestCase):
    def test_legacy_key_set_and_roundtrip(self):
        row = {'pl_orbper': 1.27, 'pl_bmassj': 1.3, 'pl_radj': 1.8,
               'ra': 180.0, 'dec': -20.0, 'pl_trandur': 2.3,
               'st_spectype': 'M2.5V', 'tran_flag': 1}
        mapping = parameters_from_row(row)
        self.assertEqual(list(mapping), EXPECTED_KEYS)
        original = Config.from_mapping(mapping)
        for key in ATMOSPHERE_KEYS:
            original.set(key, 'some_value')
        text = original.to_text()
        reread = Config.from_text(text)
        self.assertEqual(reread.to_text(), text)
        self.assertEqual(set(reread), set(EXPECTED_KEYS + ATMOSPHERE_KEYS))
        self.assertEqual(reread['stellar_type'], 'M2.5V')
        self.assertEqual(reread['transit_flag'], 'True')
        self.assertAlmostEqual(reread.get_float('P'), 1.27)
        self.assertAlmostEqual(reread.get_float('observation_window_hours'), 4.3)

    def test_archive_query_and_output(self):
        row = {'pl_orbper': 1.27, 'pl_bmassj': 1.3, 'pl_radj': 1.8,
               'ra': 180.0, 'dec': -20.0, 'pl_trandur': 2.3}
        with TemporaryDirectory() as tmp:
            # Stub only the network boundary; everything else runs for real.
            with patch('astroquery.ipac.nexsci.nasa_exoplanet_archive.NasaExoplanetArchive.query_criteria',
                       return_value=[row]) as query:
                config = find_parameters('WASP-121 b', print_params=False, data_dir=tmp)
                query.assert_called_once_with(table='pscomppars', select='*',
                                              where="pl_name='WASP-121 b'")
            path = Path(tmp) / 'WASP-121_b' / 'config'
            self.assertEqual(Config.from_file(path).to_text(), config.to_text())
            self.assertTrue(set(EXPECTED_KEYS).issubset(set(config)))
            self.assertEqual(config['chemistry'], 'fastchem')
            self.assertEqual(config['instrument'], 'ESPRESSO')
            with self.assertRaises(FileExistsError):
                with patch('astroquery.ipac.nexsci.nasa_exoplanet_archive.NasaExoplanetArchive.query_criteria',
                           return_value=[row]):
                    find_parameters('WASP-121 b', print_params=False, data_dir=tmp)


if __name__ == '__main__':
    unittest.main()
