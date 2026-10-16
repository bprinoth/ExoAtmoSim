import math
import unittest
from exoatmosim import Config, Star, Planet

class TestPhysicalObjects(unittest.TestCase):
    def test_existing_config_to_star_and_planet(self):
        config = Config.from_text('Rs 1.5\nMs 1.2\nTeff 6500\nstellar_type F7V\n'
                                  'Rp 1.8\nMp 1.3\nTeq 2400\nloggp 3.0\np0 0.01\n'
                                  'P 1.27\nRA 12:00:00\ncustom_tayph_key kept\n')
        star, planet = Star.from_config(config), Planet.from_config(config)
        self.assertEqual(star.radius_solar, 1.5)
        self.assertEqual(star.spectral_type, 'F7V')
        self.assertEqual(planet.reference_pressure_bar, 0.01)
        self.assertTrue(math.isfinite(planet.calculated_logg_cgs))
        self.assertAlmostEqual(planet.gravity_discrepancy_dex,
                               3.0 - planet.calculated_logg_cgs)
        self.assertEqual(config['custom_tayph_key'], 'kept')
        self.assertEqual(config['RA'], '12:00:00')

    def test_missing_nan_and_invalid_physics(self):
        star = Star.from_config(Config.from_text('Rs nan\nTeff 5700\n'))
        self.assertIsNone(star.radius_solar)
        planet = Planet.from_config(Config.from_text('Rp 1.2\n'))
        self.assertIsNone(planet.calculated_logg_cgs)
        with self.assertRaisesRegex(ValueError, 'radius_jupiter'):
            Planet.from_config(Config.from_text('Rp -1\n'))

if __name__ == '__main__':
    unittest.main()
