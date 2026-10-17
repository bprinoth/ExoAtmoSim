import tempfile
import unittest
from pathlib import Path
import numpy as np
from exoatmosim import AbundanceProfile, Atmosphere, Config
from exoatmosim.converters import KalteneggerConverter, FetukinesConverter


class TestAbundanceProfiles(unittest.TestCase):
    def test_isochemical_roundtrip_and_bulk_species(self):
        pressures = np.logspace(-6, 0, 7)
        a = AbundanceProfile(pressures, {
            'H2': np.full(7, 0.85),
            'He': np.full(7, 0.149),
            'H2O': np.full(7, 0.001),
        }, ('H2', 'He'))
        with tempfile.TemporaryDirectory() as tmp:
            path = a.to_file(Path(tmp) / 'abundances.dat')
            b = AbundanceProfile.from_file(path, bulk_species='H2 He')
            self.assertEqual(b.bulk_species, ('H2', 'He'))
            for species in a.mixing_ratios:
                np.testing.assert_allclose(a.mixing_ratios[species], b.mixing_ratios[species])

    def test_missing_bulk_rejected(self):
        with self.assertRaisesRegex(ValueError, 'bulk_species'):
            AbundanceProfile([0.1, 1], {'H2O': [0.1, 0.1]}, ('H2',))

    def test_total_above_one_rejected(self):
        with self.assertRaisesRegex(ValueError, 'exceeds one'):
            AbundanceProfile([0.1, 1], {'N2': [0.9, 0.9], 'O2': [0.3, 0.3]}, ('N2',))

    def test_interp_log_pressure_and_refuse_extrapolation(self):
        a = AbundanceProfile([1e-6, 1e-2], {'H2': [0.8, 0.6]}, ('H2',))
        interpolated = a.on_atmosphere(Atmosphere([1e-6, 1e-4, 1e-2], [600, 600, 600]))
        np.testing.assert_allclose(interpolated.mixing_ratios['H2'], [0.8, 0.7, 0.6])
        with self.assertRaisesRegex(ValueError, 'extrapolate'):
            a.on_atmosphere(Atmosphere([1e-7, 1e-4], [600, 600]))

    def test_kaltenegger_converter_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / 'pt').write_text('Z,T,P 0 300 1000\nZ,T,P 10 250 100\n')
            (d / 'chem').write_text('Molecule N2\nN2_MR 0.8, ^\n 0.7\nN2_ALT 0, 10\nMolecule O2\nO2_MR 0.2, 0.3\nO2_ALT 0, 10\n')
            converter = KalteneggerConverter()
            a, b = converter.convert(d / 'pt', d / 'chem', bulk_species='N2')
            np.testing.assert_allclose(a.pressure_bar, [0.1, 1.0])
            np.testing.assert_allclose(b.mixing_ratios['N2'], [0.7, 0.8])
            _, out = converter.export(d / 'pt', d / 'chem', bulk_species='N2', output_dir=d / 'output')
            loaded = AbundanceProfile.from_file(out, bulk_species='N2')
            np.testing.assert_allclose(loaded.mixing_ratios['O2'], [0.3, 0.2])

    def test_fetukines_converter_sort_all_columns_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            source = d / 'network.dat'
            source.write_text('! altitude temperature[K] pressure[bar] MMW[g] N2 O2\n'
                              '0 300 1.0 -23.6 -0.096910013  -0.698970004\n'
                              '10 250 0.1 -23.6 -0.154901960 -0.522878745\n')
            converter = FetukinesConverter()
            atm, abundances = converter.convert(source, bulk_species='N2')
            np.testing.assert_allclose(atm.pressure_bar, [0.1, 1])
            np.testing.assert_allclose(atm.temperature_k, [250, 300])
            np.testing.assert_allclose(abundances.mixing_ratios['N2'], [0.7, 0.8], rtol=1e-7)
            _, out = converter.export(source, bulk_species='N2', output_dir=d / 'out')
            restored = AbundanceProfile.from_file(out, bulk_species='N2')
            np.testing.assert_allclose(restored.mixing_ratios['O2'], [0.3, 0.2], rtol=1e-7)

    def test_plot_optional(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest('optional matplotlib unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            a = AbundanceProfile([1e-5, 1e-3], {'N2': [0.9, 0.9]}, ('N2',))
            path = a.plot(Path(tmp) / 'diagnostics' / 'abundances.png')
            self.assertTrue(path.is_file())

if __name__ == '__main__':
    unittest.main()
