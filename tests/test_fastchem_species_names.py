import unittest
import sys
from types import SimpleNamespace
from unittest.mock import patch
from exoatmosim.fastchem import resolve_fastchem_species_index, fastchem_hill_name


class FakeEngine:
    def __init__(self, names):
        self.names = {name: i for i, name in enumerate(names)}

    def getGasSpeciesIndex(self, name):
        return self.names.get(name, 2**32-1)


class TestFastChemSpeciesNames(unittest.TestCase):
    def setUp(self):
        # Stub the dependency in this offline test environment while exercising
        # the real element sorting, charge formatting and database lookup.
        import re
        class StubFormula:
            def __init__(self, formula, charge=0):
                self.element = {}
                tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
                if ''.join(a+b for a,b in tokens) != formula:
                    raise ValueError(formula)
                for a,b in tokens:
                    self.element[a] = self.element.get(a, 0) + int(b or 1)
        self.patch = patch.dict(sys.modules, {'chemformula': SimpleNamespace(ChemFormula=StubFormula)})
        self.patch.start()
        self.addCleanup(self.patch.stop)

    def test_hill_notation_general(self):
        self.assertEqual(fastchem_hill_name('H2O'), 'H2O1')
        self.assertEqual(fastchem_hill_name('TiO'), 'O1Ti1')
        self.assertEqual(fastchem_hill_name('CO2'), 'C1O2')
        self.assertEqual(fastchem_hill_name('Fe+'), 'Fe1+')
        self.assertEqual(fastchem_hill_name('Ti2+'), 'Ti1++')
        self.assertEqual(fastchem_hill_name('e-'), 'e-')

    def test_isomer_requires_exact_identifier(self):
        with self.assertRaisesRegex(ValueError, 'isomer-ambiguous'):
            resolve_fastchem_species_index(FakeEngine(['C1H1N1']), 'HNC', 1)

    def test_neutral_and_charged_metals(self):
        engine = FakeEngine(['Fe', 'Fe1+', 'Ti1+', 'e-', 'H1-', 'H'])
        self.assertEqual(resolve_fastchem_species_index(engine, 'Fe', 6), (0, 'Fe'))
        self.assertEqual(resolve_fastchem_species_index(engine, 'Fe+', 6), (1, 'Fe1+'))
        self.assertEqual(resolve_fastchem_species_index(engine, 'Ti+', 6), (2, 'Ti1+'))
        self.assertEqual(resolve_fastchem_species_index(engine, 'e-', 6), (3, 'e-'))
        self.assertEqual(resolve_fastchem_species_index(engine, 'H-', 6), (4, 'H1-'))

    def test_molecular_formula_to_fastchem_stoichiometry(self):
        engine = FakeEngine(['H2O1', 'O1Ti1', 'H2', 'He'])
        self.assertEqual(resolve_fastchem_species_index(engine, 'H2O', 4), (0, 'H2O1'))
        self.assertEqual(resolve_fastchem_species_index(engine, 'TiO', 4), (1, 'O1Ti1'))

    def test_missing_molecule_reports_backend_candidate(self):
        with self.assertRaisesRegex(ValueError, 'H2O1'):
            resolve_fastchem_species_index(FakeEngine(['H2', 'He']), 'H2O', 2)

    def test_h_minus_never_falls_back_to_electron(self):
        with self.assertRaisesRegex(ValueError, 'H-'):
            resolve_fastchem_species_index(FakeEngine(['e-']), 'H-', 1)

    def test_missing_ion_raises(self):
        with self.assertRaisesRegex(ValueError, r'Fe1\+'):
            resolve_fastchem_species_index(FakeEngine(['Fe']), 'Fe+', 1)

if __name__ == '__main__':
    unittest.main()
