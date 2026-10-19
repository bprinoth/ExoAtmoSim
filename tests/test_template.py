"""Tests for isolated species transmission templates (no POSEIDON installation required)."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from exoatmosim import cli, workflow
from exoatmosim.profiles import AbundanceProfile


class TemplateTests(unittest.TestCase):
    def test_cli_routes_species(self):
        with patch.object(workflow, 'calculate_template', return_value=Path('O3.fits')) as fn:
            self.assertEqual(cli.main(['template', 'config', '--species', 'O3']), 0)
            self.assertEqual(fn.call_args.kwargs['species'], 'O3')

    def test_retains_bulk_and_target_without_other_opacities(self):
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / 'planet' / 'config'
            config.parent.mkdir()
            config.write_text('instrument CRIRES_K\n')
            p = np.array([1e-5, 1.0])
            abundance = AbundanceProfile(p, {
                'H2': np.array([0.8, 0.8]),
                'He': np.array([0.19, 0.19]),
                'O3': np.array([1e-7, 2e-7]),
                'H2O': np.array([1e-5, 1e-5]),
            }, ('H2', 'He'))
            recorded = {}
            saved = {}

            def compute_transmission(**kwargs):
                recorded.update(kwargs)
                return SimpleNamespace(wavelength_um=np.array([1., 1.1]),
                                       transmission=np.array([0.99, 0.98]))

            run = SimpleNamespace(abundances=abundance, model=SimpleNamespace(
                compute_transmission=compute_transmission), name='planet',
                atmosphere=object(), config=object())
            def capture_fits(spectrum, destination, **kwargs):
                saved['spectrum'] = spectrum
                return destination

            with patch('exoatmosim.poseidon.TransmissionRun.from_config', return_value=run), \
                 patch('exoatmosim.spectrum_fits.write_spectrum_fits', side_effect=capture_fits), \
                 patch('exoatmosim.poseidon.TransmissionSpectrum.plot', return_value=Path('O3.png')) as plot:
                result = workflow.calculate_template(config, species='O3')
                self.assertEqual(result.name, 'O3.fits')
                np.testing.assert_allclose(saved['spectrum'].transmission, [0.01, 0.02])
                self.assertEqual(plot.call_args.args[0].name, 'O3.png')
                self.assertEqual(plot.call_args.kwargs['ylabel'], 'Transit depth (T)')
                self.assertEqual(result.parent.name, 'templates')
                composition = recorded['composition'].mixing_ratios
                self.assertEqual(set(composition), {'H2', 'He', 'O3'})
                np.testing.assert_array_equal(composition['O3'], abundance.mixing_ratios['O3'])
                self.assertEqual(recorded['bulk_species'], ('H2', 'He'))
                with self.assertRaisesRegex(ValueError, 'bulk species'):
                    workflow.calculate_template(config, species='H2')
                with self.assertRaisesRegex(ValueError, 'not in the atmospheric abundances'):
                    workflow.calculate_template(config, species='CH4')

    def test_existing_template_requires_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / 'planet' / 'config'
            config.parent.mkdir()
            config.write_text('instrument CRIRES_K\n')
            folder = config.parent / 'simulations' / 'CRIRES_K' / 'seed_0042' / 'templates'
            folder.mkdir(parents=True)
            (folder / 'O3.fits').touch()
            run = SimpleNamespace(abundances=SimpleNamespace(
                mixing_ratios={'H2': np.ones(2) * .9, 'O3': np.ones(2) * 1e-6},
                bulk_species=('H2',)))
            with patch('exoatmosim.poseidon.TransmissionRun.from_config', return_value=run):
                with self.assertRaises(FileExistsError):
                    workflow.calculate_template(config, species='O3')
                with self.assertRaisesRegex(ValueError, 'valid single species'):
                    workflow.calculate_template(config, species='../O3')


if __name__ == '__main__':
    unittest.main()
