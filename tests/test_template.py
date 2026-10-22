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

    def test_preserves_full_composition_for_isolated_opacity(self):
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
                self.assertEqual(set(composition), {'H2', 'He', 'O3', 'H2O'})
                self.assertEqual(recorded['opacity_species'], 'O3')
                np.testing.assert_array_equal(composition['H2O'], abundance.mixing_ratios['H2O'])
                np.testing.assert_array_equal(composition['O3'], abundance.mixing_ratios['O3'])
                self.assertEqual(recorded['bulk_species'], ('H2', 'He'))
                with self.assertRaisesRegex(ValueError, 'bulk species'):
                    workflow.calculate_template(config, species='H2')
                with self.assertRaisesRegex(ValueError, 'not in the atmospheric abundances'):
                    workflow.calculate_template(config, species='CH4')

    def test_spectrum_create_templates_uses_one_shared_model_call(self):
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / 'planet' / 'config'
            config.parent.mkdir()
            config.write_text('instrument CRIRES_K\n')
            p = np.array([1e-5, 1.])
            abundance = AbundanceProfile(p, {
                'H2': np.array([.8, .8]),
                'He': np.array([.19, .19]),
                'H2O': np.array([1e-5, 1e-5]),
                'O3': np.array([1e-7, 1e-7]),
            }, ('H2', 'He'))
            full = SimpleNamespace(wavelength_um=np.array([1., 1.1]),
                                   transmission=np.array([.992, .991]), plot=lambda *a, **k: None)
            water = SimpleNamespace(wavelength_um=np.array([1., 1.1]),
                                    transmission=np.array([.992, .990]))
            ozone = SimpleNamespace(wavelength_um=np.array([1., 1.1]),
                                    transmission=np.array([.991, .989]))
            call = []
            def compute(**kwargs):
                call.append(kwargs)
                return full, {'H2O': water, 'O3': ozone}
            run = SimpleNamespace(abundances=abundance, name='planet', config=object(),
                                  atmosphere=object(), save=lambda *a, **k: None, model=SimpleNamespace(
                                      compute_transmission=compute))
            written = {}
            def save_fits(spectrum, destination, **kwargs):
                written[Path(destination).name] = np.asarray(spectrum.transmission)
                return destination
            with patch('exoatmosim.poseidon.TransmissionRun.from_config', return_value=run), \
                 patch('exoatmosim.spectrum_fits.write_spectrum_fits', side_effect=save_fits), \
                 patch('exoatmosim.poseidon.TransmissionSpectrum.plot'), \
                 patch('exoatmosim.workflow.publication_style'):
                workflow.calculate_spectrum(config, create_templates=True)
            self.assertEqual(len(call), 1)
            self.assertEqual(call[0]['template_species'], ('H2O', 'O3'))
            self.assertEqual(set(written), {'planet.fits', 'H2O.fits', 'O3.fits'})
            np.testing.assert_allclose(written['H2O.fits'], [.008, .010])

    def test_spectrum_cli_create_templates(self):
        with patch.object(workflow, 'calculate_spectrum', return_value=Path('planet.fits')) as fn:
            self.assertEqual(cli.main(['spectrum', 'config', '--create-templates']), 0)
            self.assertTrue(fn.call_args.kwargs['create_templates'])

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
