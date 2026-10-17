"""Regression checks for the public modular interface and example files."""
import ast
import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from exoatmosim import cli, workflow

ROOT = Path(__file__).resolve().parents[1]


class ModularInterfaceTests(unittest.TestCase):
    def test_stage_help_is_available(self):
        for stage in ('config', 'atmosphere', 'spectrum', 'observe', 'simulate'):
            with self.subTest(stage=stage):
                with redirect_stdout(io.StringIO()) as output:
                    with self.assertRaises(SystemExit) as exc:
                        cli.main([stage, '--help'])
                self.assertEqual(exc.exception.code, 0)
                self.assertIn('usage:', output.getvalue())

    def test_config_routes_to_workflow(self):
        with patch.object(workflow, 'make_config', return_value=Path('config')) as fn:
            with redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(['config', '--planet', 'Test b', '--instrument', 'CRIRES_K']), 0)
            self.assertEqual(fn.call_args.kwargs['instrument'], 'CRIRES_K')

    def test_observe_missing_model_does_not_recompute(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as directory:
            planet = Path(directory) / 'Test_b'
            planet.mkdir()
            config = planet / 'config'
            config.write_text('instrument CRIRES_K\n')
            with self.assertRaisesRegex(FileNotFoundError, 'run `exoatmosim spectrum` first'):
                workflow.simulate_observation(config, seed=71)

    def test_examples_are_lightweight_public_api_calls(self):
        names = ('crires_k', 'espresso', 'igrins', 'igrins2', 'spirou', 'andes', 'andes_k')
        for name in names:
            with self.subTest(name=name):
                tree = ast.parse((ROOT / 'examples' / (name + '.py')).read_text())
                imports = [node for node in tree.body if isinstance(node, ast.ImportFrom)]
                self.assertTrue(any(node.module == 'exoatmosim' and any(alias.name == 'simulate' for alias in node.names)
                                    for node in imports))
