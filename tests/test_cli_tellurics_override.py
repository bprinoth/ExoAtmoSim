"""Telluric CLI overrides must never change the reusable planet config."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from exoatmosim.config import Config
from exoatmosim.full_crires import run


class TelluricOverrideTests(unittest.TestCase):
    def test_override_is_effective_and_template_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'config'
            original = 'instrument CRIRES_K\ntellurics False\n'
            source.write_text(original)
            def fake_pipeline(config_path, **kwargs):
                self.assertTrue(Config.from_file(config_path).get_bool('tellurics'))
                self.assertEqual(Path(config_path).parent, source.parent)
                self.assertEqual(source.read_text(), original)
                return 'ok'
            with patch('exoatmosim.full_crires._run_pipeline', side_effect=fake_pipeline):
                self.assertEqual(run(source, tellurics=True), 'ok')
            self.assertEqual(source.read_text(), original)
            self.assertEqual(list(Path(tmp).glob('.exoatmosim_effective_*')), [])

    def test_default_uses_original_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'config'
            source.write_text('tellurics False\n')
            with patch('exoatmosim.full_crires._run_pipeline', return_value='ok') as runner:
                self.assertEqual(run(source), 'ok')
                self.assertEqual(runner.call_args.args[0], source.resolve())
