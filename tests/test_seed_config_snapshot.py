"""Regression tests for config snapshots in seed-numbered output."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from exoatmosim.full_crires import run


class SeedConfigSnapshotTests(unittest.TestCase):
    def test_existing_snapshot_with_changed_config_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'Test_b'
            root.mkdir()
            source = root / 'config'
            source.write_text('instrument CRIRES_K\nchemistry fastchem\n')
            snapshot = root / 'simulations' / 'CRIRES_K' / 'seed_0042' / 'config'
            snapshot.parent.mkdir(parents=True)
            snapshot.write_text('instrument CRIRES_K\nchemistry fastchem\ntellurics True\n')
            with self.assertRaisesRegex(FileExistsError, 'different configuration'):
                run(source, seed=42)

    def test_existing_run_cannot_be_overwritten_implicitly(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'Test_b'
            root.mkdir()
            source = root / 'config'
            source.write_text('instrument CRIRES_K\nchemistry fastchem\n')
            out = root / 'simulations' / 'CRIRES_K' / 'seed_0042'
            out.mkdir(parents=True)
            (out / 'config').write_bytes(source.read_bytes())
            (out / 'order_0.fits').write_bytes(b'placeholder')
            with self.assertRaises(FileExistsError):
                run(source, seed=42)


if __name__ == '__main__':
    unittest.main()
