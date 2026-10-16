import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
import numpy as np
from exoatmosim.tayph_mock import run_tayph_mock

class MockIntegratedRV:
    def __init__(self, path, observatory):
        cfg=(Path(path)/'config').read_text()
        entries=dict(line.split(None,1) for line in cfg.splitlines() if line.strip())
        assert observatory == 'Paranal'
        # ExoAtmoSim preserves the source configuration and lets the integrated
        # estimator apply the circular omega convention in memory.
        self.e=float(entries['ecc'])
        self.omega=float(entries['omega'])
        self.bjds=np.loadtxt(Path(path)/'obs_times', ndmin=1)+2400000.5
        self.true_anomaly=np.linspace(-.2,.2,len(self.bjds))
    def plot_doppler_trace(self, RV_ext=200):
        import matplotlib.pyplot as plt
        assert RV_ext==200
        plt.figure()
        plt.plot(self.true_anomaly,self.true_anomaly)

class TestIntegratedBridge(unittest.TestCase):
    def _case(self,e):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);d=root/'data'/'WASP-189_b';d.mkdir(parents=True)
            config=d/'config'
            config.write_text(f'instrument ESPRESSO\nP 2.724\nTc 2458926.5\necc {e}\nomega {"nan" if e==0 else "35"}\nK .182\nMs 2.03\nMp 1.99\n')
            (d/'WASP-189_b.fits').write_bytes(b'fits placeholder')
            bjds=2458926.5+np.linspace(-.025,.025,6)
            mock=Mock()
            module=types.ModuleType('exoatmosim.rvtrace')
            module.RVTraceEstimator=MockIntegratedRV
            old=os.getcwd()
            try:
                os.chdir(root)
                with patch.dict(sys.modules, {'exoatmosim.rvtrace':module}):
                    prepared,kwargs=run_tayph_mock(config,project_root=root,bjds=bjds,mock_module=mock)
            finally:
                os.chdir(old)
            self.assertTrue(kwargs['eccentric_orbit'])
            self.assertFalse(kwargs['tellurics'])
            np.testing.assert_allclose(kwargs['bjds'],bjds,atol=1e-8,rtol=0)
            np.testing.assert_allclose(kwargs['phase_range'],np.linspace(-.2,.2,len(bjds)))
            self.assertTrue((root/'diagnostics'/'WASP-189_b'/'rv_doppler_trace.png').is_file())
            mock.create_mock_obs.assert_called_once()
            self.assertEqual(mock.create_mock_obs.call_args.kwargs['mode'],'ESPRESSO')
            self.assertTrue(prepared.runfile.exists())
            self.assertIn('omega nan' if e==0 else 'omega 35',config.read_text())
    def test_telluric_flag_in_dry_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); d=root/'data'/'p'; d.mkdir(parents=True)
            config=d/'config'
            config.write_text('instrument ESPRESSO\ntellurics True\n')
            (d/'p.fits').write_bytes(b'fake')
            _, kwargs=run_tayph_mock(config,project_root=root,dry_run=True)
            self.assertTrue(kwargs['tellurics'])

    def test_circular(self):self._case(0.)
    def test_eccentric(self):self._case(.3)

if __name__=='__main__': unittest.main()
