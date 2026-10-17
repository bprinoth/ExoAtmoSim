"""Unit-level TelFit injection tests, independent of installed LBLRTM."""
import unittest
from types import SimpleNamespace
import numpy as np
from exoatmosim.observing.tellurics import zenith_angles, telfit_transmission, air_to_vacuum_nm, vacuum_to_air_nm
from exoatmosim.observing.order_mock import simulate_orders


class TelluricTests(unittest.TestCase):
    def test_real_zenith_angles(self):
        np.testing.assert_allclose(zenith_angles([90, 60, 30]), [0, 30, 60])
        with self.assertRaisesRegex(ValueError, 'visible'):
            zenith_angles([60, -2])

    def test_model_called_for_each_actual_exposure(self):
        calls = []
        class Modeler:
            def MakeModel(self, **params):
                calls.append(params)
                x = np.linspace(1948., 1962., 100)
                return SimpleNamespace(x=x, y=np.full(x.size, 1-params['angle']/200.))
        schedule = SimpleNamespace(altitude_deg=np.array([75., 45.]))
        grid, transmissions = telfit_transmission(schedule, [1950., 1960.], modeler_factory=Modeler)
        self.assertEqual(transmissions.shape, (2,2))
        self.assertEqual([p['angle'] for p in calls], [15.,45.])
        self.assertTrue(all(p['vac2air'] for p in calls))

    def test_air_vacuum_round_trip_k_band(self):
        air = np.array([1925., 1930., 1935., 2470.])
        vac = air_to_vacuum_nm(air)
        np.testing.assert_allclose(vacuum_to_air_nm(vac), air, rtol=0, atol=1e-9)
        self.assertTrue(0.50 < vac[2] - air[2] < 0.55)

    def test_telfit_bounds_are_vacuum_and_cover_air_grid(self):
        calls = []
        class Modeler:
            def MakeModel(self, **kwargs):
                calls.append(kwargs)
                # Mimics TelFit: vacuum bounds -> air-nm output.
                vacuum = np.linspace(1e7 / kwargs['highfreq'],
                                     1e7 / kwargs['lowfreq'], 500)
                wavelength_air = vacuum_to_air_nm(vacuum)
                return SimpleNamespace(x=wavelength_air,
                                       y=np.ones(wavelength_air.size))
        schedule = SimpleNamespace(altitude_deg=np.array([65., 60.]))
        air = np.linspace(1925., 1935., 50)
        w, transmissions = telfit_transmission(schedule, air, modeler_factory=Modeler)
        np.testing.assert_array_equal(w, air)
        self.assertEqual(len(calls), 2)
        self.assertEqual(transmissions.shape, (2, 50))
        np.testing.assert_allclose(vacuum_to_air_nm(1e7 / calls[0]['highfreq']),
                                   air[0] - max(.02, 2 * np.median(np.diff(air))), atol=1e-8)

    def test_small_telfit_ringing_is_clipped_with_warning(self):
        class Modeler:
            def MakeModel(self, **kwargs):
                x = np.linspace(1948., 1962., 100)
                y = np.ones_like(x)
                y[3] = -1.865762169e-4
                return SimpleNamespace(x=x, y=y)
        schedule = SimpleNamespace(altitude_deg=np.array([75.]))
        with self.assertWarns(RuntimeWarning):
            _, t = telfit_transmission(schedule, [1950., 1960.], modeler_factory=Modeler)
        self.assertTrue(np.all((t >= 0) & (t <= 1)))
        with self.assertWarns(RuntimeWarning):
            _, strict = telfit_transmission(schedule, [1950., 1960.], modeler_factory=Modeler,
                                             ringing_tolerance=0)
        self.assertTrue(np.all((strict >= 0) & (strict <= 1)))

    def test_large_telfit_negative_is_clipped(self):
        class Modeler:
            def MakeModel(self, **kwargs):
                x = np.linspace(1948., 1962., 100)
                y = np.ones_like(x)
                y[3] = -.001
                return SimpleNamespace(x=x, y=y)
        schedule = SimpleNamespace(altitude_deg=np.array([75.]))
        with self.assertWarns(RuntimeWarning):
            _, clipped = telfit_transmission(schedule, [1950., 1960.], modeler_factory=Modeler)
        self.assertTrue(np.all((clipped >= 0) & (clipped <= 1)))

    def test_above_unity_beyond_tolerance_is_rejected(self):
        class Modeler:
            def MakeModel(self, **kwargs):
                x = np.linspace(1948., 1962., 100)
                y = np.ones_like(x)
                y[3] = 1.002
                return SimpleNamespace(x=x, y=y)
        schedule = SimpleNamespace(altitude_deg=np.array([75.]))
        with self.assertRaisesRegex(ValueError, 'outside'):
            telfit_transmission(schedule, [1950., 1960.], modeler_factory=Modeler)

    def test_telluric_application_pre_noise(self):
        waves = np.linspace(1950, 1960, 101)
        orders = np.tile(waves[10:90], (1,1))
        out = simulate_orders(orders,waves,np.ones_like(waves),waves,np.ones_like(waves),
                              [0.,0.], [0.,0.], [0.,0.],add_noise=False,
                              telluric_wavelength_nm=waves,
                              telluric_transmission=np.array([np.full(101,.9),np.full(101,.7)]))
        np.testing.assert_allclose(out[0][1][:,0], [0.9,0.7])

    def test_missing_telluric_input_rejected(self):
        w = np.linspace(1950, 1960, 101)
        with self.assertRaisesRegex(ValueError, 'both'):
            simulate_orders(w[None,:],w,w,w,w,[0.],[0.],[0.],
                            telluric_wavelength_nm=w)

if __name__=='__main__': unittest.main()
