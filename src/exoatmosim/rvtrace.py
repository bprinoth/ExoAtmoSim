"""Integrated RVTraceEstimator based on the original user-supplied implementation.

Requires astropy, radvel, PyAstronomy and matplotlib for full observational traces.
"""
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from PyAstronomy import pyasl
import astropy.units as u
import radvel
from PyAstronomy import modelSuite as ms




# def occulted_fraction(d, Rp, zp):
#     """
#     Compute the fraction of the stellar disk occulted by a planet.
    
#     Parameters
#     ----------
#     d : float
#         Projected separation between centers (in units of stellar radius).
#     Rp : float
#         Planet radius (in units of stellar radius).
#     zp : float
#         Planet position in z direction (in units of stellar radius)
    
#     Returns
#     -------
#     f : float
#         Fraction of stellar disk occulted by the planet.
#     """
#     if d >= 1 + Rp or zp < 0: # outside of stellar disc, or behind star!
#         return 0.0  # No overlap
#     elif d <= abs(1 - Rp):
#         return Rp**2  # Full overlap (planet entirely within star)
#     else:
#         # Partial overlap
#         acos_arg1 = (d**2 + Rp**2 - 1) / (2 * d * Rp)
#         acos_arg2 = (d**2 + 1 - Rp**2) / (2 * d)

#         # Clamp values to [-1, 1] to avoid numerical issues
#         acos_arg1 = np.clip(acos_arg1, -1, 1)
#         acos_arg2 = np.clip(acos_arg2, -1, 1)

#         part1 = Rp**2 * np.arccos(acos_arg1)
#         part2 = np.arccos(acos_arg2)

#         # Compute the area of overlap using the circle-circle intersection formula
#         part3 = 0.5 * np.sqrt(
#             (-d + Rp + 1) *
#             (d + Rp - 1) *
#             (d - Rp + 1) *
#             (d + Rp + 1)
#         )

#         A = part1 + part2 - part3
#         return A / np.pi



# import numpy as np

from exoatmosim.observing.geometry import occulted_fraction
from exoatmosim.config import Config
from exoatmosim.observing.obs_times import read_obs_times_mjd

class RVTraceEstimator:
    
    
    def __init__(self, dp, obs='paranal'):
        
        self.dp = dp
        # The estimator accepts a planet directory, not a external mock dependency.
        self.config = Config.from_file(Path(dp) / 'config')
        self.obs = obs
        self.omega = None  # Define these variables in __init__ or set them elsewhere
        self.transitC = None
        self.period = None
        self.ecc = None
        self.bjds = None
        self.aRs = None
        self.Ms = None
        self.Mp = None
        self.transit = None
        self.RpRs = None
        self.T14 = None
        
        self.load_observation()
        self.read_orbital_configuration()
        self.calculate_planet_position()
        self.calculate_berv()
        self.RV_RM()
        self.RV_star()
        self.RV_planet()
        
    def _parameter(self, key):
        """Read legacy numeric orbital parameters from ExoAtmoSim's config.

        RA/DEC may be sexagesimal strings; all other inputs used by this
        estimator are physical numeric quantities.
        """
        if key in ('RA', 'DEC'):
            return self.config[key]
        return self.config.get_float(key)

    def load_observation(self):
        print('[INFO] Loading observation')
           
        # obstimes need to be in MJDs
        
        obs_file = Path(self.dp) / 'obs_times'
        if obs_file.exists():
            obstimes = read_obs_times_mjd(obs_file)
        else:
            print("[INFO] Obstime doesn't exist, I will create obstimes for you.")
            obstimes = []
        
        if len(obstimes) == 0: # so in case this file is empty, we need to produce them from scratch
            self.transitC = self._parameter('Tc')
            self.period = self._parameter('P')
            self.T14 = self.config.get_float('T14') if 'T14' in self.config else self.config.get_float('t', 4.0) 
            
            obstimes = np.linspace(self.transitC-(self.T14)/24/2, self.transitC+(self.T14)/24/2, 2000) - 2400000.5
            
        
        self.obstimes=obstimes
        self.mjds = Time(obstimes, format='mjd')
        self.bjds = self.mjds.jd
        
    
    def transit_ecc(self):
        print('[INFO] Determining in-transit expsoures')
        
        self.calculate_planet_position()
        
        # rho is the distance between the *centres* of the two discs.
        # Passing rho - RpRs would shift the overlap boundaries and is wrong.
        rho = np.hypot(self.xp, self.yp)
        obscured = np.array([occulted_fraction(rho[i], self.RpRs, self.zp[i])
                             for i in range(len(rho))])
        transit = 1.0 - obscured
        if not np.all(np.isfinite(transit)) or np.any((transit < 0) | (transit > 1)):
            raise ValueError('Non-physical transit light curve from orbital geometry')

        self.transit = transit        
    
    def read_orbital_configuration(self):
        print('[INFO] Reading orbital configuration')
        #self.load_observation()
        self.RA = self._parameter('RA')
        self.DEC = self._parameter('DEC')
        self.a = self._parameter('a')
        self.aRs = self._parameter('aRstar')
        self.orbinc = self._parameter('inclination')
        self.vsini = self._parameter('vsini')
        self.pob = self._parameter('lampoo')
        self.transitC = self._parameter('Tc')
        self.period = self._parameter('P')
        self.ecc = self._parameter('ecc')
        self.RpRs = self._parameter('RpRstar')
        self.K = self._parameter('K')
        self.vsys = self._parameter('vsys')        
        self.omega = self._parameter('omega')  # degrees; circular periastron undefined
        if self.ecc == 0.0 and not np.isfinite(self.omega):
            self.omega = 90.0
        elif not np.isfinite(self.omega):
            raise ValueError('A finite omega is required for eccentric orbits')
        omega_bar = np.radians(self.omega)
        self.Ms = self._parameter('Ms') * u.Msun
        self.Mp = self._parameter('Mp') * u.Mjup
        self.T_per = radvel.orbit.timetrans_to_timeperi(self.transitC, self.period, self.ecc, omega_bar)
        self.transit_ecc()          
    
    def calculate_planet_position(self):
        
        print('[INFO] Calculate planet position')
            
        # degree to radians
        omega_bar = np.radians(self.omega)
        inclin_bar = np.radians(self.orbinc)


        ke = pyasl.KeplerEllipse(
            a=self.aRs,
            per=self.period,
            e=self.ecc,
            i=self.orbinc,  # in degrees
            Omega=0.,  # We don't care about the orientation in the sky
            w=self.omega,  # in degrees,
            tau=self.T_per
        )

        ecc_anomaly = ke.eccentricAnomaly(self.bjds)
        node = np.radians(self.pob)

        # Use arctan2 for quadrant safety, then unwrap along observation time.
        # A 2*pi wrap is not a physical jump in radial velocity.
        E = np.asarray(ecc_anomaly, float)
        anomaly = 2. * np.arctan2(np.sqrt(1. + self.ecc) * np.sin(E / 2.),
                                 np.sqrt(1. - self.ecc) * np.cos(E / 2.))
        self.true_anomaly = np.unwrap(anomaly)

        # in matrix form

        self.r0_p = self.aRs * np.array([
                                   (np.cos(ecc_anomaly) - self.ecc),
                                   np.sqrt(1. - self.ecc * self.ecc) * np.sin(ecc_anomaly),
                                    np.zeros_like(ecc_anomaly)
                                   ])


        rotation_angle = omega_bar - np.pi/2
        r1_rotation_matrix =  np.array([
                                              [np.cos(rotation_angle), -np.sin(rotation_angle), 0],
                                              [np.sin(rotation_angle), np.cos(rotation_angle), 0],
                                              [0,0,1]
                                   ])


        # # create the eccentric orbit
        # self.X0_p = self.r0_p[0]
        # self.Y0_p = self.r0_p[1]

        self.r0_1 = r1_rotation_matrix @ self.r0_p

        self.X1_p = self.r0_1[0]
        self.Y1_p = self.r0_1[1]

        # # turn it w.r.t the argument of periastron
        # self.X1_p = self.X0_p * np.sin(omega_bar) + self.Y0_p * np.cos(omega_bar)
        # self.Y1_p = -self.X0_p * np.cos(omega_bar) + self.Y0_p * np.sin(omega_bar)

        # turn it w.r.t to the inclination

        r_pl_rotation_matrix = np.array([
                    [0, 1, 0],
                    [-np.cos(inclin_bar), 0, 0],
                    [np.sin(inclin_bar), 0, 0]

        ])

        self.r_pl = r_pl_rotation_matrix @ self.r0_1

        rp_rotation_matrix = np.array([
                    [np.cos(node), -np.sin(node), 0],
                    [np.sin(node), np.cos(node), 0],
                    [0, 0, 1]
        ])

        # x_pl = self.Y1_p
        # y_pl = -self.X1_p * np.cos(inclin_bar)
        # z_pl = self.X1_p * np.sin(inclin_bar)

        self.rp = rp_rotation_matrix @ self.r_pl

        self.xp = self.rp[0]
        self.yp = self.rp[1]
        self.zp = self.rp[2]

        # self.xp = x_pl * np.cos(node) - y_pl * np.sin(node)
        # self.yp = x_pl * np.sin(node) + y_pl * np.cos(node)
        # self.zp = z_pl
        self.ecc_anomaly =  ecc_anomaly
        

    def calculate_berv(self):
        print('[INFO] Calculating BERV')
        berv = []
        
    
        # Armazones is not available in every Astropy site registry. Prefer
        # explicit config coordinates for this site, and use the same values
        # for RV corrections as for the observation schedule.
        site = str(self.obs).strip().lower()
        if site in ('armazones', 'cerro armazones'):
            observatory = EarthLocation.from_geodetic(
                lon=self.config.get_float('long') * u.deg,
                lat=self.config.get_float('lat') * u.deg,
                height=self.config.get_float('elev') * u.m)
        else:
            try:
                observatory = EarthLocation.of_site(self.obs)
            except Exception as exc:
                if all(key in self.config for key in ('lat', 'long', 'elev')):
                    import warnings
                    warnings.warn(
                        f'Astropy site {self.obs!r} unavailable; using explicit '
                        'config lat/long/elev for barycentric correction.',
                        RuntimeWarning)
                    observatory = EarthLocation.from_geodetic(
                        lon=self.config.get_float('long') * u.deg,
                        lat=self.config.get_float('lat') * u.deg,
                        height=self.config.get_float('elev') * u.m)
                else:
                    raise ValueError(
                        f'Unknown observatory {self.obs!r}; provide finite '
                        'lat, long and elev in the planet config.') from exc
        
        try:
            sc = SkyCoord(self.RA+' '+self.DEC, unit=(u.hourangle, u.deg))
        except:
            sc = SkyCoord(ra=self.RA * u.deg, dec=self.DEC*u.deg)
        for date in self.obstimes :
            barycorr = sc.radial_velocity_correction(obstime=Time(date,format='mjd'), location=observatory).to(u.km/u.s)
            berv.append(barycorr.value)
            
        self.berv = np.asarray(berv)
        
        
    def RV_RM(self):
        print('[INFO] Calculating the RV extent of the RM effect')
        self.calculate_planet_position()
        self.vstar = self.xp * self.vsini
        
        
    def RV_star(self):
        #calculate_planet_position(self)
        print('[INFO] Calculating the RV extent of the star')
        
        rv = ms.KeplerRVModel()  # Changed from ms.KeplerRVModel()
        rv.assignValue({
            "per1": self.period,
            'K1': self.K,
            'e1': self.ecc,
            'tau1': self.T_per,
            'w1': self.omega,
            'mstar': self.Ms,
            'c0': 0.,
            'a1': self.a,
            'msini1': self.Mp * np.sin(np.radians(self.orbinc))
        })

        RVs = rv.evaluate(self.bjds)  # evaluate the RV of the star at a given time
        self.RVstar = RVs
            
    def RV_planet(self):
        print('[INFO] Calculating the RV extent of the planet')
        #calculate_planet_position(self)
        
        self.RV_star() 
        self.RVplanet = -1. * (self.RVstar) * (self.Ms / self.Mp).decompose()
              
    def plot_doppler_trace(self, RF='system', RV_ext=100, save=True):
        """Plot velocity tracks against transit-centered orbital phase.

        The original chart divided true anomaly by 2*pi and called it phase.
        In eccentric orbits, true anomaly at transit generally is NOT zero.
        The right axis reports the corresponding true anomaly in degrees.
        """
        from exoatmosim.plotting import publication_style

        if RF == 'all':
            for frame in ('star', 'planet', 'system', 'obs', 'berv'):
                self.plot_doppler_trace(RF=frame, RV_ext=RV_ext, save=save)
            return

        time = np.asarray(self.bjds, float)
        order = np.argsort(time)
        time = time[order]
        phase = ((time - self.transitC) / self.period + 0.5) % 1.0 - 0.5
        # The observing window should be shorter than one orbit; prevent
        # rendering artificial lines across the phase-wrap discontinuity.
        if np.any(np.diff(phase) < 0):
            raise ValueError('Observation window crosses the orbital phase wrap; '
                             'plot the transit-centered observing window instead')
        true_anomaly = np.unwrap(np.asarray(self.true_anomaly, float)[order])
        anomaly_deg = np.degrees(true_anomaly)
        # Choose a continuous branch near 180 deg for conjunction, instead
        # of jumping between -180 and +180 around the plotted transit.
        anomaly_deg += 360. * np.round((180. - np.median(anomaly_deg)) / 360.)
        self.true_phases = np.asarray(self.true_anomaly)[order] / (2 * np.pi)
        self.obtimes = np.asarray(self.obstimes)[order]
        rv_star = np.asarray(self.RVstar, float)[order]
        rv_planet = np.asarray(self.RVplanet, float)[order]
        rm = np.asarray(self.vstar, float)[order] + rv_star
        berv = np.asarray(self.berv, float)[order]
        telluric_rv = -self.vsys + berv

        if RF == 'star':
            curves = (rm-rv_star, rv_planet-rv_star, np.zeros_like(rv_star))
        elif RF == 'planet':
            curves = (rm-rv_planet, np.zeros_like(rv_planet), rv_star-rv_planet)
        elif RF == 'obs':
            curves = (rm-telluric_rv, rv_planet-telluric_rv,
                      rv_star-telluric_rv)
        elif RF == 'berv':
            curves = (rm+self.vsys, rv_planet+self.vsys,
                      rv_star+self.vsys)
        else:
            curves = (rm, rv_planet, rv_star)
            if RF != 'system':
                print(f'[WARNING] Unknown rest frame {RF!r}; using system')

        self.RVs_RM, self.RVs_planet, self.RVs_star = curves
        self.RVs_tel = telluric_rv
        inside = np.asarray(self.transit, float)[order] < 1.0

        with publication_style():
            fig, ax = plt.subplots(figsize=(7.2, 6.2), layout='constrained')
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            ax.plot(curves[1], phase, color=colors[1], lw=2.3, label='Planet')
            ax.plot(curves[2], phase, color=colors[2], lw=1.7, label='Star')
            # Mask RM outside the transit so it cannot connect separate segments.
            ax.plot(np.where(inside, curves[0], np.nan), phase, color=colors[0],
                    lw=1.8, ls='--', label='RM')
            for index in np.flatnonzero(np.diff(inside.astype(int))):
                ax.axhline(phase[index], color='0.55', ls=':', lw=1)
            ax.axhline(0, color='0.5', ls=':', lw=1)
            ax.set(xlabel='RV (km/s)', ylabel='Orbital phase relative to transit',
                   title=f'{RF.capitalize()} rest frame')
            if RV_ext is not None and RV_ext > 0:
                ax.set_xlim(-RV_ext, RV_ext)
            ax.legend(loc='best', fontsize=10)
            right = ax.twinx()
            right.set_ylim(ax.get_ylim())
            right.set_yticks(ax.get_yticks())
            right.set_yticklabels([
                f'{np.interp(tick, phase, anomaly_deg):.1f}'
                if phase[0] <= tick <= phase[-1] else ''
                for tick in ax.get_yticks()
            ])
            right.set_ylabel('True anomaly (deg)')
            if save:
                fig.savefig(f'{self.dp}/RVTraceEstimator_{RF}.png', dpi=180)
            plt.close(fig)
            return fig
