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

        self.true_anomaly = 2. * np.arctan(np.sqrt((1. + self.ecc) / (1. - self.ecc)) * np.tan(ecc_anomaly / 2.))

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
        
    
        observatory = EarthLocation.of_site(self.obs)
        
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
        #print('[INFO] Calculating the RV extent of the star')
        
        # plotting in different restframes
        
        true_anomaly_sorted_idx = np.argsort(self.true_anomaly)
        true_anomaly_sorted = np.take_along_axis(self.true_anomaly, true_anomaly_sorted_idx, axis=0)
        RVstar_sorted = np.take_along_axis(self.RVstar, true_anomaly_sorted_idx, axis=0)
        vstar_sorted = np.take_along_axis(self.vstar, true_anomaly_sorted_idx, axis=0)  + RVstar_sorted
        vplanet_sorted = np.take_along_axis(self.RVplanet, true_anomaly_sorted_idx, axis=0)
        berv_sorted = np.take_along_axis(self.berv, true_anomaly_sorted_idx, axis=0)
        tel = - self.vsys + berv_sorted
        self.true_phases = true_anomaly_sorted / (2 * np.pi)
        transit_sorted =  np.take_along_axis(self.transit, true_anomaly_sorted_idx, axis=0)
        self.obtimes = np.take_along_axis(self.obstimes, true_anomaly_sorted_idx, axis=0)
            
        #else:
        # true_anomaly_sorted = self.true_anomaly
        # RVstar_sorted = self.RVstar
        # vstar_sorted = self.vstar.T +  RVstar_sorted
        # vplanet_sorted = self.RVplanet
        # berv_sorted = self.berv
        # tel = - self.vsys + berv_sorted
        # true_phases = self.true_phases * 2 * np.pi
        # transit_sorted =  self.transit
        
        
        if RF=='star':
            
            self.RVs_RM = vstar_sorted - RVstar_sorted
            self.RVs_planet = vplanet_sorted - RVstar_sorted
            self.RVs_star = RVstar_sorted - RVstar_sorted
            self.RVs_tel = tel - RVstar_sorted
            
            print(f"[INFO] Plot in stellar restframe")
            
        
        elif RF=='planet':
            self.RVs_RM = vstar_sorted - vplanet_sorted
            self.RVs_planet = vplanet_sorted - vplanet_sorted
            self.RVs_star = RVstar_sorted - vplanet_sorted
            self.RVs_tel = tel - vplanet_sorted
            
            print(f"[INFO] Plot in planetary restframe")
        
        
        elif RF=='obs':
            self.RVs_RM = vstar_sorted - tel
            self.RVs_planet = vplanet_sorted - tel
            self.RVs_star = RVstar_sorted - tel
            self.RVs_tel = tel - tel
            
            print(f"[INFO] Plot in observatory restframe")
        
        elif RF=='berv':
      
            self.RVs_RM = vstar_sorted + self.vsys
            self.RVs_planet = vplanet_sorted + self.vsys
            self.RVs_star = RVstar_sorted + self.vsys
            self.RVs_tel = tel + self.vsys
            
            print(f"[INFO] Plot in berv restframe")

        elif RF=='system':
            self.RVs_RM = vstar_sorted
            self.RVs_planet = vplanet_sorted
            self.RVs_star = RVstar_sorted
            self.RVs_tel = tel

            print(f"[INFO] Plot in system restframe")
            
        
        
        elif RF=='all':
            
            print(f"[INFO] Plotting in all restframes")
            
            self.plot_doppler_trace(RF='star')
            self.plot_doppler_trace(RF='planet')
            self.plot_doppler_trace(RF='system')
            self.plot_doppler_trace(RF='obs')
            self.plot_doppler_trace(RF='berv')
            
        else:
            self.RVs_RM = vstar_sorted.T - RVstar_sorted
            self.RVs_planet = vplanet_sorted
            self.RVs_star = RVstar_sorted
            self.RVs_tel = tel

            print(f"[INFO] I don't know this restframe, I am plotting in the system restframe")


        
        #RV_ext = np.max(np.abs(np.concatenate()))
        
        if RF != 'all':
            plt.figure(figsize=(6,6))
            plt.title(f'Restframe {RF}')


            colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:4]

            # RM
            plt.plot(self.RVs_RM[transit_sorted < 1], true_anomaly_sorted[transit_sorted < 1] / (2 * np.pi), c=colors[0], label='RM', lw=3, ls='dashed')  # Use self.transit
            # planet
            plt.plot(self.RVs_planet, true_anomaly_sorted / (2 * np.pi),c=colors[1],  label='planet', lw=4)
            # star
            plt.plot(self.RVs_star, true_anomaly_sorted / (2 * np.pi), c=colors[2], ls='solid', label='star', lw=3)
            # tellurics:
            #plt.plot(self.RVs_tel, true_anomaly_sorted / (2 * np.pi), lw=4, c=colors[3], ls='dotted', label='tellurics')  # Use self.vsys


                
            # Contact markers are only meaningful when both states are sampled.
            in_transit = transit_sorted < 1.0
            changes = np.flatnonzero(np.diff(in_transit.astype(int)))
            for change in changes:
                plt.axhline(true_anomaly_sorted[change] / (2 * np.pi),
                            c='k', lw=1, ls='dashed')
                

            #plt.xlim(-RV_ext, RV_ext)
            plt.ylabel('Phase')
            plt.xlabel('RV [km/s]')
            plt.legend(ncols=2)
            plt.tight_layout()
            

            #plt.show()
            if save:
                plt.savefig(f'{self.dp}/RVTraceEstimator_{RF}.png', dpi=180)