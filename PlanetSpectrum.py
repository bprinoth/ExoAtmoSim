import pyfastchem # requirement that FastChemCond is installed
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
import astropy.units as u
import pdb
import seaborn as sns
import pandas as pd
from astropy.io import fits
import petitRADTRANS as pRT
#from petitRADTRANS.radtrans import Radtrans
from petitRADTRANS import nat_cst as nc
from petitRADTRANS.physics import guillot_global

######### HELPER FUNCTIONS - DEFINED OUTSIDE
# # paramget allows you to get parameters from a config file instead of typing all the stuff into a list

def paramget(keyword,dp,full_path=False,force_string = False):
    """This code queries a planet system parameter from a config file located in the folder
    specified by the path dp; or run configuration parameters from a file speciefied by the full
    path dp, if full_path is set to True. It is taken from tayph @Hoeijmakers

    Parameters
    ----------
    keyword : str
        A keyword present in the cofig file.

    dp : str, Path
        Output filename/path.

    full_path: bool
        If set, dp refers to the actual file, not the location of a folder with a config.dat;
        but the actual file itself.

    Returns
    -------
    value : int, float, bool, str
        The value corresponding to the requested keyword.

    """
    import pathlib
    import distutils.util
    import pdb


    #This is probably the only case where I don't need obs_times and config to exist together...
    #dp=check_path(dp)
    #typetest(keyword,str,'keyword in paramget()')

    if isinstance(dp,str) == True:
        dp=pathlib.Path(dp)
    try:
        if full_path == False:
            dp = dp/'config'
        f = open(dp, 'r')
    except FileNotFoundError:
        raise FileNotFoundError('Parameter file does not exist at %s' % str(dp)) from None
    
    x = f.read().splitlines()
    f.close()
    n_lines=len(x)
    keywords={}
    
    for i in range(0,n_lines):
        line=x[i].split()
        if len(line) > 1:
            if force_string:
                value=(line[1])
            else:
                try:
                    value=float(line[1])
                except ValueError:
                    try:
                        value=bool(distutils.util.strtobool(line[1]))
                    except ValueError:
                        value=(line[1])
            keywords[line[0]] = value
    try:
        return(keywords[keyword])
    
    except KeyError:
        # print(keywords)
        raise Exception('Keyword %s is not present in parameter file at %s' % (keyword,dp)) from None

def BB(temperature): # BB function from frei @ Brett Morris.
    import numpy as np
    
    """
    Compute the blackbody flux

    Parameters
    ----------
    temperature : ~astropy.units.Quantity
        Temperature of the blackbody

    Returns
    -------
    bb : function
        Returns a function which takes the wavelength as an
        `~astropy.units.Quantity` and returns the Planck flux
    """
    # h = 6.62607015e-34  # J s
    # c = 299792458.0  # m/s
    # k_B = 1.380649e-23  # J/K
    return lambda wavelength: (
            2 * const.h * const.c ** 2 / np.power(wavelength, 5) /
            np.expm1(const.h * const.c / (wavelength * const.k_B * temperature))
    )
    

class Planet:

    def __init__(self, dp, wl_range=[0.4, 0.8], pressure_limits=[-9, 1]): # wl_range in micron, pressure_limits in 10** bar
        
        self.dp = dp # Data path to where your config file is
        
        self.Rp = paramget('Rp', self.dp)  # Rp in Rjup
        self.Rs = paramget('Rs', self.dp)  # Rs in Rsun
        
        
        try:
            self.gp = 10**paramget('logp', self.dp) # logp in cgs units
        except:
            self.Mp = paramget('Mp', self.dp) # Mp in Mjup
            self.gp = (const.G * self.Mp * u.Mjup / (self.Rs * u.Rjup)**2).to('cm/s2') # computes the surface gravity from the mass instead
            
        self.pRT_species = paramget('pRT', self.dp).split(',')
        self.FastChem_species = paramget('FastChem', self.dp).split(',')
        self.wl_range = wl_range
        
        # Condensation:
        
        self.cond = paramget('condensation', self.dp) # can be rainout, eq. or False
        
        try:
            self.tp_profile = np.loadtxt(paramget('tp_profile', self.dp))
            
        except:
            self.tp_profile = paramget('tp_profile', self.dp)
            
            if self.tp_profile == 'isothermal':
            
                self.temp = paramget('T', self.dp) # temperature in Kelvin
                self.layers = paramget('layers', self.dp) # layers in your atmosphere
                
            elif self.tp_profile == 'guillot':
                
                # try figuring out the values, otherwise we default back
                try: self.kappa_IR = paramget('kappa_IR', self.dp)
                except: self.kappa_IR = 0.01
                
                try: self.gamma = paramget('gamma', self.dp)
                except: self.gamma = 0.4

                try: self.T_int = paramget('T_int', self.dp)
                except: self.T_int = 200 
                
                self.temp = paramget('T', self.dp) # eq. temperature in Kelvin
                
            
            
            
        self.FastChem_input_path = paramget('FastChem_input', self.dp)
        self.p0 = paramget('p0', self.dp) # reads ref pressure in bar
        
        # some other optional inputs:
        
        # metallicity
        try: 
            self.metallicity = paramget('metallicity', self.dp)
        except: 
            self.metallicity = 0.0
            
        # C/O ratio, because people care about this it seems
        
        try: 
            self.C_O = paramget('C/O', self.dp)
        except: 
            self.C_O = 'solar'
        
        # clouds?
        #
        
            
    def compute_chemistry(self):
        
        if self.tp_profile == 'isothermal':
            print(f"[INFO] Assuming an isothermal temperature pressure profile with T={self.temp} K over {self.layers} layers.")
            self.temperature = np.full(self.layers, self.temp)
            self.pressure = np.logspace(-9, 1, num=self.layers)
            
        elif self.tp_profile == 'guillot':
            
            print(f"[INFO] Assuming Guillot temperature pressure profile over {self.layers} layers.")
            print(f"[INFO] kappa_IR = {self.kappa_IR}, gamma = {self.gamma}, T (int) = {self.T_int} K, T (eq) = {self.temp} K")
            
            self.pressure = np.logspace(-9, 1, num=self.layers)     
            self.temperature = guillot_global(self.pressure, self.kappa_IR, self.gamma, self.gp, self.T_int, self.T_equ)     
            
        else:  
            try:
                print(f"[INFO] Assuming your provided temperature pressure profile")
                self.temperature = self.tp_profile[:,0]
                self.pressure = self.tp_profile[:,1] # This assumes you have a txt file with your temperatures in col 0 in K and pressures in col 1 in bar
                
            except: 
                print( "[ERROR] No temperature pressure profile found. Make sure it is in the right format.")
                print( "[ERROR] --- In case you want to use a built-in profile, pick isothermal or guillot.")
                print(f"[ERROR] --- Your current setting is {self.tp_profile}.")
                
        
        # This generates the fastchem input for your run
        input_data = pyfastchem.FastChemInput()
        output_data = pyfastchem.FastChemOutput()
        
        input_data.temperature = self.temperature
        input_data.pressure = self.pressure
        
        if self.cond == False:
            print(f"[INFO] Condensation is set to {self.cond}.")
            print(f"[INFO] Initialising FastChem object without condensation.")
            
            self.fastchem = pyfastchem.FastChem(f'{self.FastChem_input_path}/element_abundances/asplund_2009.dat',
                               f'{self.FastChem_input_path}logK/logK.dat', 0)
            
        elif self.cond == 'rainout':
            print(f"[INFO] Condensation is set to {self.cond}.")
            print(f"[INFO] Initialising FastChem object with rainout condensation.")
            
            self.fastchem = pyfastchem.FastChem(f'{self.FastChem_input_path}/element_abundances/asplund_2009.dat',
                               f'{self.FastChem_input_path}logK/logK.dat', 
                               f'{self.FastChem_input_path}logK/logK_condensates.dat', 1)
            input_data.rainout_condensation = True
            self.fastchem.setParameter('accuracyChem', 1e-5)
            
        elif self.cond == 'equilibrium':
            print(f"[INFO] Condensation is set to {self.cond}.")
            print(f"[INFO] Initialising FastChem object with equilibrium condensation.")
            
            self.fastchem = pyfastchem.FastChem(f'{self.FastChem_input_path}/element_abundances/asplund_2009.dat',
                               f'{self.FastChem_input_path}logK/logK.dat', 
                               f'{self.FastChem_input_path}logK/logK_condensates.dat', 1)
            input_data.equilibrium_condensation = True
            self.fastchem.setParameter('accuracyChem', 1e-5)

        else:
            print(f"[INFO] Condensation is set to {self.cond}, which is not supported by FastChem. Assuming no condensation instead")
            print(f"[INFO] Initialising FastChem object without condensation.")
            
            self.fastchem = pyfastchem.FastChem(f'{self.FastChem_input_path}/element_abundances/asplund_2009.dat',
                               f'{self.FastChem_input_path}logK/logK.dat', 0)
        
        
        
        # If the metallicity is not solar, we are going to scale the abundances
        
        if self.metallicity != 0.0:
            print(f"[INFO] Metallicity is {self.metallicity}. Scaling your abundances.")
            self.abundances = np.array(fastchem.getElementAbundances())
            
            for j in range(0, fastchem.getElementNumber()):
                if self.fastchem.getElementSymbol(j) != 'H' and self.fastchem.getElementSymbol(j) != 'He':
                    self.abundances[j] *= self.metallicity # all except He and H are scaled by the metallicity. 
        
        
        else: 
            print(f"[INFO] Assuming solar metallicity.")
            self.abundances = np.array(self.fastchem.getElementAbundances())
         
         
         
        if self.C_O != 'solar':
            print(f'[INFO] Changing your C/O ratio to {self.C_O}')
            index_C = self.fastchem.getGasSpeciesIndex('C')
            index_O = self.fastchem.getGasSpeciesIndex('O')
        
            self.abundances[index_C] = self.abundances[index_O] * self.C_O 
            
            
        print("[INFO] Running FastChem")
        fastchem_flag = self.fastchem.calcDensities(input_data, output_data)
        
        print("[INFO] FastChem reports:")
        print("  -", pyfastchem.FASTCHEM_MSG[fastchem_flag])

        if np.amin(output_data.element_conserved[:]) == 1:
            print("  --- element conservation: ok")
        else:
            print("  --- element conservation: fail")
            
            
        #convert the output into a numpy array
        self.number_densities = np.array(output_data.number_densities)
        
        if self.cond != False:
            self.number_densities_cond = np.array(output_data.number_densities_cond)

        #total gas particle number density from the ideal gas law
        #used later to convert the number densities to mixing ratios
        
        self.gas_number_density = self.pressure*1e6 / (const.k_B.cgs * self.temperature)

        self.total_element_density = np.array(output_data.total_element_density)

        # extract mean molecular weight
        self.mean_molecular_weight = np.array(output_data.mean_molecular_weight)
        
        return 0
        
    def compute_transmission_spectrum(self, spectrum_species):
        
        try:
            index = self.fastchem.getSpeciesIndex(spectrum_species[0])
        
        except:
            print("[WARNING] I cannot find your chemistry. Computing chemistry from scratch.")
            self.compute_chemistry()
            index = self.fastchem.getSpeciesIndex(spectrum_species[0])
        
        return 0
    
    def compute_single_species_template(self, template_species, save_name, mode='transmission', clouds=None, hazes=None): 
        
        if len(template_species) == 1:
            print(f'[INFO] Computing single {mode} species template for {template_species[0]}.')
            #print(f'[INFO] ')
            self.compute_spectrum(template_species, save_name, mode='transmission', clouds=None, hazes=None, template=True)
            
        else: 
            print(f'[INFO] You should only have one template species for your template. Otherwise consider using the self.compute_spectrum function to produce a spectrum.')
    
        
    def compute_spectrum(self, template_species, save_name, mode='transmission', clouds=None, hazes=None, template=False): 
        # can be transmission or emission, if Pcloud is None, we ignore it, otherwise we set the cloud at that pressure.
        

        try:
            index = self.fastchem.getSpeciesIndex(template_species)
        
        except:
            print("[WARNING] I cannot find your chemistry. Computing chemistry from scratch.")
            self.compute_chemistry()
            index = self.fastchem.getSpeciesIndex(template_species)
            
        
        template_components = ['H2', 'He'] + template_species
        self.mass_fractions = np.zeros((len(template_components), len(self.pressure)))
        self.indices = []
        
        print(f"[INFO] Computing mass fractions for {template_components}")
        for i, species enumerate(template_components):
            
            index = self.fastchem.getSpeciesIndex(species)
            self.indices.append(index)
            
            VMR = self.number_densities[:, self.indices[i]]/self.gas_number_density
            
            molecular_weight = self.fastchem.getSpeciesMolecularWeight(index)
            
            self.mass_fraction[i] = (VMR * molecular_weight / self.mean_molecular_weight)
        
        
        print(f'[INFO] Setting up radiative transfer object with pRT')
        print(f'[INFO] ----- Wavelength coverage from {self.wl_range[0]} to {self.wl_range[-1]} micron')
        print(f'[INFO] ----- Assuming the following parameters:')
        print(f'[INFO] ----- Rp = {self.Rp} Jupiter radii')
        print(f'[INFO] ----- g = {self.gp} cm/ s^2')
        print(f'[INFO] ----- Rs = {self.Rs} solar radii')
        print(f'[INFO] ----- p0 = {self.p0} bar')
        
        if clouds is None:
            print("[INFO] Ignoring clouds.")
            kappa_zero = None
            gamma_scat = None 
            Pcloud = None
            
        elif clouds == 'Rayleigh':
            print('[INFO] Assuming Rayleigh-like scattering')
            kappa_zero = 0.01
            gamma_scat = -4
            Pcloud = None
            
        elif clouds == 'weak': # powerlaw weak
            print('[INFO] Assuming weak scattering')
            kappa_zero = 0.01
            gamma_scat = -2.
            Pcloud = None
            
        elif clouds == 'flat': # flat opacity
            print('[INFO] Assuming a flat opacity')
            kappa_zero = 0.01
            gamma_scat = 0.
            Pcloud = None
            
        elif clouds == 'positiv': # 1, positive
            print('[INFO] Assuming the exotic case of a positiv opacity slope')
            kappa_zero = 0.01
            gamma_scat = 1. 
            Pcloud = None
            
        elif type(clouds) == float: # if it is set to a bar value
            print(f'[INFO] Assuming a cloud deck at P = {clouds} bar.')
            kappa_zero = None
            gamma_scat = None
            Pcloud = clouds
        
        else:
            print("[INFO] Ignoring clouds.")
            kappa_zero = None
            gamma_scat = None 
            Pcloud = None
          
        if mode == 'transmission':  
            if hazes in not None:
                print(f'[INFO] Assuming hazes with a factor {hazes}')
                haze_factor = hazes
            
            else:
                print('[INFO] Ignoring hazes')
                haze_factor = None
        else:
            print('[INFO] Ignoring hazes for emission.')
        
        if mode == 'transmission':
            
            atmosphere = pRT.Radtrans(line_species = [template_species],
                        rayleigh_species = ['H2', 'He'],
                        continuum_opacities = ['H2-H2', 'H2-He'],
                        wlen_bords_micron = self.wl_range, mode='lbl')
        
            atmosphere.setup_opa_structure(self.pressure)

            mass_fractions_dict = {}
            
            for i in range(len(self.mass_fraction)):
                mass_fractions_dict[template_components[i]] = self.mass_fraction[i]
                MMW = self.mean_molecular_weight
        
            
            R_pl = self.Rp * nc.r_jup_mean #weird pRT components
            gravity = self.gp
            P0 = self.p0
    
            print('[INFO] Calculating transmission spectrum in planetary radii')
            atmosphere.calc_transm(self.temperature, mass_fractions_dict, gravity, MMW, R_pl=R_pl, P0_bar=P0, kappa_zero = kappa_zero, gamma_scat = gamma_scat, haze_factor=haze_factor, Pcloud=Pcloud)

            self.wl_nm = (nc.c/atmosphere.freq/1e-4 * u.micron).to(u.nm)
            self.transit_radius = atmosphere.transm_rad/nc.r_jup_mean
            
            print(f'[INFO] Calculating transit depth')
            self.transit_depth = ((self.transit_radius * u.Rjup)**2 / (self.Rs * u.Rsun)**2).decompose()

            if template:
                print(f'[INFO] Saving transmission template to output/{save_name}.fits')
                print(f'[INFO] The template is in nm, and vacuum, in transit radii.')
                
            else: 
                print(f'[INFO] Saving transmission spectrum to output/{save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in transit radii.')
                
            
            fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, 1.- self.transit_depth.value)), overwrite=True)
            # in nm, in vac
            
        elif mode == 'emission_no_scat':
            
            atmosphere = pRT.Radtrans(line_species = [template_species],
                        rayleigh_species = ['H2', 'He'],
                        continuum_opacities = ['H2-H2', 'H2-He'],
                        wlen_bords_micron = self.wl_range, mode='lbl',
                        do_scat_emis = False
                        )
        
            atmosphere.setup_opa_structure(self.pressure)

            mass_fractions_dict = {}
            
            for i in range(len(self.mass_fraction)):
                mass_fractions_dict[template_components[i]] = self.mass_fraction[i]
                MMW = self.mean_molecular_weight
        
            
            R_pl = self.Rp * nc.r_jup_mean 
            gravity = self.gp
            P0 = self.p0
            
            print('[INFO] Calculating emission spectrum without scattering in planetary flux units')
            atmosphere.calc_flux(self.temperature, mass_fractions_dict, gravity, MMW, kappa_zero = kappa_zero, gamma_scat = gamma_scat, Pcloud=Pcloud)
            
            self.wl_nm = (nc.c/atmosphere.freq/1e-4 * u.micron).to(u.nm)
            self.flux = atmosphere.flux 
            
            if template==False:
                print(f'[INFO] Saving emission spectrum to output/{save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in planetary flux [erg / cm2 / s / Hz].')
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, self.flux)), overwrite=True)
                
            else:
                print(f'[INFO] Saving emission template to output/{save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in Fp/Fs (flux contrast).')
                
                self.Teff = paramget('Teff', self.dp)
                print(f"[INFO] Assuming a blackbody for the star at T_eff = {self.Teff} K.")
                
                # Unit conversion from hell
                prt_unit = u.Unit('erg/cm2/s/Hz')
                target_unit = u.Unit('erg/cm3/s')
                
                self.flux_frei = (self.flux * prt_unit).to(target_unit, u.spectral_density(self.wl_nm))
                
                Fs = BB(self.Teff * u.K)(self.wl_nm)
                self.FpFs = (self.flux_frei * (self.Rp * u.Rjup)** 2 / (Fs * (self.Rs * u.Rsun)**2)).decompose()
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, self.FpFs)), overwrite=True)
                
            
        elif mode == 'emission_scat':
            
            atmosphere = pRT.Radtrans(line_species = [template_species],
                        rayleigh_species = ['H2', 'He'],
                        continuum_opacities = ['H2-H2', 'H2-He'],
                        wlen_bords_micron = self.wl_range, mode='lbl',
                        do_scat_emis = True
                        )
        
            atmosphere.setup_opa_structure(self.pressure)

            mass_fractions_dict = {}
            
            for i in range(len(self.mass_fraction)):
                mass_fractions_dict[template_components[i]] = self.mass_fraction[i]
                MMW = self.mean_molecular_weight
        
            
            R_pl = self.Rp * nc.r_jup_mean 
            gravity = self.gp
            P0 = self.p0
            
            self.wl_nm = (nc.c/atmosphere.freq/1e-4 * u.micron).to(u.nm)
            self.flux = atmosphere.flux
            
            print('[INFO] Calculating emission spectrum with scattering in planetary flux units')
            atmosphere.calc_flux(self.temperature, mass_fractions_dict, gravity, MMW, kappa_zero = kappa_zero, gamma_scat = gamma_scat, Pcloud=Pcloud)
        
            
            self.wl_nm = (nc.c/atmosphere.freq/1e-4 * u.micron).to(u.nm)
            self.flux = atmosphere.flux 
            
            if template==False:
                print(f'[INFO] Saving emission spectrum to output/{save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in planetary flux [erg / cm2 / s / Hz].')
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, self.flux)), overwrite=True)
                
            else:
                print(f'[INFO] Saving emission template to output/{save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in Fp/Fs (flux contrast).')
                
                self.Teff = paramget('Teff', self.dp)
                print(f"[INFO] Assuming a blackbody for the star at T_eff = {self.Teff} K.")
                
                # Unit conversion from hell
                prt_unit = u.Unit('erg/cm2/s/Hz')
                target_unit = u.Unit('erg/cm3/s')
                
                self.flux_frei = (self.flux * prt_unit).to(target_unit, u.spectral_density(self.wl_nm))
                
                Fs = BB(self.Teff * u.K)(self.wl_nm)
                self.FpFs = (self.flux_frei * (self.Rp * u.Rjup)** 2 / (Fs * (self.Rs * u.Rsun)**2)).decompose()
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, self.FpFs)), overwrite=True)