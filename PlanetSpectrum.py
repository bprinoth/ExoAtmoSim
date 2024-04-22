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
from utilities import paramget, BB

### THE BIG DICTIONARY

species_dictionary = { # from FastChem notation to pRT notation
    "Al": "Al",
    "Fe": "Fe",
    "H2S1": "H2S_main_iso",
    "H3N1": "NH3_main_iso",
    "Ti": "Ti",		      
    "B": "B",
    "H1Fe1": "FeH_main_iso",
    "Be": "Be",
    "C1O1": "CO_main_iso",
    "Fe1+": "FeII",
    "H1C1N1_1": "HCN_main_iso",
    "O": "O", 
    "V": "V",
    "V1+": "VII",
    "H2C2": "C2H2_main_iso",
    "Na": "Na",
    "O2": "O2",
    "K": "K",
    "C": "C",
    "H4C1": "CH4_main_iso",
    "C1O2": "CO2_main_iso",
    "O1Ti1": "TiO_48_Exomol_McKemmish",
    "O3": "O3_main_iso",
    "Mg": "Mg",
    "Li": "Li", 
    "Y": "Y",
    "Ca": "Ca",
    "H3P1": "PH3_main_iso",
    "Ca1+": "CaII",
    "Mg1+": "MgII",
    "H2O1": "H2O_main_iso",
    "Si": "Si",
    "Cr": "Cr",
    "N": "N",  
    "O1Si1": "SiO_main_iso"
}

network_dictionary = {
    'O2': 'O2',
    'NH3': "H3N1",
    'HCN':"H1C1N1_1",
    'CH4':"H4C1",
    'CO':"C1O1",
    'H2O':"H2O1",
    'CO2':"C1O2"
} 


# Molecular formulas in the network
molecular_formulas = ['CH3COOOH', 'C4H9O', 'C3H7O', 'NO3', 'CH3COOO', 'C2H5OO', 'C2H4OOH', 'HONO2',
 'C2H5OOH', 'CH3ONO', 'C3H8CO', 'CH3NO2', '1C4H9', '2C4H9', 'C4H10', 'C3H7OH',
 'CH3OO', 'C4H8Y', 'CH3OOH', 'HNO2', 'CH3OCO', 'C2H5CHO', 'C2H6CO', 'C2H5O',
 'CH3NO', '2C2H4OH', 'NO2', '2C3H7', '1C3H7', '1C2H4OH', 'HONO', 'C3H8', 'HCNN',
 'cC2H4O', 'HCNO', 'C2H5OH', 'N2O', 'C2H3CHOZ', 'OOH', 'CH2CHO', 'H2O2', 'CH3CO',
 'NCO', 'CH3O', 'O2', 'CH3CHO', 'HNO', 'C', 'CHCO', 'CO2H', 'HOCN', 'C2H5', 'C2H',
 'CH2OH', 'CH', 'C2H6', 'C2H3', 'CH2CO', 'NNH', 'H2CN', 'CH3OH', 'N4S', 'N2D', 'CN',
 '1CH2', 'HNCO', 'NO', 'O3P', 'O1D', 'C2H4', 'NH', '3CH2', 'HCO', 'C2H2', 'H2CO',
 'NH2', 'CO2', 'OH', 'CH3', 'HCN', 'NH3', 'CH4', 'N2', 'CO', 'H2O', 'H', 'He', 'H2',
 'N2O4', 'N2O3', 'N2H2', 'N2H3', 'N2H4', 'HNNO', 'HNOH', 'HNO3', 'NH2OH', 'H2NO',
 'CNN', 'H2CNO', 'C2N2', 'HCNH', 'HON', 'NCN', 'HCOH', 'HOCHO', 'HOCH2O']

# Molecular weights belonging to the network
molecular_weights = [75.044, 73.113, 59.080, 62.004, 73.052, 59.068, 60.052, 63.012,
 60.052, 61.037, 74.123, 61.042, 57.084, 57.084, 58.123, 61.084, 47.040, 56.106,
 48.041, 47.013, 59.044, 58.079, 59.072, 45.066, 45.041, 61.083, 46.005, 43.064,
 43.064, 45.067, 47.013, 44.097, 27.026, 44.053, 27.010, 46.069, 44.013, 55.041,
 33.006, 42.037, 34.014, 43.030, 42.018, 31.019, 31.999, 44.053, 29.992, 12.011,
 41.028, 45.018, 43.018, 29.055, 24.021, 31.033, 13.019, 30.070, 25.029, 42.037,
 30.021, 27.018, 32.042, 92.094, 30.008, 26.018, 15.023, 43.025, 30.006, 79.977,
 31.997, 28.054, 15.015, 15.034, 29.018, 26.038, 30.026, 16.018, 44.010, 17.007,
 15.035, 27.026, 17.031, 16.043, 28.013, 28.010, 18.015, 1.008, 4.003, 2.016,
 92.011, 76.011, 30.028, 31.021, 32.028, 61.022, 33.018, 63.019, 33.033, 33.018,
 27.022, 41.018, 52.025, 27.022, 17.008, 42.024, 29.029, 45.032, 46.040]

# Create dictionary
molecular_weights_dict = dict(zip(molecular_formulas, molecular_weights))
inv_network_dictionary  = {v: k for k, v in network_dictionary.items()}

class Planet:

    def __init__(self, dp, wl_range=[0.4, 0.8], pressure_limits=[-9, 1]): # wl_range in micron, pressure_limits in 10** bar
        
        self.dp = dp # Data path to where your config file is
        
        self.Rp = paramget('Rp', self.dp)  # Rp in Rjup
        self.Rs = paramget('Rs', self.dp)  # Rs in Rsun
        
        
        try:
            self.gp = 10**paramget('loggp', self.dp) # logp in cgs units
        except:
            self.Mp = paramget('Mp', self.dp) # Mp in Mjup
            self.gp = (const.G * self.Mp * u.Mjup / (self.Rs * u.Rjup)**2).to('cm/s2').value # computes the surface gravity from the mass instead
            
        #self.pRT_species = paramget('pRT', self.dp).split(',')
        #self.FastChem_species = paramget('FastChem', self.dp).split(',')
        self.wl_range = wl_range
        
        # Condensation:
        
        self.cond = paramget('condensation', self.dp) # can be rainout, eq. or False
        
        try:
            self.tp_profile = np.loadtxt(paramget('tp_profile', self.dp))
            
        except:
            self.tp_profile = paramget('tp_profile', self.dp)
            
            if self.tp_profile == 'isothermal':
            
                self.temp = paramget('Teq', self.dp) # temperature in Kelvin
                self.layers = int(paramget('layers', self.dp)) # layers in your atmosphere
                
            elif self.tp_profile == 'guillot':
                
                # try figuring out the values, otherwise we default back
                try: self.kappa_IR = paramget('kappa_IR', self.dp)
                except: self.kappa_IR = 0.01
                
                try: self.gamma = paramget('gamma', self.dp)
                except: self.gamma = 0.4

                try: self.T_int = paramget('T_int', self.dp)
                except: self.T_int = 200 
                
                self.temp = paramget('Teq', self.dp) # eq. temperature in Kelvin
                self.layers = int(paramget('layers', self.dp))
            
            elif self.tp_profile == 'network':
                self.temperature = None
                self.pressure = None    
                
            
        self.chemistry = paramget('chemistry', self.dp)
        
        if self.chemistry == 'FastChem':
            self.FastChem_input_path = paramget('FastChem_input', self.dp)
            
        elif self.chemistry == 'network':
            print('[INFO] Loading chemistry from network.')
            self.file_path_network = paramget('network_path', self.dp)
            
        
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
            
        try: 
            self.Ti_H = paramget('Ti/H', self.dp)
        except: 
            self.Ti_H = 'solar'
        # # clouds?
        # try:
        #     self.clouds = paramget('clouds', self.dp)
            
        # except:
        #     self.clouds = None
        
            
    def compute_chemistry(self):
        
        if self.tp_profile == 'isothermal':
            print(f"[INFO] Assuming an isothermal temperature pressure profile with T={self.temp} K over {self.layers} layers.")
            self.temperature = np.full(self.layers, self.temp)
            self.pressure = np.logspace(-9, 1, num=self.layers)
            
        elif self.tp_profile == 'guillot':
            
            print(f"[INFO] Assuming Guillot temperature pressure profile over {self.layers} layers.")
            print(f"[INFO] kappa_IR = {self.kappa_IR}, gamma = {self.gamma}, T (int) = {self.T_int} K, T (eq) = {self.temp} K")
            
            self.pressure = np.logspace(-9, 1, num=self.layers)     
            self.temperature = guillot_global(self.pressure, self.kappa_IR, self.gamma, self.gp, self.T_int, self.temp)     
            
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
            self.abundances = np.array(self.fastchem.getElementAbundances())
            
            for j in range(0, self.fastchem.getElementNumber()):
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
            
        if self.Ti_H != 'solar':
            print(f'[INFO] Changing your Ti_H ratio to {self.Ti_H}')
            index_Ti = self.fastchem.getGasSpeciesIndex('Ti')
            index_H = self.fastchem.getGasSpeciesIndex('H')
        
            self.abundances[index_Ti] = self.abundances[index_H] * self.Ti_H 
            
        print("[INFO] Running FastChem")
        fastchem_flag = self.fastchem.calcDensities(input_data, output_data)
        
        print("[INFO] FastChem reports:")
        print("[INFO]  ---", pyfastchem.FASTCHEM_MSG[fastchem_flag])

        if np.amin(output_data.element_conserved[:]) == 1:
            print("[INFO]  --- element conservation: ok")
        else:
            print("[INFO]  --- element conservation: fail")
            
            
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
    
    def extract_molecular_weights(self):
        self.molecular_weights_dict = molecular_weights_dict 
   
    def load_network(self):
        
        # Read the dataformat
        comment_rows = 0
        with open(self.file_path_network, 'r') as file:
            for line in file:
                if line.startswith('!'):
                    comment_rows += 1
                else:
                    break  # Exit the loop as soon as non-comment lines are encountered

        print("Number of comment rows:", comment_rows)

        # reading it once to get the comment 
        data = pd.read_csv(self.file_path_network, delim_whitespace=True, skiprows=comment_rows-1, header=None)
        # Assign the first row as the column headers
        header = data.iloc[0][1::]

        data = pd.read_csv(self.file_path_network, delim_whitespace=True, comment='!', header=None)
        data.columns = header
        
        # we need to flip the axis!!!
        self.temperature = np.asarray(data['temperature[K]'])[::-1]
        self.pressure = np.asarray(data['pressure[bar]'])[::-1]
        self.mean_molecular_weight = (10**np.asarray(data['MMW[g]']) * u.g).to('u').value # translation from grams to atomic u.
        
        self.network_data = data
     
    # def plot_abundances(self):
    #     if self.chemistry == 'network':
            
    
     
        
    def compute_single_species_template(self, template_species, save_name, mode='transmission', clouds=None, hazes=None): 
        
        if len(template_species) == 1:
            print(f'[INFO] Computing single {mode} species template for {template_species[0]}.')
            #print(f'[INFO] ')
            self.compute_spectrum(template_species, save_name, mode=mode, clouds=clouds, hazes=hazes, template=True)
            
        else: 
            print(f'[INFO] You should only have one template species for your template. Otherwise consider using the self.compute_spectrum function to produce a spectrum.')
    
    #def compute_continuum(self, savename, mode='transmission', clouds=None, template=False):
        
        
    def compute_spectrum(self, template_species, save_name, mode='transmission', clouds=None, hazes=None, template=False): 
        # can be transmission or emission, if Pcloud is None, we ignore it, otherwise we set the cloud at that pressure.
        template_components = ['H2', 'He'] + template_species

        try:
            if self.chemistry == 'FastChem':
                index = self.fastchem.getGasSpeciesIndex(template_species[0])
                template_components = ['H2', 'He'] + template_species
                template_species_pRT_notation = [species_dictionary[template_species[i]] for i in range(len(template_species))]
                #template_components_pRT_notation = ['H2', 'He'] + template_species_pRT_notation
       
            elif self.chemistry == 'network':
                size = np.shape(self.network_data)
                # first to FastChem notation
                
                template_components_network_notation = ['H2', 'He'] + [inv_network_dictionary[template_species[i]] for i in range(len(template_species))] 
                
                # this translates all the components to the pRT notation, so we can use them in the 
                template_species_pRT_notation = [species_dictionary[template_species[i]] for i in range(len(template_species))]
                
                # all the template components
                template_components = ['H2', 'He'] + template_species
        
        except:
            print(f"[WARNING] I cannot find your chemistry. Building with {self.chemistry} from scratch.")
            
            if self.chemistry == 'FastChem':
                self.compute_chemistry()
                index = self.fastchem.getGasSpeciesIndex(template_species[0])
                template_species_pRT_notation = [species_dictionary[template_species[i]] for i in range(len(template_species))]
                #template_components_pRT_notation = ['H2', 'He'] + template_species_pRT_notation

            elif self.chemistry == 'network':
                self.load_network()
                # this translates all the components to the network notation, so we can access the file
                template_components_network_notation = ['H2', 'He'] + [inv_network_dictionary[template_species[i]] for i in range(len(template_species))] 
                
                # this translates all the components to the pRT notation, so we can use them in the 
                template_species_pRT_notation = [species_dictionary[template_species[i]] for i in range(len(template_species))]
                
                # all the template components
                template_components = ['H2', 'He'] + template_species
                
            else:
                print(f"[ERROR] Something is terribly wrong. Chemistry is set to {self.chemistry}, but I only know FastChem or network." )
                return -1
        
        self.mass_fractions = np.zeros((len(template_components), len(self.pressure)))
        self.indices = []
        
        print(f"[INFO] Computing mass fractions for {template_components} with {self.chemistry} chemistry.")
        if self.chemistry == 'FastChem':
            
            for i, species in enumerate(template_components):
                
                index = self.fastchem.getGasSpeciesIndex(species)
                self.indices.append(index)
                
                VMR = self.number_densities[:, self.indices[i]]/self.gas_number_density
                
                molecular_weight = self.fastchem.getGasSpeciesWeight(index)
                
                self.mass_fractions[i] = (VMR * molecular_weight / self.mean_molecular_weight)
                
        elif self.chemistry == 'network':
            
            for i, species in enumerate(template_components):
                
                
                self.mass_fractions[i] = 10**np.asarray(self.network_data[f'{template_components_network_notation[i]}']) * molecular_weights_dict[f'{template_components_network_notation[i]}'] / self.mean_molecular_weight  #np.asarray(10**data)[:, 4::] * molecular_weights).T  / mean_molecular_weight
        
        #import pdb
        #pdb.set_trace()
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
            print('[INFO] Assuming Rayleigh-like scattering.')
            kappa_zero = 0.01
            gamma_scat = -4
            Pcloud = None
            
        elif clouds == 'weak': # powerlaw weak
            print('[INFO] Assuming weak scattering.')
            kappa_zero = 0.01
            gamma_scat = -2.
            Pcloud = None
            
        elif clouds == 'flat': # flat opacity
            print('[INFO] Assuming a flat opacity.')
            kappa_zero = 0.01
            gamma_scat = 0.
            Pcloud = None
            
        elif clouds == 'positive': # 1, positive
            print('[INFO] Assuming the exotic case of a positiv opacity slope.')
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
            if hazes is not None:
                print(f'[INFO] Assuming hazes with a factor {hazes}.')
                haze_factor = hazes
            
            else:
                print('[INFO] Ignoring hazes.')
                haze_factor = None
        else:
            print('[INFO] Ignoring hazes for emission.')
            
        template_components_pRT = ['H2', 'He'] +  template_species_pRT_notation
        mass_fractions_dict = {}
            
        for i in range(len(self.mass_fractions)):
            mass_fractions_dict[template_components_pRT[i]] = self.mass_fractions[i]
            MMW = self.mean_molecular_weight
        #pdb.set_trace()        
        
        if mode == 'transmission':
            
            atmosphere = pRT.Radtrans(line_species = template_species_pRT_notation,
                        rayleigh_species = ['H2', 'He'],
                        continuum_opacities = ['H2-H2', 'H2-He'],
                        wlen_bords_micron = self.wl_range, mode='lbl')
        
            atmosphere.setup_opa_structure(self.pressure)
        
            
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
            
            atmosphere = pRT.Radtrans(line_species = template_species_pRT_notation,
                        rayleigh_species = ['H2', 'He'],
                        continuum_opacities = ['H2-H2', 'H2-He'],
                        wlen_bords_micron = self.wl_range, mode='lbl',
                        #do_scat_emis = False
                        )
        
            atmosphere.setup_opa_structure(self.pressure)
            
            R_pl = self.Rp * nc.r_jup_mean 
            gravity = self.gp
            P0 = self.p0
            
            print('[INFO] Calculating emission spectrum without scattering in planetary flux units')
            atmosphere.calc_flux(self.temperature, mass_fractions_dict, gravity, MMW, kappa_zero = kappa_zero, gamma_scat = gamma_scat, Pcloud=Pcloud)
            
            self.wl_nm = (nc.c/atmosphere.freq/1e-4 * u.micron).to(u.nm)
            self.flux = atmosphere.flux 
            
            if template==False:
                print(f'[INFO] Saving emission spectrum to {save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in fp/fs (flux density contrast).')
                
                self.Teff = paramget('Teff', self.dp)
                print(f"[INFO] Assuming a blackbody for the star at T_eff = {self.Teff} K.")
                
                # Unit conversion from hell
                prt_unit = u.Unit(u.erg/u.cm**2/u.s/u.Hz)
                target_unit = u.Unit(u.erg / u.cm**3/u.s)
                
                self.flux_frei = (self.flux * prt_unit).to(target_unit, u.spectral_density(self.wl_nm))
                
                Fs = BB(self.Teff * u.K)(self.wl_nm)
                self.FpFs = (self.flux_frei / Fs).decompose().value
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, self.FpFs)), overwrite=True)
                
            else:
                print(f'[INFO] Saving emission template to {save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in Fp/Fs (flux contrast).')
                
                self.Teff = paramget('Teff', self.dp)
                print(f"[INFO] Assuming a blackbody for the star at T_eff = {self.Teff} K.")
                
                # Unit conversion from hell
                prt_unit = u.Unit(u.erg/u.cm**2/u.s/u.Hz)
                target_unit = u.Unit(u.erg / u.cm**3/u.s)
                
                self.flux_frei = (self.flux * prt_unit).to(target_unit, u.spectral_density(self.wl_nm))
                
                Fs = BB(self.Teff * u.K)(self.wl_nm)
                self.FpFs = (self.flux_frei * (self.Rp * u.Rjup)** 2 / (Fs * (self.Rs * u.Rsun)**2)).decompose().value
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, 1.+self.FpFs)), overwrite=True)
                
            
        elif mode == 'emission_scat':
            
            atmosphere = pRT.Radtrans(line_species = template_species_pRT_notation,
                        rayleigh_species = ['H2', 'He'],
                        continuum_opacities = ['H2-H2', 'H2-He'],
                        wlen_bords_micron = self.wl_range, mode='lbl',
                        do_scat_emis = True
                        )
        
            atmosphere.setup_opa_structure(self.pressure)        
            
            R_pl = self.Rp * nc.r_jup_mean 
            gravity = self.gp
            P0 = self.p0
            
            #self.wl_nm = (nc.c/atmosphere.freq/1e-4 * u.micron).to(u.nm)
            #self.flux = atmosphere.flux
            
            print('[INFO] Calculating emission spectrum with scattering in planetary flux units')
            atmosphere.calc_flux(self.temperature, mass_fractions_dict, gravity, MMW, kappa_zero = kappa_zero, gamma_scat = gamma_scat, Pcloud=Pcloud)
        
            
            self.wl_nm = (nc.c/atmosphere.freq/1e-4 * u.micron).to(u.nm)
            self.flux = atmosphere.flux 
            
            if template==False:
                print(f'[INFO] Saving emission spectrum to output/{save_name}.fits')
                print(f'[INFO] The spectrum is in nm, and vacuum, in fp/fs (flux density contrast).')
                
                self.Teff = paramget('Teff', self.dp)
                print(f"[INFO] Assuming a blackbody for the star at T_eff = {self.Teff} K.")
                
                # Unit conversion from hell
                prt_unit = u.Unit(u.erg/u.cm**2/u.s/u.Hz)
                target_unit = u.Unit(u.erg / u.cm**3/u.s)
                
                self.flux_frei = (self.flux * prt_unit).to(target_unit, u.spectral_density(self.wl_nm))
                
                Fs = BB(self.Teff * u.K)(self.wl_nm)
                self.FpFs = (self.flux_frei / Fs).decompose().value
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, self.FpFs)), overwrite=True)
                
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
                
                fits.writeto(save_name+'.fits', np.vstack((self.wl_nm.value, 1.+self.FpFs)), overwrite=True)