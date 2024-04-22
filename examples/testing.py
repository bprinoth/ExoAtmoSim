import sys
import os

sys.path.insert(0,'/home/bibi/GitHub/petitRADTRANS/')
sys.path.insert(0, '/home/bibi/.local/lib/python3.8/site-packages/lib/python3.8/site-packages/pyfastchem-3.0-py3.8-linux-x86_64.egg/')

os.environ['pRT_input_data_path'] = '/usr/local/src/petitRADTRANS/petitRADTRANS/input_data/'
os.environ["OMP_NUM_THREADS"] = "1"

#sys.path.insert(0, '/home/bibi/GitHub/ExoSim/')

import PlanetSpectrum

## PLOTTING
import matplotlib.pyplot as plt
#plt.rcParams.update({'font.size': 16})
plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["ytick.major.width"] = 2
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5
plt.rcParams["xtick.minor.size"] = 3.5
plt.rcParams["ytick.minor.size"] = 3.5
plt.rcParams["ytick.minor.width"] = 1
plt.rcParams["xtick.minor.width"] = 1
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['legend.facecolor'] = 'white'
plt.rcParams['legend.edgecolor'] = 'white'
plt.rcParams['legend.framealpha'] = 1


W77Ab = PlanetSpectrum.Planet(
    dp='testing_the_module/', 
    wl_range=[0.3, 0.8], 
)

W77Ab.compute_chemistry()

template_species = [
    'O1Ti1',
    #'O1V1',
    'H2O1',
    'C1O1',
    'Fe',
    'Ti'
                   ]

W77Ab.compute_spectrum(
   template_species=template_species, 
   save_name='W77Ab',
   mode='emission_no_scat'
)


# W77Ab.compute_single_species_template(
#     template_species=['O1Ti1'], 
#     save_name='TiO_template',
#     mode='emission_no_scat'
# )

# W77Ab.compute_single_species_template(
#     template_species=['Fe'], 
#     save_name='Fe_template',
#     mode='emission_no_scat'
# )
