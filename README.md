# ExoSim
THIS IS WORK IN PROGRESS. USE AT YOUR OWN RISK ;) 

Proposals to observe exoplanet atmospheres often require the presence of some model observations to show that the expected atmospheric features can indeed be detected using the cross-correlation technique. 
Many assumptions have to be made, but the principle is simple: Do we have a high enough SNR to see the feature we expect? Is the feature even expected?

ExoSim provides a short simulation package to simulate the detection of atmospheric features in exoplanet atmospheres, linking <a href="https://github.com/exoclime/FastChem">FastChem Cond</a> to compute abundance profiles,  <a href="https://gitlab.com/mauricemolli/petitRADTRANS.git">petitRADTRANS</a> to simulate the spectra, and <a href="https://github.com/Hoeijmakers/tayph">tayph</a> to calculate the cross-correlation maps. 

The requirement is that each of these packages is already successfully installed on your machine, together with your opacities for petitRADTRANS. We refer you to the corresponding installation manuals.

## Supported instruments (as of February 2024)

- ESPRESSO: with constant SNR
- SPIRou: with constant SNR
- ANDES: with wavelength-dependent SNR, see ANDES generator
- HARPS: with constant SNR
- HARPS-N: with constant SNR, see HARPS


## Additional files

- Order definition: Defines the wavelength solution of the instruments. For convenience, these are provided with the installation of this module. You need to make sure you add them to your folder before creating your observations.
- Blaze functions: If known, we provide the blaze functions. These are theoretically not necessary an can be just an array of ones.
  
## Generating a simple atmospheric model

This is a bit of trial and error tbh. So don't get too disappointed if you run into problems. This may be on me and not on you. 

To generate a simple atmosphere, you need to make a few choices:

1) **Chemistry:** We currently only support FastChem Cond as the chemistry setting (update for kinetics around the corner.
2) **Condensation:** Using FastChem Cond, you can opt for condensation through rainout, equilibrium condensation or no condensation.
3) **T-p profile:** You can use an isothermal profile, a Guillot profile or a custom profile.

These choices are all taken in the config file, in this example located in a folder called 'testing_the_module'. See config_example for an example file.

Let us first read in the paths and generate an object for the planet spectrum.

```python

import sys

# this is just me adding my petitRADTRANS & FastChem to the path

sys.path.insert(0,'/home/bibi/GitHub/petitRADTRANS/') 
sys.path.insert(0, '/home/bibi/.local/lib/python3.8/site-packages/lib/python3.8/site-packages/pyfastchem-3.0-py3.8-linux-x86_64.egg/')

# and this is just me setting the pRT path.
import os
os.environ['pRT_input_data_path'] = '/usr/local/src/petitRADTRANS/petitRADTRANS/input_data/'
os.environ["OMP_NUM_THREADS"] = "1"

import PlanetSpectrum # this is where all the magic happens.

# Reading in the config file located at testing_the_module/ and initiating the Planet object with a spectrum over the wavelength range from 0.3 to 0.8 micron
 W77Ab = PlanetSpectrum.Planet(
    dp='testing_the_module/', 
    wl_range=[0.3, 0.8],
    pressure_limits=[-9, 1] # from 10^-9 to 1 bar in log space
)

```
Now we can compute the chemistry according to your settings using:

```python

# Computing the chemistry according to your settings. This computes the chemistry for all the molecules and atoms in FastChem Cond (see their documentation).
W77Ab.compute_chemistry()

```

Now we need to define the species in our atmosphere. 
Note that there is a dictionary at the beginning of PlanetSpectrum.py that translates between the FastChem Cond (Hill notation) and pRT notation, see dic.txt file for current dict.
These may be subject to change, so make sure you check those out before you run it!

```python

template_species = [
    'O1Ti1',
    #'O1V1',
    'H2O1',
    'C1O1',
    'Fe',
    'Ti'
                   ]
```
This now produces the spectrum with your species and saves it. Make sure you check out the different modes. You can also add clouds and hazes following the pRT definitions. Check out their documentation.
For clouds, you need a clouds keyword in your config file, as well as the scattering coefficients following pRT definition. The following implementations exist:

- clouds == 'Rayleigh': Rayleigh-like scattering
- clouds == 'weak': power law weak (kappa_0 = 0.01, gamma = -2)
- clouds == 'flat': flat opacity (kappa_0 = 0.01, gamma = 0)
- clouds == 'positive': exotic case of a positive opacity slope (kappa_0 = 0.01, gamma = 1.)
- clouds == float: This means we want a cloud layer at a given pressure level in bar.
- hazes == float: haze factor according to pRT

```python
W77Ab.compute_spectrum(
   template_species=template_species, # including those species
   save_name='W77Ab', #saves here .fits
   mode='emission_no_scat' # could also be emission_scat, transmission
)
```

If you now want to produce the template for your cc, you can run for example:

```python
W77Ab.compute_single_species_template(
    template_species=['O1Ti1'], 
    save_name='TiO_template',
    mode='emission_no_scat'
)
```


This has now produced an emission spectrum in fp/fs (flux density contrast) assuming an underlying black body for the star. 

## So how do I predict my cross-correlation signal now?

Hold on a second. We first need to produce our simulated data!
This step is easily done if you have a bit of an idea how tayph works. If not, go check out the documentation. You'll need a standard tayph runfile (for us 'W77.dat', make sure to turn off berv correction!!) with one added line:

<code> mockdatapath      datapath/to/your/fake/datafolder</code>

Simple as that. But now you need to produce the simulated data. In your data folder, you want to make sure you have the order definition file, and if needed, the blaze function file. Additionally, you need the tayph config file. 

```python
import sys
sys.path.insert(0, '/home/bibi/GitHub/tayph/') 
import tayph.mock as mock
import tayph.run as run

savename = 'W77Ab' # this is the savename you had for your fits file.
mock_dp = f'data/WASP-77Ab/' # this is the datapath to your fake data

```

To bring the spectra into the right format, you need:

```python
mock.translate_phi_spectra(f'{savename}.fits', mock_dp, savename, phi0=0.35)
```

This assigns this spectrum to phase phi0=0.35. This is important if you e.g. use emission spectroscopy and want to centre around phi0. 
For transmission spectra, this becomes phi0 = 0.0

The mock observations are then generated through:

```python
mock.create_mock_obs(f'W77.dat', create_obs_file=True, mode='ESPRESSO', real_data=False, spec='flux', rot=True)
```

For the cross-correlation you can then use, see tayph for full documentation. You can of course use these orders also for your own cc code. No need to use tayph for it. Be aware though that the mock tool is part of tayph.

```python
run.start_run(f'W77.dat', xcor_parallel=False)
```


## Eccentric orbits

We all know this: Eccentric orbits are special. Everything above works like a charm for circular orbits, but as soon as we introduce eccentricity, it all becomes significantly more difficult. The simulation tool allows for creating mock data with eccentric orbits, but tayph is not equipped (yet) with the same functionality. So you need to turn off the Keplerian correction that corrects the orbital motion and do that outside, moving your spectra to the rest-frame of the star. This then produces your ccf (cleaning steps). From there you need to recalculate the out and in transit exposures (if that is what you do!) because tayph intrinsically assumes a circular orbit again. 

When generating your exposures for the eccentric orbit, you need to be a bit careful because 'orbital phase' is not a thing anymore. Much more, you need the true anomaly. So you want to figure out first at which phases your planet is where you think it should be. You can use the RVEstimator for this (rv.bibianaprinoth.ch for an interactive tool). Then you produce your phases and tell the code which phases to consider.

```python
mock.create_mock_obs(f'TOI-3362-mock-ANDES.dat', create_obs_file=True, mode='ANDES', real_data=False, spec='depth', rot=True, eccentric_orbit=True, phase_range=phases) # this is an example for a transiting planet. Phases are not the true anomaly.
```
