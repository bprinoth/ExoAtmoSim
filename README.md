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

```python
 W77Ab = PlanetSpectrum.Planet(
    dp='testing_the_module/', 
    wl_range=[0.3, 0.8], 
)

```
```python

W77Ab.compute_chemistry()

template_species = [
    'O1Ti1',
    #'O1V1',
    'H2O1',
    'C1O1',
    'Fe',
    'Ti'
                   ]

```
