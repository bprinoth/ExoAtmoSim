# ExoSim

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
  
