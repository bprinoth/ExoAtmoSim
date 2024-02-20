# ExoSim

Proposals to observe exoplanet atmospheres often require the presence of some model observations to show that the expected atmospheric features can indeed be detected using the cross-correlation technique. 
Many assumptions have to be made, but the principle is simple: Do we have a high enough SNR to see the feature we expect? Is the feature even expected?

ExoSim provides a short simulation package to simulate the detection of atmospheric features in exoplanet atmospheres, linking FastChem Cond (https://github.com/exoclime/FastChem) to compute abundance profiles, 
petitRADTRANS (https://gitlab.com/mauricemolli/petitRADTRANS.git) to compute the spectra, and tayph (https://github.com/Hoeijmakers/tayph) to calculate the cross-correlation maps. 


## Supported instruments (as of February 2024)


## Additional files

- Order definition: Defines the wavelength solution of the instruments. For convenience, these are provided with the installation of this module. You need to make sure you add them to your folder before creating your observations.
- Blaze functions: If known, we provide the blaze functions. It is theoretically not necessary, but makes your data look more 'realistic'.
  
