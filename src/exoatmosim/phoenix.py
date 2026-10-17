"""PHOENIX HiRes stellar FITS retrieval using Astropy's downloadable cache.

The wavelength file and model fluxes originate from the PHOENIX archive.
Returned wavelength units are Angstrom, matching the native FITS files.
"""
import numpy as np

phoenix_model_temps = np.array([3200,6400,14500,4100,12500,5000,6700,10400,2700,11600,3600,7000,7600,4500,10600,5400,5700,3100,7400,10000,4000,6000,4900,8400,7200,2600,6300,15000,3500,4400,6600,9600,13500,5300,7800,3000,6900,9400,3900,12000,4800,5600,2500,9200,8800,3400,5900,4300,8000,9000,5200,6200,9800,2900,13000,3800,6500,11400,14000,4700,2400,6800,8600,3300,4200,5500,11800,5100,11000,10200,2800,5800,3700,8200,10800,4600,6100,2300,11200])
phoenix_model_loggs = np.arange(0., 6., 0.5)
_ROOT = 'ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'

def get_url(T_eff, log_g):
    t = int(phoenix_model_temps[np.argmin(abs(phoenix_model_temps-T_eff))])
    g = float(phoenix_model_loggs[np.argmin(abs(phoenix_model_loggs-log_g))])
    return (_ROOT + 'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/' +
            f'lte{t:05d}-{g:1.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

def _fetch(url, cache=True):
    from astropy.utils.data import download_file
    from astropy.io import fits
    return np.asarray(fits.getdata(download_file(url, cache=cache, timeout=30)))

def get_phoenix_model_spectrum(T_eff, log_g=4.5, cache=True):
    return _fetch(get_url(T_eff, log_g), cache=cache)

def get_phoenix_wavelengths(cache=True, vacuum=True):
    wavelengths_vacuum = _fetch(_ROOT + 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits', cache=cache)
    if vacuum:
        return wavelengths_vacuum
    # Native PHOENIX FITS wavelengths are Angstrom; conversion function uses nm.
    from .wavelengths import vactoair
    return vactoair(wavelengths_vacuum / 10.) * 10.
