"""Observer-frame RVs and physically specified rotational/LSF broadening.

All input wavelength arrays here are AIR nanometres. This is independent of
``previous external mock implementation``; PyAstronomy remains an optional scientific dependency for
Gray-profile rotational broadening, not a external mock dependency.
"""
from __future__ import annotations
import numpy as np

C_KMS = 299792.458


def observer_frame_velocities(star_orbit_kms, planet_orbit_kms, systemic_kms, berv_kms):
    """Match frozen previous implementation's observer frame: orbital + systemic - BERV.

    RVTraceEstimator orbital velocities have zero systemic offset; its BERV is
    the Astropy barycentric *correction*, to subtract from stellar motion.
    """
    star = np.asarray(star_orbit_kms, dtype=float)
    planet = np.asarray(planet_orbit_kms, dtype=float)
    berv = np.asarray(berv_kms, dtype=float)
    systemic = float(systemic_kms)
    if star.ndim != 1 or planet.shape != star.shape or berv.shape != star.shape:
        raise ValueError('Orbital RV and BERV must be equal-length 1-D arrays')
    if not np.isfinite(np.r_[star, planet, berv, systemic]).all():
        raise ValueError('Non-finite orbital RV, systemic RV, or BERV')
    return star + systemic - berv, planet + systemic - berv


def broaden_spectrum(wavelength_air_nm, flux, *, resolution=None, vsini_kms=0.0,
                     limb_darkening=0.0, oversampling=3.0):
    """Apply Gray rotation (star) then Gaussian instrumental line spread.

    The Gaussian's FWHM is c/R *in velocity*, NOT its sigma. Internally
    convolution uses a uniformly sampled log-wavelength grid, where the
    velocity step is constant. Planet spectra use vsini_kms=0. No flux is
    extrapolated onto wider wavelength coverage. Assumes flux is per unit
    wavelength; over small resolving elements this approximation is safe.
    """
    wave = np.asarray(wavelength_air_nm, float)
    f = np.asarray(flux, float)
    if wave.ndim != 1 or f.shape != wave.shape or len(wave) < 8 or not np.isfinite(wave).all() or not np.isfinite(f).all() or np.any(wave <= 0) or np.any(np.diff(wave) <= 0):
        raise ValueError('Expected increasing finite air-nm wavelengths and flux')
    if resolution is not None and (not np.isfinite(resolution) or resolution <= 0):
        raise ValueError('Resolving power must be positive')
    if not np.isfinite(vsini_kms) or vsini_kms < 0 or not 0 <= limb_darkening <= 1:
        raise ValueError('Invalid vsini or limb darkening')
    if not np.isfinite(oversampling) or oversampling < 2:
        raise ValueError('oversampling must be >= 2')
    if vsini_kms == 0 and resolution is None:
        return wave.copy(), f.copy()
    # Reproduce legacy previous implementation's rotational broadening on a linear Angstrom grid.
    if vsini_kms > 0:
        from PyAstronomy import pyasl
        w_linear = np.linspace(wave[0], wave[-1], 2 * len(wave))
        f_linear = np.interp(w_linear, wave, f)
        f_rot = pyasl.fastRotBroad(w_linear * 10., f_linear, limb_darkening, vsini_kms)
        wave, f = w_linear, np.asarray(f_rot, float)
    if resolution is None:
        return wave, f
    from scipy.ndimage import gaussian_filter1d
    # Match at least 3 samples per instrumental FWHM and avoid downsampling
    # a finer input wavelength grid before convolution.
    dv_input = np.min(np.diff(np.log(wave))) * C_KMS
    dv_grid = min(C_KMS / resolution / oversampling, dv_input)
    n = int(np.ceil((np.log(wave[-1]) - np.log(wave[0])) * C_KMS / dv_grid)) + 1
    if n > 5_000_000:
        raise ValueError('Broadening would require >5 million wavelength samples; crop the model first')
    loggrid = np.linspace(np.log(wave[0]), np.log(wave[-1]), n)
    grid = np.exp(loggrid)
    sampled = np.interp(grid, wave, f)
    velocity_step = np.diff(loggrid)[0] * C_KMS
    sigma_pixels = C_KMS / resolution / (2 * np.sqrt(2 * np.log(2))) / velocity_step
    smoothed = gaussian_filter1d(sampled, sigma=sigma_pixels, mode='nearest')
    return grid, smoothed
