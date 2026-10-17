"""Write the planetary model into a three-HDU air-wavelength FITS file."""
from pathlib import Path
def write_spectrum_fits(spectrum, destination: str | Path, *, overwrite: bool = False,
                             config_path: str | Path | None = None) -> Path:
    """Write the validated three-HDU planetary model FITS format.

    Schema established from the user's working TRAPPIST-1_b.fits example:
    PrimaryHDU([0.0]), ImageHDU(wavelength_nm), ImageHDU(transmission).
    Input wavelength: vacuum microns (POSEIDON). FITS wavelength: air nm
    (matches the supplied ESPRESSO order definitions). Conversion uses the
    internal vactoair exactly once. The transmission array
    is already 1 - POSEIDON output and must not be inverted a second time.
    """
    import numpy as np
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError('Writing planetary FITS requires astropy') from exc

    vacuum_nm = np.asarray(spectrum.wavelength_um, dtype=np.float64) * 1000.0
    if vacuum_nm.ndim != 1 or vacuum_nm.size < 2 or not np.isfinite(vacuum_nm).all():
        raise ValueError('Input vacuum wavelength grid must be finite and one-dimensional')
    if not np.all(np.diff(vacuum_nm) > 0) or vacuum_nm[0] <= 0:
        raise ValueError('Input vacuum wavelength grid must be positive and strictly ascending')
    from .wavelengths import vactoair
    wavelength_nm = np.asarray(vactoair(vacuum_nm), dtype=np.float64)
    transmission = np.asarray(spectrum.transmission, dtype=np.float64)
    if wavelength_nm.ndim != 1 or wavelength_nm.size < 2 or transmission.shape != wavelength_nm.shape:
        raise ValueError('Expected matching 1-D wavelength and transmission arrays')
    if not np.isfinite(wavelength_nm).all() or not np.all(np.diff(wavelength_nm) > 0):
        raise ValueError('FITS wavelength grid must be finite, positive and strictly ascending')
    if wavelength_nm[0] <= 0 or not np.isfinite(transmission).all():
        raise ValueError('FITS transmission and wavelengths must be finite with positive wavelengths')
    destination = Path(destination)
    if destination.exists() and not overwrite:
        raise FileExistsError(f'Refusing to overwrite existing planetary model FITS: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    hdus = fits.HDUList([
        fits.PrimaryHDU(np.array([0.0], dtype=np.float64)),
        fits.ImageHDU(wavelength_nm),
        fits.ImageHDU(transmission),
    ])
    hdus.writeto(destination, overwrite=overwrite, output_verify='exception')
    return destination

