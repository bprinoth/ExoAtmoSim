"""Air/vacuum conversions, wavelength input and output in nanometres.

Dispersion relations follow the established astronomical formulas used by
the prior model exporter. Input must be finite and strictly positive.
"""
import numpy as np

def _wavelength(value):
    arr = np.asarray(value, dtype=float)
    if not np.isfinite(arr).all() or np.any(arr <= 0):
        raise ValueError("Wavelengths must be positive and finite")
    return arr

def airtovac(wlnm):
    w = _wavelength(wlnm)
    s = 1000.0 / w
    n = 1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    return w * n

def vactoair(wlnm):
    w = _wavelength(wlnm)
    s = 1000.0 / w
    f = 1.0 + 5.792105e-2 / (238.0185 - s**2) + 1.67917e-3 / (57.362 - s**2)
    return w / f
