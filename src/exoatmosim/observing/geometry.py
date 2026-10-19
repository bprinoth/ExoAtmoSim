"""Dependency-free projected uniform-disc occultation geometry."""
import numpy as np

def occulted_fraction(d, Rp, zp):
    """Fraction of the *stellar disc* blocked by a foreground planet.

    Parameters
    ----------
    d : float
        Projected centre-to-centre separation in stellar-radius units.
        Do NOT subtract the planet radius before calling this function.
    Rp : float
        Planet-to-star radius ratio (positive).
    zp : float
        Line-of-sight position in the RVTraceEstimator convention:
        zp >= 0 means the planet is in front of the star.

    Returns
    -------
    float
        Fractional loss of stellar light, between 0 and min(Rp**2, 1).
        A planet behind the star contributes *zero* loss in a transit
        light curve; secondary eclipses need a separate emission model.

    Assumes uniform stellar surface brightness (no limb darkening).
    """
    if not np.all(np.isfinite([d, Rp, zp])) or d < 0 or Rp <= 0:
        raise ValueError(f'Occultation geometry requires finite d >= 0 [is {d}] and Rp > 0 [is {Rp}]')
    if zp < 0 or d >= 1.0 + Rp:
        return 0.0
    # One disc lies completely inside the other.
    if d <= abs(1.0 - Rp):
        return float(min(Rp * Rp, 1.0))

    # Standard circle-circle intersection area (in units of Rstar^2).
    arg1 = np.clip((d*d + Rp*Rp - 1.0)/(2.0*d*Rp), -1.0, 1.0)
    arg2 = np.clip((d*d + 1.0 - Rp*Rp)/(2.0*d), -1.0, 1.0)
    radical = max(0.0, (-d + Rp + 1.0)*(d + Rp - 1.0)*
                   (d - Rp + 1.0)*(d + Rp + 1.0))
    area = Rp*Rp*np.arccos(arg1) + np.arccos(arg2) - 0.5*np.sqrt(radical)
    return float(np.clip(area / np.pi, 0.0, min(Rp*Rp, 1.0)))

