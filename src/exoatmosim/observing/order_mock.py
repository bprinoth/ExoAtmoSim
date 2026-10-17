"""Independent, legacy-compatible sampled spectral-order synthesis.

This is a deterministic-first reference implementation: users supply the
stellar flux and planetary transmission on the SAME wavelength convention as
the order definition, plus per-exposure Doppler velocities and occultation
weights. It does not use external_mock. The caller may apply rotation and instrumental convolution to the input
spectra before calling this function. Exposure smearing and
Rossiter-McLaughlin distortions are NOT included. Optional telluric transmission may be applied before noise.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

C_KMS = 299792.458


def simulate_orders(order_wavelengths_nm, stellar_wave_nm, stellar_flux,
                    planet_wave_nm, planet_transmission, rv_star_kms,
                    rv_planet_kms, transit_weight, *, snr=100.0, seed=42,
                    blaze=None, add_noise=True, stellar_observer_frame=False, telluric_wavelength_nm=None,
                    telluric_transmission=None):
    """Return (wavelength_2d, flux_2d, error_2d) per order.

    ``planet_transmission`` is a dimensionless planet+star factor near unity,
    with 1 outside absorption. For transmission, F = F_star *
    [1 + transit_weight * (F_planet - 1)]. The model is Doppler-shifted using
    λ_observed = λ_rest (1+v/c) for small v/c, matching previous implementation's convention.
    ``transit_weight`` is an explicit 0..1 occultation mask, NOT a phase angle.
    SNR applies to the final normalized flux (absolute noise = flux/SNR),
    consistent with legacy previous external mock implementation.make_orders.
    """
    orders = [np.asarray(w, dtype=float) for w in order_wavelengths_nm]
    sw, sf = np.asarray(stellar_wave_nm, float), np.asarray(stellar_flux, float)
    pw, pf = np.asarray(planet_wave_nm, float), np.asarray(planet_transmission, float)
    vs = np.asarray(rv_star_kms, float)
    vp = np.asarray(rv_planet_kms, float)
    weights = np.asarray(transit_weight, float)
    if not orders or any(w.ndim != 1 or len(w) < 2 for w in orders):
        raise ValueError('Expected per-segment wavelength arrays with at least two samples')
    if stellar_observer_frame:
        if sf.shape != (len(vs), len(sw)) or not np.isfinite(sf).all():
            raise ValueError('Observer-frame stellar flux must be finite Nexp x Nwave')
        if telluric_wavelength_nm is not None or telluric_transmission is not None:
            raise ValueError('Observer-frame VIOLA stellar input already contains tellurics; do not apply TelFit again')
        inputs = ((pw, pf),)
        if sw.ndim != 1 or sw.size < 2 or not np.isfinite(sw).all() or np.any(sw <= 0) or np.any(np.diff(sw) <= 0):
            raise ValueError('Observer-frame stellar wavelengths must be finite and increasing')
    else:
        inputs = ((sw, sf), (pw, pf))
    for w, flux in inputs:
        if w.ndim != 1 or flux.shape != w.shape or len(w) < 2 or not np.isfinite(w).all() or not np.isfinite(flux).all() or np.any(w <= 0) or np.any(np.diff(w) <= 0):
            raise ValueError('Stellar and planetary spectra must be finite, increasing, matching 1D arrays')
    if vs.ndim != 1 or len(vs) < 1 or vp.shape != vs.shape or weights.shape != vs.shape:
        raise ValueError('Velocities and occultation mask must be matching 1D exposure arrays')
    if not np.isfinite(np.r_[vs, vp, weights]).all() or np.any((weights < 0) | (weights > 1)):
        raise ValueError('Non-finite velocities or invalid 0..1 transit weights')
    if not np.isfinite(snr) or snr <= 0:
        raise ValueError('SNR must be finite and positive')
    if blaze is not None:
        blaze = np.asarray(blaze, float)
        if len({len(w) for w in orders}) != 1 or blaze.shape != (len(orders), len(orders[0])) or not np.isfinite(blaze).all():
            raise ValueError('Blaze must have same shape as wavelength orders')
    tw = None
    tt = None
    if (telluric_wavelength_nm is None) != (telluric_transmission is None):
        raise ValueError('Provide both telluric_wavelength_nm and telluric_transmission')
    if telluric_wavelength_nm is not None:
        tw = np.asarray(telluric_wavelength_nm, float)
        tt = np.asarray(telluric_transmission, float)
        if (tw.ndim != 1 or len(tw) < 2 or not np.isfinite(tw).all()
                or np.any(np.diff(tw) <= 0) or np.any(tw <= 0)
                or tt.shape != (len(vs), len(tw)) or not np.isfinite(tt).all()
                or np.any(tt < 0) or np.any(tt > 1)):
            raise ValueError('Telluric model must be finite, air-nm increasing grid and Nexp x Nwave transmission within [0,1]')
    rng = np.random.default_rng(seed)
    results = []
    for i, wl in enumerate(orders):
        if not np.isfinite(wl).all() or np.any(wl <= 0) or np.any(np.diff(wl) <= 0):
            raise ValueError(f'Order {i} has invalid wavelength grid')
        star_at = wl[None, :] if stellar_observer_frame else wl[None, :] / (1 + vs[:, None]/C_KMS)
        planet_at = wl[None, :] / (1 + vp[:, None]/C_KMS)
        # Reject extrapolation rather than silently creating fake edge pixels.
        if star_at.min() < sw[0] or star_at.max() > sw[-1] or planet_at.min() < pw[0] or planet_at.max() > pw[-1]:
            raise ValueError(f'Order {i} is not fully covered by star/planet spectra at Doppler-shifted wavelengths')
        if stellar_observer_frame:
            stellar = np.stack([np.interp(wl, sw, sf[j]) for j in range(len(vs))])
        else:
            stellar = np.interp(star_at, sw, sf)
        planet = np.interp(planet_at, pw, pf)
        flux = stellar * (1.0 + weights[:, None] * (planet - 1.0))
        if tw is not None:
            if wl[0] < tw[0] or wl[-1] > tw[-1]:
                raise ValueError(f'Telluric model does not cover order {i}; refusing extrapolation')
            telluric_on_order = np.stack([np.interp(wl, tw, curve) for curve in tt])
            flux *= telluric_on_order
        if blaze is not None:
            flux *= blaze[i][None, :]
        sigma = np.abs(flux) / float(snr)
        if add_noise:
            flux = flux + rng.normal(size=flux.shape) * sigma
        results.append((np.broadcast_to(wl, flux.shape).copy(), flux, sigma))
    return results


def write_spectral_orders(results, directory: str | Path, *, overwrite=False):
    """Write order_N.fits / wave_N.fits; reject missing or stale outputs."""
    from astropy.io import fits
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    planned = []
    for i, (wave, flux, sigma) in enumerate(results):
        for kind, data in (('order', flux), ('wave', wave)):
            filename = directory / f'{kind}_{i}.fits'
            if filename.exists() and not overwrite:
                raise FileExistsError(filename)
            planned.append((filename, np.asarray(data)))
    for filename, data in planned:
        fits.writeto(filename, data, overwrite=overwrite)
    for filename, data in planned:
        with fits.open(filename) as hdul:
            if hdul[0].data.shape != data.shape:
                raise RuntimeError(f'Mock output verification failed: {filename}')
    return [filename for filename, _ in planned]
