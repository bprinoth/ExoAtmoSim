"""Keplerian orbital traces for synthetic observations, no optional packages.

Keplerian RV here is a point-mass approximation. Tc is the inferior-
conjunction epoch for an edge-on orbit; obliquity/transit contact geometry,
BERV and RM require separate treatment. Angles omega are in degrees.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class RVTraceEstimator:
    period_days: float
    tc_bjd: float
    eccentricity: float
    omega_deg: float
    stellar_k_kms: float
    mass_ratio_star_to_planet: float

    def __post_init__(self):
        if not np.all(np.isfinite([self.period_days, self.tc_bjd, self.eccentricity, self.stellar_k_kms, self.mass_ratio_star_to_planet])):
            raise ValueError('Orbital parameters must be finite')
        if self.eccentricity > 0 and not np.isfinite(self.omega_deg):
            raise ValueError('An eccentric orbit requires a finite omega in degrees')
        if self.period_days <= 0 or not 0 <= self.eccentricity < 1:
            raise ValueError('Require period > 0 days and 0 <= ecc < 1')
        if self.mass_ratio_star_to_planet <= 0 or self.stellar_k_kms < 0:
            raise ValueError('Mass ratio must be positive and stellar K non-negative')

    def evaluate(self, bjd: np.ndarray) -> dict[str, np.ndarray]:
        t = np.asarray(bjd, dtype=float)
        if t.ndim != 1 or not np.isfinite(t).all():
            raise ValueError('BJD times must be a finite one-dimensional array')
        e = self.eccentricity
        omega = np.deg2rad(90.0 if e == 0.0 else self.omega_deg)
        f_conj = np.pi / 2 - omega
        ecc_at_conj = 2 * np.arctan2(np.sqrt(1-e) * np.sin(f_conj/2),
                                      np.sqrt(1+e) * np.cos(f_conj/2))
        m_at_conj = ecc_at_conj - e * np.sin(ecc_at_conj)
        m = m_at_conj + 2*np.pi*(t - self.tc_bjd)/self.period_days
        ecc_anomaly = np.mod(m + np.pi, 2*np.pi) - np.pi
        for _ in range(30):
            # use wrapped mean anomaly to avoid loss of precision across long baselines
            target = (m + np.pi) % (2*np.pi) - np.pi
            update = (ecc_anomaly - e*np.sin(ecc_anomaly) - target)/(1-e*np.cos(ecc_anomaly))
            ecc_anomaly -= update
            if np.max(np.abs(update), initial=0) < 1e-13:
                break
        anomaly = 2*np.arctan2(np.sqrt(1+e)*np.sin(ecc_anomaly/2),
                               np.sqrt(1-e)*np.cos(ecc_anomaly/2))
        star = self.stellar_k_kms * (np.cos(anomaly+omega) + e*np.cos(omega))
        planet = -star * self.mass_ratio_star_to_planet
        return {'true_anomaly_rad': anomaly, 'phase': ((t-self.tc_bjd)/self.period_days + 0.5) % 1 - 0.5,
                'rv_star_kms': star, 'rv_planet_kms': planet}

    @classmethod
    def from_config(cls, config):
        # Existing configured stellar K is in km/s; planet/star masses are Jupiter/Solar masses.
        ms = config.get_float('Ms')
        mp = config.get_float('Mp')
        return cls(config.get_float('P'), config.get_float('Tc'),
                   config.get_float('ecc', 0.), config.get_float('omega', 90.),
                   config.get_float('K'), ms / mp * 1047.348644)
