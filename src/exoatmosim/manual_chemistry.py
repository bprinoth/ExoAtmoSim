"""Constant manual VMR profiles from a single planet config."""
from __future__ import annotations
import numpy as np
from .profiles import AbundanceProfile

def manual_profile(config, atmosphere):
    """Keys: chemistry manual; bulk_species H2 He; manual_abundances 'H2=0.85 He=0.15'.

    User supplies all species, including bulk. No inference or normalisation.
    """
    if 'manual_abundances' not in config:
        raise ValueError('chemistry manual requires manual_abundances (species=VMR ...)')
    if 'bulk_species' not in config:
        raise ValueError('chemistry manual requires explicit bulk_species')
    entries = config['manual_abundances'].split()
    vmrs = {}
    for entry in entries:
        if entry.count('=') != 1:
            raise ValueError('manual_abundances entries must have species=VMR form')
        species, val = entry.split('=', 1)
        if not species or species in vmrs:
            raise ValueError(f'Missing or duplicate manual species {species!r}')
        try:
            vmr = float(val)
        except ValueError as exc:
            raise ValueError(f'Invalid VMR for {species}') from exc
        if not np.isfinite(vmr) or not 0 <= vmr <= 1:
            raise ValueError(f'VMR for {species} must be finite in [0,1]')
        vmrs[species] = np.full_like(atmosphere.pressure_bar, vmr, dtype=float)
    return AbundanceProfile(atmosphere.pressure_bar, vmrs, tuple(config['bulk_species'].split()))
