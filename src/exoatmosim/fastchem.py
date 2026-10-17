"""Optional FastChem equilibrium chemistry bridge.

No chemistry values, bulk fractions or elemental composition are fabricated.
Input abundance and thermochemical files are supplied by the user.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from .config import Config
from .atmosphere import Atmosphere
from .profiles import AbundanceProfile


# FastChem gas identifiers use Hill-order elements and explicit atom counts.
# The user-facing names remain unchanged for POSEIDON and abundance tables.
# Electrons are not atoms; database lookup is exact for this species.
def fastchem_hill_name(name: str) -> str:
    """Convert molecular formulae to FastChem Hill notation, preserving charge.

    ChemFormula parses elemental composition; FastChem requires count=1 explicitly.
    This conversion cannot distinguish structural isomers or electronic states.
    """
    import re
    if name == 'e-':
        return name
    match = re.fullmatch(r'([A-Z][A-Za-z0-9()]*?)([+-]+)?', name)
    if not match:
        raise ValueError(f'Unsupported chemical formula {name!r}; use an exact FastChem identifier')
    formula, suffix = match.groups()
    # Monatomic spectroscopic ion notation Fe2+ means Fe with charge +2,
    # not a two-atom iron molecular ion.
    mono = re.fullmatch(r'([A-Z][a-z]?)([2-9][0-9]*)', formula)
    if mono and suffix in ('+', '-'): 
        formula = mono.group(1)
        suffix = suffix * int(mono.group(2))
    charge = (suffix or '').count('+') - (suffix or '').count('-')
    try:
        from chemformula import ChemFormula
    except ImportError as exc:
        raise ImportError('FastChem species conversion requires chemformula: pip install chemformula') from exc
    try:
        parsed = ChemFormula(formula, charge=charge)
        elements = parsed.element
        if not elements or any(int(n) != n or n <= 0 for n in elements.values()):
            raise ValueError('invalid atom counts')
    except Exception as exc:
        raise ValueError(f'Could not parse chemical formula {name!r}') from exc
    if 'C' in elements:
        ordering = [a for a in ('C', 'H') if a in elements] + sorted(a for a in elements if a not in ('C', 'H'))
    else:
        ordering = sorted(elements)
    return ''.join(f'{el}{int(elements[el])}' for el in ordering) + ('+' * charge if charge > 0 else '-' * -charge)


def resolve_fastchem_species_index(engine, name: str, number_of_species: int) -> tuple[int, str]:
    """Find an exact or Hill-form gas species in the active FastChem database.

    Exact database names win. Refuse structural-isomer ambiguity instead of
    pretending a Hill composition uniquely identifies a molecular structure.
    """
    exact = int(engine.getGasSpeciesIndex(name))
    if 0 <= exact < number_of_species:
        return exact, name
    candidate = fastchem_hill_name(name)
    # HCN and HNC share the same Hill composition but differ structurally.
    # They require exact identifiers, not composition-based resolution.
    if candidate == 'C1H1N1':
        raise ValueError(f'{name!r} has an isomer-ambiguous Hill composition ({candidate}); '
                         'specify the exact FastChem gas database identifier')
    index = int(engine.getGasSpeciesIndex(candidate))
    if 0 <= index < number_of_species:
        return index, candidate
    raise ValueError(f'FastChem has no gas species for {name!r} (tried {name!r}, {candidate!r}); '
                     'check the installed logK database')


def equilibrium_profile(config: Config, atmosphere: Atmosphere, *, config_dir: str | Path,
                        engine_module=None) -> AbundanceProfile:
    """Calculate linear VMR on the atmospheric P--T grid via PyFastChem.

    Supported Python API: FastChem 3-style ``FastChem``, ``FastChemInput`` and
    ``FastChemOutput``. The installed version must be independently validated.
    """
    if engine_module is None:
        try:
            import pyfastchem as engine_module
        except ImportError as exc:
            raise ImportError('FastChem requested: install pyfastchem and supply chemistry data files') from exc
    root = Path(config_dir)
    inputs = []
    for key in ('fastchem_elements_file', 'fastchem_species_file'):
        if key not in config:
            raise ValueError(f'{key} required when chemistry fastchem')
        path = Path(config[key]).expanduser()
        if not path.is_absolute():
            path = root / path
        if not path.is_file():
            raise FileNotFoundError(f'{key}: {path}')
        inputs.append(str(path.resolve()))
    if 'bulk_species' not in config:
        raise ValueError('Explicit bulk_species required for FastChem')
    if 'fastchem_species' not in config:
        raise ValueError('Explicit fastchem_species required, e.g. H2 He H- Fe Ti')
    species = config['fastchem_species'].split()
    if not species or len(set(species)) != len(species):
        raise ValueError('fastchem_species must contain unique names')
    bulk = tuple(config['bulk_species'].split())
    if not set(bulk).issubset(species):
        raise ValueError('bulk_species must be included in fastchem_species')
    condensation = config.get_bool('fastchem_condensation', False)
    rainout = config.get_bool('fastchem_rainout', False)
    if rainout:
        raise NotImplementedError('FastChem Cond rainout needs validated layer-ordering and condensate API')
    if condensation:
        if 'fastchem_condensates_file' not in config:
            raise ValueError('fastchem_condensates_file required for equilibrium condensation')
        condensates = Path(config['fastchem_condensates_file']).expanduser()
        if not condensates.is_absolute():
            condensates = root / condensates
        if not condensates.is_file():
            raise FileNotFoundError(f'fastchem_condensates_file: {condensates}')

    try:
        engine = (engine_module.FastChem(inputs[0], inputs[1], str(condensates.resolve()), 0)
                  if condensation else engine_module.FastChem(inputs[0], inputs[1], 0))
        fc_input = engine_module.FastChemInput()
        fc_output = engine_module.FastChemOutput()
    except Exception as exc:
        raise RuntimeError('PyFastChem API incompatible with this adapter; validate your installed version') from exc
    fc_input.temperature = list(np.asarray(atmosphere.temperature_k, float))
    fc_input.pressure = list(np.asarray(atmosphere.pressure_bar, float))
    fc_input.equilibrium_condensation = condensation
    result = engine.calcDensities(fc_input, fc_output)
    success = getattr(engine_module, 'FASTCHEM_SUCCESS', 0)
    if result != success:
        raise RuntimeError(f'FastChem did not converge successfully (status={result})')
    density = np.asarray(fc_output.number_densities, float)
    n_levels = len(atmosphere.pressure_bar)
    if density.ndim != 2 or density.shape[0] != n_levels:
        raise RuntimeError(f'Unexpected FastChem number_densities shape: {density.shape}')
    # P/(k_B T) yields physical total gas number density [cm^-3], P in bar.
    # This includes all particles, rather than renormalising selected species.
    boltzmann_erg_per_k = 1.380649e-16
    n_gas = np.asarray(atmosphere.pressure_bar) * 1e6 / (boltzmann_erg_per_k * atmosphere.temperature_k)
    ratios = {}
    for name in species:
        index, _backend_name = resolve_fastchem_species_index(engine, name, density.shape[1])
        ratios[name] = density[:, index] / n_gas
    return AbundanceProfile(atmosphere.pressure_bar, ratios, bulk)
