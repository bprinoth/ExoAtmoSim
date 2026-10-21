"""Generate ExoAtmoSim planet configurations using NASA Exoplanet Archive.

Astroquery/Astropy are optional dependencies, only needed for archive queries.
This module intentionally retains the keys and default values of the user's
existing configuration generator.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
from typing import Any

from .config import Config

ALLOWED_SPECTRAL_TYPES = ('F0V', 'F7V', 'G8V', 'K2V', 'K7V', 'M1A', 'M3A', 'M8A', 'Sun')


def closest_spectral_type(spectral_type: Any) -> str:
    """Map a catalogue spectral type to an existing profile template."""
    if spectral_type is None:
        return 'Sun'
    value = str(spectral_type).strip()
    if not value or value.lower() in ('nan', 'none', '--', 'unknown'):
        return 'Sun'
    for allowed in ALLOWED_SPECTRAL_TYPES:
        if value.lower() == allowed.lower():
            return allowed
    if 'sun' in value.lower():
        return 'Sun'
    match = re.search(r'([OBAFGKM])\s*(\d+(?:\.\d+)?)', value.upper())
    if not match:
        return 'Sun'
    base = {'F': 0, 'G': 10, 'K': 20, 'M': 30}
    if match.group(1) not in base:
        return 'Sun'
    position = base[match.group(1)] + float(match.group(2))
    candidates = {'F0V': 0, 'F7V': 7, 'Sun': 12, 'G8V': 18, 'K2V': 22,
                  'K7V': 27, 'M1A': 31, 'M3A': 33, 'M8A': 38}
    return min(candidates, key=lambda name: abs(candidates[name] - position))


def _value(row: Any, key: str, unit: Any = None, default: Any = math.nan) -> Any:
    try:
        value = row[key]
    except (KeyError, IndexError):
        return default
    if value is None:
        return default
    # Astropy masked columns/quantities; import only on archive use.
    import numpy as np
    if np.ma.is_masked(value):
        return default
    if unit is not None and hasattr(value, 'to'):
        try:
            value = value.to(unit)
        except (ValueError, TypeError):
            return default
    if hasattr(value, 'value'):
        value = value.value
    if np.ma.is_masked(value):
        return default
    if isinstance(value, bytes):
        value = value.decode('utf-8')
    if isinstance(value, (int, float, np.number)):
        value = float(value)
        return value if math.isfinite(value) else default
    return value


def _flag(row: Any, key: str, default: bool = False) -> bool:
    value = _value(row, key)
    try:
        return bool(value) if math.isfinite(float(value)) else default
    except (TypeError, ValueError):
        return default


def parameters_from_row(row: Any, *, resolution: int = 300_000,
                        air: bool = True) -> dict[str, Any]:
    """Convert one Archive row to the full legacy set of configuration keys."""
    import astropy.constants as const
    import astropy.units as u
    from astropy.coordinates import Angle

    def val(key: str, unit: Any = None, default: Any = math.nan) -> Any:
        return _value(row, key, unit, default)

    mass, radius = val('pl_bmassj'), val('pl_radj')
    if math.isfinite(mass) and math.isfinite(radius) and mass > 0 and radius > 0:
        gravity = (const.G * (mass * u.M_jup) / (radius * u.R_jup)**2)
        loggp = math.log10(gravity.to_value(u.cm / u.s**2))
    else:
        loggp = math.nan

    ra, dec = val('ra'), val('dec')
    ra_text = (Angle(ra, unit=u.deg).to_string(unit=u.hourangle, sep=':',
               precision=2, pad=True) if math.isfinite(ra) else 'nan')
    dec_text = (Angle(dec, unit=u.deg).to_string(unit=u.deg, sep=':',
                precision=2, pad=True, alwayssign=True) if math.isfinite(dec) else 'nan')
    vsys = val('st_radv', u.km / u.s)
    if not math.isfinite(vsys):
        vsys = 0.0  # Preserved legacy default
    duration = val('pl_trandur')
    return {
        'P': val('pl_orbper'), 'a': val('pl_orbsmax'), 'aRstar': val('pl_ratdor'),
        'Rp': radius, 'Mp': mass, 'K': val('pl_rvamp', u.km / u.s),
        'vsini': val('st_vsin', u.km / u.s), 'RpRstar': val('pl_ratror'),
        'vsys': vsys, 'inclination': val('pl_orbincl'), 'Tc': val('pl_tranmid'),
        'Tp': val('pl_orbtper'), 'T14': duration, 'RA': ra_text, 'DEC': dec_text,
        'Teq': val('pl_eqt'), 'Teff': val('st_teff'), 'ecc': val('pl_orbeccen'),
        'impact': val('pl_imppar'), 'omega': val('pl_orblper'),
        'Rs': val('st_rad'), 'Ms': val('st_mass'), 'logg': val('st_logg'),
        'metallicity': val('st_met'),
        'stellar_type': str(val('st_spectype', default='nan') or 'nan').strip() or 'nan',
        'Hmag': val('sy_hmag'), 'Kmag': val('sy_kmag'), 'Jmag': val('sy_jmag'),
        'Gyr_age': val('st_age'), 'loggp': loggp,
        'lampoo': val('pl_projobliq'), 'resolution': resolution,
        'observatory': 'paranal', 'air': air,
        'transit_flag': _flag(row, 'tran_flag'),
        'lat': -24.6276, 'long': -70.4051, 'elev': 2635,
        'observation_window_hours': duration + 2.0 if math.isfinite(duration) and duration > 0 else math.nan,
        'readout': 60, 'exposure_time': 150,
    }


# Comments are parsed by ExoAtmoSim Config.
_GROUPS = (
    ('Planet and orbit (catalogue values; verify missing nan values)',
     'P a aRstar Rp Mp K RpRstar inclination Tc Tp T14 Teq ecc impact omega loggp'.split()),
    ('Host star and system',
     'RA DEC Rs Ms Teff logg metallicity stellar_type vsini vsys Hmag Kmag Jmag Gyr_age lampoo'.split()),
    ('Observatory and observing defaults (adapt to planned observations)',
     'resolution observatory air transit_flag lat long elev observation_window_hours readout exposure_time'.split()),
)


# Only defaults that describe the workflow belong in every planet config.
# Local installation paths are intentionally supplied by --site-config.
_WORKFLOW_GROUPS = (
    ('Atmospheric structure', {
        'tp_profile': 'isothermal', 'pressure_min': '1e-15',
        'pressure_max': '10', 'n_pressure_layers': '200',
    }),
    ('Chemistry: FastChem gas equilibrium or manual abundance file', {
        'chemistry': 'fastchem', 'bulk_species': 'H2 He',
        'fastchem_species': 'H2 He Fe', 'fastchem_condensation': 'False',
    }),
    ('Radiative transfer and instrument', {
        'instrument': 'ESPRESSO', 'instrument_margin_nm': '5',
    }),
    ('Orbital velocities and transit geometry', {
        'rv_plot_extent': '200',
    }),
    ('Observing atmosphere and tellurics', {
        'tellurics': 'False',
        'pwv': '2.0',
        'pwv_sigma': '0.1',
    }),
)

# Site-level values are copied only when explicitly provided. The original
# config reader is deliberately whitespace based, so no quoted paths/spaces.
_SITE_KEYS = frozenset({
    'fastchem_elements_file', 'fastchem_species_file', 'fastchem_condensates_file',
    'poseidon_path', 'poseidon_input_data', 'pysyn_cdbs',
    'instrument_order_file', 'instrument_library_path',
    'kaltenegger_profiles_dir', 'telluric_grid_dir', 'viola_skycalc_dir',
})

# Local profile path must not be embedded in a planet-specific science config.
_SITE_LOCAL_ONLY_KEYS = frozenset({'kaltenegger_profiles_dir', 'telluric_grid_dir', 'viola_skycalc_dir'})


def load_site_config(path: str | Path) -> dict[str, str]:
    """Read optional installation-specific backend paths, with strict key validation."""
    site = Config.from_file(path)
    unknown = set(site) - _SITE_KEYS
    if unknown:
        raise ValueError('Unknown site-config keys: ' + ', '.join(sorted(unknown)))
    return dict(site)


def format_science_config(params: dict[str, Any], *,
                          site_values: dict[str, str] | None = None,
                          include_workflow: bool = True) -> str:
    """Generate a readable complete config without changing catalogue values."""
    lines = ['# ExoAtmoSim planet configuration',
             '# NASA Exoplanet Archive parameters; verify missing/nan entries.',
             '# Physical and observing parameters; edit workflow defaults below.', '']
    written = set()
    for heading, keys in _GROUPS:
        lines.append(f'# --- {heading} ---')
        for key in keys:
            if key in params:
                lines.append(Config.from_mapping({key: params[key]}).to_text().strip())
                written.add(key)
        lines.append('')
    if include_workflow:
        for heading, defaults in _WORKFLOW_GROUPS:
            lines.append(f'# --- {heading} ---')
            if heading == 'Observing atmosphere and tellurics':
                lines.append('# Set tellurics True to use the external HDF5 telluric library.')
                lines.append('# PWV mean and sigma are in mm; the per-exposure series is reproducible.')
            if heading.startswith('Chemistry'):
                lines.append('# For manual profiles: chemistry file; abundances_file abundances_manual.dat')
                lines.append('# For FastChem Cond: set fastchem_condensation True and supply condensates file')
            for key, value in defaults.items():
                if key not in written and key not in params:
                    lines.append(f'{key}\t{value}')
                    written.add(key)
            lines.append('')
    if site_values:
        unknown = set(site_values) - _SITE_KEYS
        if unknown:
            raise ValueError('Unknown site-config keys: ' + ', '.join(sorted(unknown)))
        lines.append('# --- Local backend and database paths (from site config) ---')
        for key, value in site_values.items():
            if key in _SITE_LOCAL_ONLY_KEYS:
                continue
            if key not in params:
                lines.append(f'{key}\t{value}')
                written.add(key)
        lines.append('')
    for key, value in params.items():
        if key not in written:
            lines.append(Config.from_mapping({key: value}).to_text().strip())
    return '\n'.join(lines).rstrip() + '\n'


def find_parameters(planetname: str, resolution: int = 300_000,
                    air: bool = True, print_params: bool = True,
                    *, data_dir: str | Path = 'data',
                    overwrite: bool = False, preserve_existing: bool = True,
                    site_config: str | Path | None = None,
                    instrument: str = 'ESPRESSO') -> Config:
    """Generate a complete ExoAtmoSim planet config from NASA parameters.

    Existing files are protected by default. --overwrite updates only archive
    values, preserving user edits. --fresh intentionally regenerates everything.
    """
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

    safe_name = planetname.replace("'", "''")
    rows = NasaExoplanetArchive.query_criteria(
        table='pscomppars', select='*', where=f"pl_name='{safe_name}'")
    if rows is None or len(rows) == 0:
        raise ValueError(f'No NASA Exoplanet Archive result for {planetname!r}')
    if len(rows) > 1:
        print(f'[WARNING] Found {len(rows)} rows for {planetname!r}; using first')
    params = parameters_from_row(rows[0], resolution=resolution, air=air)
    if print_params:
        for key, value in params.items():
            print(f'{key:15} {value}')
    path = Path(data_dir) / planetname.replace(' ', '_') / 'config'
    if path.exists() and not overwrite:
        raise FileExistsError(f'{path} already exists; pass overwrite=True to regenerate')
    if path.exists() and overwrite and preserve_existing:
        config = Config.from_file(path)
        for key, value in params.items():
            config.set(key, value)
        # Do not add workflow defaults during a catalogue refresh: an existing
        # hand-edited configuration must be preserved exactly apart from NASA keys.
    else:
        from .instruments import instrument_resource
        instrument_resource(instrument)  # Reject missing/unsupported instrument orders.
        site_values = load_site_config(site_config) if site_config else None
        text = format_science_config(params, site_values=site_values)
        config = Config.from_text(text)
        config.set('instrument', instrument)
    config.write(path)
    print(f'[INFO] Config written to {path}')
    return config


def main() -> None:
    parser = argparse.ArgumentParser(description='Generate key/value config from NASA Exoplanet Archive')
    parser.add_argument('planet', help="Planet name, e.g. 'WASP-121 b'")
    parser.add_argument('--resolution', type=int, default=300_000)
    parser.add_argument('--no-air', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--instrument', default='ESPRESSO',
                        help='Instrument with bundled order grids (e.g. CRIRES_K)')
    parser.add_argument('--site-config', type=Path,
                        help='Optional reusable file of machine-specific backend paths')
    parser.add_argument('--overwrite', action='store_true',
                        help='Refresh NASA keys while retaining existing custom settings')
    parser.add_argument('--fresh', action='store_true',
                        help='With --overwrite, replace entire config (removes custom settings)')
    args = parser.parse_args()
    find_parameters(args.planet, resolution=args.resolution, air=not args.no_air,
                    print_params=args.verbose, overwrite=args.overwrite,
                    preserve_existing=not args.fresh, site_config=args.site_config,
                    instrument=args.instrument)


if __name__ == '__main__':
    main()
