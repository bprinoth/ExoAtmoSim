"""Experimental POSEIDON transmission adapter based on the legacy run_poseidon.py.

No POSEIDON import occurs until compute_transmission() is called. Native
POSEIDON-dependent numerical behavior must be validated against a legacy run.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np

from .config import Config
from .objects import Planet, Star
from .atmosphere import Atmosphere
from .backend_paths import configure_backend_path, configure_poseidon_environment

# Opacity names from the user's existing POSEIDON script; this list reflects
# the user's database, and should not be interpreted as universal support.
SUPPORTED_SPECIES = frozenset('''H2O CO2 CH4 CO Na K NH3 HCN SO2 H2S PH3 C2H2 OCS TiO VO
AlO SiO CaO MgO NaO LaO ZrO SO NO PO TiH CrH FeH ScH AlH SiH BeH CaH
MgH LiH NaH OH OH+ CH NH SH PN PS CS C2 CH3 H3+ N2O NO2 C2H4 C2H6
CH3CN CH3OH CH3Cl GeH4 CS2 O2 O3 C2H6S Al Ba Ba+ Ca Ca+ Cr Cs Fe
Fe+ Li Mg Mg+ Mn Ni O Rb Sc Ti Ti+ V V+ N2 12C-16O 13C-16O 12C-18O
12C-17O 13C-18O 13C-17O H H2 He H- e-'''.split())


@dataclass(frozen=True)
class Composition:
    """Linear volume mixing ratios on the ascending-pressure atmosphere grid."""
    mixing_ratios: dict[str, np.ndarray]

    def __post_init__(self):
        if not self.mixing_ratios:
            raise ValueError('No chemistry species provided')
        checked = {}
        lengths = set()
        for name, values in self.mixing_ratios.items():
            a = np.array(values, dtype=float, copy=True)
            if a.ndim != 1 or not np.all(np.isfinite(a)) or np.any(a < 0):
                raise ValueError(f'Invalid volume mixing ratios for {name}')
            lengths.add(len(a))
            a.setflags(write=False)
            checked[name] = a
        if len(lengths) != 1:
            raise ValueError('Chemistry profiles must share a pressure grid')
        object.__setattr__(self, 'mixing_ratios', checked)


def read_kaltenegger_pt(path: str | Path) -> Atmosphere:
    """Read legacy `Z,T,P` records: altitude [km], T [K], P [mbar]."""
    data = []
    with open(path, encoding='utf-8') as f:
        for raw in f:
            line = raw.split('!', 1)[0].strip()
            if line.startswith('Z,T,P'):
                values = line[5:].replace(',', ' ').split()
                if len(values) >= 3:
                    data.append([float(x) for x in values[:3]])
    if not data:
        raise ValueError(f'No Kaltenegger Z,T,P records in {path}')
    a = np.asarray(data)
    order = np.argsort(a[:, 2])
    return Atmosphere(a[order, 2] / 1000, a[order, 1],
                      altitude_km=a[order, 0], tp_source=str(path))


def read_kaltenegger_chemistry(path: str | Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Parse Molecule, *_MR, *_ALT blocks; return altitude and linear VMR."""
    chunks = []
    continuation = ''
    with open(path, encoding='utf-8') as handle:
        for original in handle:
            line = original.split('!', 1)[0].strip()
            if not line:
                continue
            if line.endswith('^'):
                continuation += line[:-1].strip() + ' '
                continue
            chunks.append(continuation + line)
            continuation = ''
    if continuation:
        raise ValueError('Unterminated continuation in chemistry file')
    result = {}
    name = None
    fields = {}

    def save():
        if name is None or name.upper() == 'DUMMY':
            return
        try:
            z = np.asarray(fields[f'{name}_ALT'], dtype=float)
            mr = np.asarray(fields[f'{name}_MR'], dtype=float)
        except KeyError as exc:
            raise ValueError(f'Incomplete chemistry block for {name}') from exc
        if len(z) != len(mr) or len(z) < 2 or not np.all(np.diff(z) > 0):
            raise ValueError(f'{name}: invalid altitude grid or inconsistent dimensions')
        if not np.all(np.isfinite(z)) or not np.all(np.isfinite(mr)) or np.any(mr < 0):
            raise ValueError(f'{name}: invalid abundance values')
        result[name] = (z, mr)

    for line in chunks:
        match = re.fullmatch(r'Molecule\s+(\S+)', line)
        if match:
            save()
            name = match.group(1)
            fields = {}
            if name == '&NULL':
                name = None
                break
            continue
        match = re.fullmatch(r'(\S+_(?:MR|ALT))\s+(.+)', line)
        if match:
            fields[match.group(1)] = [float(s) for s in re.split(r'[\s,]+', match.group(2)) if s]
    save()
    if not result:
        raise ValueError(f'No chemistry species in {path}')
    return result


def composition_on_atmosphere(
    chemistry: dict[str, tuple[np.ndarray, np.ndarray]], atmosphere: Atmosphere,
) -> Composition:
    """Interpolate VMR linearly in altitude, matching legacy endpoint behavior."""
    if atmosphere.altitude_km is None:
        raise ValueError('Kaltenegger chemistry requires an atmosphere altitude grid')
    output = {}
    for name, (z, mr) in chemistry.items():
        output[name] = np.interp(atmosphere.altitude_km, z, mr,
                                 left=mr[0], right=mr[-1])
    return Composition(output)


KALTENEGGER_AGES_GA = (0.0, 0.8, 2.0, 3.9)
KALTENEGGER_EARTH_REFERENCE_AGE_GYR = 4.6


def select_kaltenegger_profiles(config: Config, profile_root: str | Path):
    """Select an existing PT/chemistry pair, without modifying NASA metadata.

    NASA Gyr_age is the measured *stellar* age, not geological lookback.
    Treat it as an Earth-analogue evolutionary clock (4.6 Gyr reference),
    an explicit modeling approximation rather than a physical equivalence.
    """
    import math
    from .nasa import closest_spectral_type

    root = Path(profile_root)
    desired_type = closest_spectral_type(config.get('stellar_type'))
    try:
        stellar_age = float(config['Gyr_age'])
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError('Kaltenegger selection requires finite NASA Gyr_age in Gyr') from exc
    if not math.isfinite(stellar_age) or stellar_age < 0:
        raise ValueError('Kaltenegger selection requires nonnegative finite NASA Gyr_age')
    lookback = max(0.0, KALTENEGGER_EARTH_REFERENCE_AGE_GYR - stellar_age)

    pairs = []
    for age in KALTENEGGER_AGES_GA:
        age_label = f'{age:.1f}'
        for name in ('F0V','F7V','G8V','K2V','K7V','M1A','M3A','M8A','Sun'):
            pt = root / 'PT_profiles' / f'AltPT_{name}_{age_label}Ga.txt'
            chem = root / 'chem_profiles' / f'Chem_{name}_{age_label}Ga.txt'
            if pt.is_file() and chem.is_file():
                pairs.append((name, age, pt, chem))
    if not pairs:
        raise FileNotFoundError(f'No matching Kaltenegger PT and chemistry pairs under {root}')

    # Keep the stellar-template preference first. Within the selected type,
    # choose the nearest geological age with both files present.
    available_types = {name for name, _, _, _ in pairs}
    if desired_type in available_types:
        selected_type = desired_type
    else:
        # Use the same ordered spectral coordinate as the NASA mapper,
        # explicitly treating Sun as G2V rather than silently picking G8V.
        coord = {'F0V':0., 'F7V':7., 'Sun':12., 'G8V':18.,
                 'K2V':22., 'K7V':27., 'M1A':31., 'M3A':33., 'M8A':38.}
        selected_type = min(available_types, key=lambda x:(abs(coord[x]-coord[desired_type]), coord[x]))
        import warnings
        warnings.warn(f'Kaltenegger template {desired_type} is unavailable; using {selected_type}',
                      UserWarning, stacklevel=2)
    name, age, pt, chem = min((x for x in pairs if x[0] == selected_type),
                              key=lambda x:(abs(x[1]-lookback), x[1]))
    print(f'[KALTENEGGER] NASA stellar_type={config.get("stellar_type")!r}, '
          f'Gyr_age={stellar_age:g} Gyr; requested lookback={lookback:.3f} Ga '
          f'(Earth-analogue reference {KALTENEGGER_EARTH_REFERENCE_AGE_GYR:g} Gyr).')
    print(f'[KALTENEGGER] Selected {name}, {age:.1f} Ga: {pt.name} / {chem.name}')
    return pt, chem


def select_legacy_profiles(config: Config, profile_root: str | Path):
    """Compatibility alias for the modern automatic Kaltenegger selector."""
    return select_kaltenegger_profiles(config, profile_root)


@dataclass(frozen=True)
class TransmissionSpectrum:
    wavelength_um: np.ndarray
    transmission: np.ndarray

    @property
    def transit_depth(self):
        """Deprecated compatibility alias; values are transmission (1 - raw POSEIDON output)."""
        return self.transmission

    def plot(self, path: str | Path) -> Path:
        """Plot transmission without saving numerical POSEIDON products."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(11, 4))
        ax.plot(self.wavelength_um, self.transmission, lw=0.7)
        ax.set_xlabel('Wavelength (µm)')
        ax.set_ylabel('Transmission T')
        ax.set_title('POSEIDON transmission spectrum')
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    def __post_init__(self):
        wl = np.asarray(self.wavelength_um, dtype=float)
        depth = np.asarray(self.transmission, dtype=float)
        if wl.ndim != 1 or depth.shape != wl.shape or not np.all(np.isfinite(wl)):
            raise ValueError('Spectrum arrays must be finite-wavelength matching 1-D arrays')
        if not np.all(np.isfinite(depth)):
            raise ValueError('Non-finite transmission spectrum')




def validate_hminus_continuum(model: dict, species) -> None:
    """Verify POSEIDON registered the continuum implied by supplied chemistry.

    POSEIDON/core.py derives H-bf from H- and H-ff from H plus e-.
    Neither electron abundance nor H- should be silently substituted for the other.
    """
    present = set(species)
    bf = set(model.get('bf_species', ()))
    ff = set(model.get('ff_pairs', ()))
    if 'H-' in present and 'H-bf' not in bf:
        raise RuntimeError('POSEIDON did not register H-bf for H- chemistry')
    if {'H', 'e-'} <= present and 'H-ff' not in ff:
        raise RuntimeError('POSEIDON did not register H-ff for H + e- chemistry')
    if ('H-' in present or 'e-' in present) and not {'H', 'H-', 'e-'} <= present:
        raise ValueError('H-minus continuum experiment requires H, H-, and e- abundances explicitly')

@dataclass
class PoseidonTransmissionModel:
    """Backend adapter, with no POSEIDON import at construction time."""
    resolution: float = 300_000
    wavelength_min_um: float = 0.6
    wavelength_max_um: float = 2.5

    @classmethod
    def from_config(cls, config: Config, *, config_path: str | Path = "config"):
        if 'instrument' in config or 'instrument_order_file' in config:
            from .instruments import instrument_wavelength_bounds_um
            lower, upper = instrument_wavelength_bounds_um(config, config_path)
        else:
            lower, upper = config.get_float('wl_min', 0.6), config.get_float('wl_max', 2.5)
        model = cls(config.get_float('resolution', 300_000), lower, upper)
        if not (model.resolution > 0 and 0 < model.wavelength_min_um < model.wavelength_max_um):
            raise ValueError('Invalid spectral resolution or wavelength bounds')
        return model

    def validate(self, atmosphere: Atmosphere, composition: Composition,
                 bulk_species: str | tuple[str, ...], planet: Planet, star: Star):
        if atmosphere.altitude_km is None:
            # Not intrinsically required by POSEIDON, but useful for the legacy profile setup.
            pass
        bulk = (bulk_species,) if isinstance(bulk_species, str) else tuple(bulk_species)
        if not bulk or len(bulk) != len(set(bulk)):
            raise ValueError('Specify at least one unique bulk species')
        missing = set(bulk) - set(composition.mixing_ratios)
        if missing:
            raise ValueError(f'Bulk species missing from chemistry: {sorted(missing)}')
        unsupported = set(composition.mixing_ratios) - SUPPORTED_SPECIES
        if unsupported:
            raise ValueError(f'Species unsupported by the current POSEIDON adapter: {sorted(unsupported)}')
        if any(len(x) != len(atmosphere.pressure_bar) for x in composition.mixing_ratios.values()):
            raise ValueError('Composition layer count differs from T-P layer count')
        if planet.radius_jupiter is None or planet.mass_jupiter is None or planet.equilibrium_temperature_k is None:
            raise ValueError('POSEIDON requires Rp, Mp, Teq')
        if star.radius_solar is None or star.teff_k is None or star.logg_cgs is None:
            raise ValueError('POSEIDON requires Rs, Teff and stellar logg')

    def compute_transmission(self, *, name: str, config: Config, atmosphere: Atmosphere,
                             composition: Composition, bulk_species: str | tuple[str, ...]) -> TransmissionSpectrum:
        """Run legacy POSEIDON transmission calculation (expensive; unverified backend integration)."""
        star = Star.from_config(config)
        planet = Planet.from_config(config)
        self.validate(atmosphere, composition, bulk_species, planet, star)
        configure_poseidon_environment(config)
        configure_backend_path(config, key='poseidon_path', package='POSEIDON')
        try:
            from POSEIDON.constants import R_Sun, R_J, M_J
            from POSEIDON.core import (create_star, create_planet, define_model,
                                       make_atmosphere, wl_grid_constant_R,
                                       read_opacities, compute_spectrum)
        except ImportError as exc:
            raise ImportError('POSEIDON is required to compute spectra; validation works without it') from exc
        poseidon_star = create_star(star.radius_solar * R_Sun, star.teff_k,
                                    star.logg_cgs, config.get_float('metallicity', 0.0))
        poseidon_planet = create_planet(name, planet.radius_jupiter * R_J,
                                        mass=planet.mass_jupiter * M_J,
                                        T_eq=planet.equilibrium_temperature_k)
        bulk = (bulk_species,) if isinstance(bulk_species, str) else tuple(bulk_species)
        param_species = [s for s in composition.mixing_ratios if s not in bulk]
        model = define_model(name, list(bulk), param_species, PT_profile='file_read',
                             X_profile='file_read', cloud_model='cloud-free', radius_unit='R_J')
        validate_hminus_continuum(model, composition.mixing_ratios)
        order = param_species + list(bulk)
        # The legacy script reverses ascending-pressure arrays before make_atmosphere.
        p = atmosphere.pressure_bar[::-1]
        t = atmosphere.temperature_k[::-1]
        x = np.stack([composition.mixing_ratios[s][::-1] for s in order])
        p_ref = config.get_float('P_ref_bar', 1.0)
        poseidon_atmosphere = make_atmosphere(poseidon_planet, model, p,
                                              p_ref, planet.radius_jupiter * R_J,
                                              T_input=t, X_input=x)
        wl = wl_grid_constant_R(self.wavelength_min_um, self.wavelength_max_um, self.resolution)
        # Keep the sampling policy from the legacy script, including the
        # configured opacity grid bounds. Review against installed POSEIDON version.
        t_min = max(float(np.min(t)) - 100.0, 0.0)
        t_max = float(np.max(t)) + 200.0
        fine_t = np.arange(t_min, t_max + 20.0, 20.0)
        # Sample the full atmospheric pressure interval. POSEIDON's underlying
        # cross-section tables have their own native bounds; outside those
        # bounds the backend applies its own boundary policy. This does not
        # establish physical accuracy of line opacities at 1e-15 bar.
        log_p_min = float(np.floor(np.log10(np.min(p))))
        log_p_max = float(np.ceil(np.log10(np.max(p))))
        fine_log_p = np.linspace(log_p_min, log_p_max,
                                 int(np.ceil((log_p_max-log_p_min)/0.2)) + 1)
        opacities = read_opacities(model, wl, 'opacity_sampling', fine_t, fine_log_p)
        spectrum = compute_spectrum(poseidon_planet, poseidon_star, model,
                                     poseidon_atmosphere, opacities, wl,
                                     spectrum_type='transmission', save_spectrum=False)
        # Legacy run_poseidon.py saved (1-spectrum); its exact scientific meaning
        # must be checked against the user's installed POSEIDON version.
        return TransmissionSpectrum(np.asarray(wl), 1.0 - np.asarray(spectrum).squeeze())


@dataclass(frozen=True)
class TransmissionRun:
    """Prepared standard inputs; validate/dry-run without importing POSEIDON."""
    config: Config
    name: str
    atmosphere: Atmosphere
    abundances: 'AbundanceProfile'
    model: PoseidonTransmissionModel

    @classmethod
    def from_config(cls, path: str | Path) -> 'TransmissionRun':
        """Read familiar key/value config plus standard file chemistry.

        Required: bulk_species, chemistry=file, abundances_file. T-P uses the
        already-supported tp_profile config choice (file/isothermal/guillot).
        Paths are relative to the config directory. Never modify the config.
        """
        from .profiles import AbundanceProfile
        path = Path(path)
        config = Config.from_file(path)
        chemistry = config['chemistry'].strip().lower() if 'chemistry' in config else 'file'
        if chemistry == 'kaltenegger':
            from .converters.kaltenegger import KalteneggerConverter
            if 'kaltenegger_profiles_dir' in config:
                root = Path(config['kaltenegger_profiles_dir']).expanduser()
                if not root.is_absolute():
                    root = path.parent / root
            else:
                from .backend_paths import discover_site_config
                site_path = discover_site_config()
                site = Config.from_file(site_path) if site_path is not None else None
                if site is None or 'kaltenegger_profiles_dir' not in site:
                    raise FileNotFoundError(
                        'Kaltenegger data are not bundled. Download from '
                        'https://zenodo.org/records/14829587 and set '
                        'kaltenegger_profiles_dir in config/site_paths.conf. '
                        'The directory must contain PT_profiles/ and chem_profiles/.'
                    )
                root = Path(site['kaltenegger_profiles_dir']).expanduser()
                if not root.is_absolute():
                    root = site_path.parent / root
            root = root.resolve()
            if not (root / 'PT_profiles').is_dir() or not (root / 'chem_profiles').is_dir():
                raise FileNotFoundError(
                    f'Invalid Kaltenegger profiles directory {root}: expected '
                    'PT_profiles/ and chem_profiles/ subdirectories'
                )
            pt_path, chem_path = select_legacy_profiles(config, root)
            atmosphere, abundances = KalteneggerConverter().load_earth_like(pt_path, chem_path)
        else:
            if 'bulk_species' not in config or not config['bulk_species'].strip():
                raise ValueError('Config must explicitly specify bulk_species')
            atmosphere = Atmosphere.from_config(config, config_dir=path.parent)
        if chemistry == 'kaltenegger':
            pass
        elif chemistry == 'fastchem':
            from .fastchem import equilibrium_profile
            abundances = equilibrium_profile(config, atmosphere, config_dir=path.parent)
        elif chemistry == 'manual':
            from .manual_chemistry import manual_profile
            abundances = manual_profile(config, atmosphere)
        elif chemistry == 'file':
            if 'abundances_file' not in config:
                raise ValueError('Config requires abundances_file for chemistry=file')
            abundance_path = Path(config['abundances_file']).expanduser()
            if not abundance_path.is_absolute():
                abundance_path = path.parent / abundance_path
            abundances = AbundanceProfile.from_file(abundance_path, bulk_species=config['bulk_species'])
        else:
            raise ValueError('chemistry must be kaltenegger, file, manual, or fastchem')
        abundances = abundances.on_atmosphere(atmosphere)
        model = PoseidonTransmissionModel.from_config(config, config_path=path)
        model.validate(atmosphere, abundances.as_poseidon_composition(),
                       abundances.bulk_species, Planet.from_config(config), Star.from_config(config))
        return cls(config, path.parent.name, atmosphere, abundances, model)

    def compute(self) -> TransmissionSpectrum:
        return self.model.compute_transmission(
            name=self.name, config=self.config, atmosphere=self.atmosphere,
            composition=self.abundances.as_poseidon_composition(),
            bulk_species=self.abundances.bulk_species,
        )

    def save(self, spectrum: TransmissionSpectrum, *, directory: str | Path) -> tuple[Path, Path]:
        """Save arrays without choosing a legacy-compatible FITS schema prematurely."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        wavelength_path = directory / 'wavelength_um.npy'
        depth_path = directory / 'transmission.npy'
        np.save(wavelength_path, spectrum.wavelength_um)
        np.save(depth_path, spectrum.transmission)
        return wavelength_path, depth_path
