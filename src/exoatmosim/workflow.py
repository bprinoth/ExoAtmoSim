"""Stable public Python entry points for the modular simulation stages.

Each stage uses the existing, independently tested scientific adapters. A full
simulation delegates to the same pipeline as the command-line interface.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

from .config import Config
from .full_crires import resolve_config, run as _simulate
from .plotting import publication_style


def make_config(planet: str, *, instrument: str = 'CRIRES_K',
                data_dir: str | Path = 'data', site_config: str | Path | None = None,
                refresh: bool = False, fresh: bool = False) -> Path:
    """Create or reuse a planet template (NASA is queried only if needed)."""
    from .backend_paths import discover_site_config
    return resolve_config(None, planet=planet, site_config=discover_site_config(site_config),
                          data_dir=Path(data_dir), refresh_config=refresh, fresh=fresh,
                          instrument=instrument.upper())


def _locations(config: str | Path, seed: int):
    path = Path(config).resolve()
    if not path.is_file():
        raise FileNotFoundError(f'Planet configuration not found: {path}')
    cfg = Config.from_file(path)
    instrument = str(cfg.get('instrument', '')).upper()
    if path.parent.name.startswith('seed_') and path.parent.parent.name == instrument:
        root = path.parent
        planet_name = path.parent.parent.parent.parent.name
    else:
        root = path.parent / 'simulations' / instrument / f'seed_{seed:04d}'
        planet_name = path.parent.name
    return path, root, planet_name


def _snapshot(config: Path, run_dir: Path, *, overwrite: bool):
    run_dir.mkdir(parents=True, exist_ok=True)
    target = run_dir / 'config'
    contents = config.read_bytes()
    if target.exists() and target.read_bytes() != contents and not overwrite:
        raise FileExistsError(f'Different configuration already saved: {target}')
    if target != config:
        target.write_bytes(contents)
    return target


def calculate_atmosphere(config: str | Path, *, seed: int = 42,
                         overwrite: bool = False) -> Path:
    """Compute T-P and chemistry diagnostics, without performing radiative transfer.

    Returns the diagnostics directory. The FastChem and profile objects are
    produced by TransmissionRun's tested setup code.
    """
    from .poseidon import TransmissionRun
    path, run_dir, _ = _locations(config, seed)
    _snapshot(path, run_dir, overwrite=overwrite)
    diagnostics = run_dir / 'diagnostics'
    diagnostics.mkdir(parents=True, exist_ok=True)
    run = TransmissionRun.from_config(path)
    with publication_style():
        run.atmosphere.plot_tp(diagnostics / 'temperature_pressure.png')
        run.abundances.plot(diagnostics / f'{str(run.config.get("chemistry", "fastchem")).lower()}_abundances.png',
                            title='Atmospheric abundances')
    run.abundances.to_file(diagnostics / f'{str(run.config.get("chemistry", "fastchem")).lower()}_abundances.dat')
    return diagnostics


def calculate_spectrum(config: str | Path, *, backend: str = 'poseidon',
                       seed: int = 42, overwrite: bool = False,
                       create_templates: bool = False) -> Path:
    """Compute a new transmission model; currently only POSEIDON is supported.

    The spectrum is saved to the same seed's instrument directory so observe
    can consume it without a re-computation. Does not call the observing backend.
    """
    if backend.lower() != 'poseidon':
        raise ValueError('Only POSEIDON is implemented; petitRADTRANS is planned')
    from .poseidon import TransmissionRun
    from .spectrum_fits import write_spectrum_fits
    path, run_dir, planet = _locations(config, seed)
    model_path = run_dir / f'{planet}.fits'
    if model_path.exists() and not overwrite:
        raise FileExistsError(f'Existing planetary model: {model_path}; use --overwrite')
    _snapshot(path, run_dir, overwrite=overwrite)
    diagnostics = run_dir / 'diagnostics'
    diagnostics.mkdir(parents=True, exist_ok=True)
    run = TransmissionRun.from_config(path)
    species = tuple(name for name in run.abundances.mixing_ratios
                    if name not in run.abundances.bulk_species) if create_templates else ()
    if create_templates:
        existing = [run_dir / 'templates' / f'{name}.fits' for name in species]
        if not overwrite and any(path.exists() for path in existing):
            raise FileExistsError('Some molecular templates already exist; use --overwrite')
        spectrum, templates = run.model.compute_transmission(
            name=run.name, config=run.config, atmosphere=run.atmosphere,
            composition=run.abundances.as_poseidon_composition(),
            bulk_species=run.abundances.bulk_species, template_species=species)
    else:
        spectrum = run.compute()
        templates = {}
    transmission = np.asarray(spectrum.transmission, float)
    if not np.isfinite(transmission).all() or not np.any(np.abs(transmission - 1) > 1e-12):
        raise ValueError('Planetary transmission must be finite and non-unity')
    with publication_style():
        spectrum.plot(diagnostics / 'poseidon_transmission.png')
    run.save(spectrum, directory=diagnostics / 'poseidon_arrays')
    # Validate all contributions before writing any new FITS products.
    depths = {}
    if create_templates:
        from .poseidon import TransmissionSpectrum
        if set(templates) != set(species):
            raise RuntimeError('POSEIDON returned an incomplete set of molecular templates')
        for name, template in templates.items():
            depth = TransmissionSpectrum(template.wavelength_um,
                                         1.0 - np.asarray(template.transmission))
            values = np.asarray(depth.transmission, float)
            amplitude = float(np.ptp(values))
            if amplitude <= 1e-12:
                raise ValueError(
                    f'{name} opacity contribution is flat (peak-to-peak '
                    f'{amplitude:.3g}); check abundance and POSEIDON '
                    'spectral_contribution before using the template')
            depths[name] = depth
    result = write_spectrum_fits(spectrum, model_path, overwrite=overwrite,
                                 config_path=path)
    for name, depth in depths.items():
        destination = run_dir / 'templates' / f'{name}.fits'
        write_spectrum_fits(depth, destination, overwrite=overwrite,
                            config_path=path)
        with publication_style():
            depth.plot(destination.with_suffix('.png'))
    if create_templates:
        print(f'Created {len(depths)} fixed-structure molecular templates: '
              + ', '.join(depths))
    return result


def calculate_template(config: str | Path, *, species: str, seed: int = 42,
                       overwrite: bool = False) -> Path:
    """Compute a species-only line-opacity transmission template using POSEIDON.

    Compute the full-composition atmospheric structure (including its original
    mean molecular weight), and isolate only the requested gas opacity.
    Uses the same wavelength grid and FITS schema as `spectrum`, but
    stores POSEIDON transit depth T rather than the usual 1 - T.
    """
    from .poseidon import TransmissionRun, TransmissionSpectrum
    from .spectrum_fits import write_spectrum_fits

    species = species.strip()
    if not species or '/' in species or chr(92) in species or species in ('.', '..'):
        raise ValueError('Provide a valid single species name (e.g. O3)')
    path, run_dir, _ = _locations(config, seed)
    destination = run_dir / 'templates' / f'{species}.fits'
    if destination.exists() and not overwrite:
        raise FileExistsError(f'Existing species template: {destination}; use --overwrite')
    run = TransmissionRun.from_config(path)
    ratios = run.abundances.mixing_ratios
    if species not in ratios:
        raise ValueError(f'Species {species!r} is not in the atmospheric abundances; '
                         f'available: {", ".join(sorted(ratios))}')
    bulk = run.abundances.bulk_species
    if species in bulk:
        raise ValueError(f'{species!r} is a bulk species; choose a trace absorber')
    # The complete composition must be retained when POSEIDON constructs
    # pressure-radius and mean-molecular-weight profiles. Opacity filtering
    # occurs *after* make_atmosphere via POSEIDON's contribution calculation.
    _snapshot(path, run_dir, overwrite=overwrite)
    spectrum = run.model.compute_transmission(
        name=run.name, config=run.config, atmosphere=run.atmosphere,
        composition=run.abundances.as_poseidon_composition(),
        bulk_species=bulk, opacity_species=species)
    transmission = np.asarray(spectrum.transmission, float)
    if not np.isfinite(transmission).all():
        raise ValueError('Species template transmission must be finite')
    # The POSEIDON adapter returns 1 - T for the standard planet model.
    # Templates require the un-inverted transit depth T instead.
    transit_depth = TransmissionSpectrum(spectrum.wavelength_um, 1.0 - transmission)
    if np.ptp(transit_depth.transmission) <= 1e-12:
        raise ValueError(f'{species} opacity contribution is flat; cannot use as CCF template')
    result = write_spectrum_fits(transit_depth, destination, overwrite=overwrite,
                                 config_path=path)
    with publication_style():
        transit_depth.plot(destination.with_suffix('.png'),
                           title=f'{species} template transit depth',
                           ylabel='Transit depth (T)')
    return result


def simulate_observation(config: str | Path, *, seed: int = 42, snr: float | None = None,
                         overwrite: bool = False, stellar_npz: str | Path | None = None):
    """Synthesize noisy orders from an EXISTING seed-scoped planetary FITS.

    Missing planetary models are an error: use `spectrum` or `simulate` first.
    """
    from .independent_crires import run as observe_run
    path, run_dir, planet = _locations(config, seed)
    model_path = run_dir / f'{planet}.fits'
    if not model_path.is_file():
        raise FileNotFoundError(f'Planet model not found: {model_path}; run `exoatmosim spectrum` first')
    instrument = str(Config.from_file(path).get('instrument', '')).upper()
    if not instrument.startswith('VIOLA_') and snr is None:
        snr = 100.0
    return observe_run(path, Path(stellar_npz) if stellar_npz else None, output=run_dir,
                       snr=snr, seed=seed, overwrite=overwrite, model_path=model_path)


def simulate(*, planet: str | None = None, config: str | Path | None = None,
             instrument: str = 'CRIRES_K', data_dir: str | Path = 'data',
             site_config: str | Path | None = None, snr: float | None = None, seed: int = 42,
             tellurics: bool | None = None, overwrite: bool = False,
             fresh: bool = False):
    """Run the existing end-to-end validated orchestration.

    Pass either planet= (NASA-backed config) or config= (existing template).
    """
    if (planet is None) == (config is None):
        raise ValueError('Pass exactly one of planet= or config=')
    path = (make_config(planet, instrument=instrument, data_dir=data_dir,
                        site_config=site_config, fresh=fresh) if planet is not None
            else Path(config).resolve())
    return _simulate(path, snr=snr, seed=seed, overwrite=overwrite, tellurics=tellurics)
