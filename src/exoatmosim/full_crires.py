"""Full FastChem -> POSEIDON -> PHOENIX -> instrument simulation.

Existing target:
    exoatmosim simulate data/TOI-1842_b/config --snr 100 --overwrite

New target (needs network and configured scientific backends):
    exoatmosim simulate --planet 'TOI-1842 b' \\
        --site-config config/site_paths.conf --snr 100

No cross-correlation or external mock package is invoked. Planetary FITS
export and PHOENIX retrieval are implemented within ExoAtmoSim.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import os
import tempfile
from pathlib import Path

from exoatmosim.config import Config
from exoatmosim.plotting import publication_style


def resolve_config(config_path: Path | None, *, planet: str | None,
                   site_config: Path | None, refresh_config: bool = False,
                   data_dir: Path = Path('data'), fresh: bool = False,
                   instrument: str | None = None) -> Path:
    """Reuse the user's science config unless explicitly asked for NASA creation."""
    supported = ('CRIRES_K', 'ANDES', 'ANDES_K', 'IGRINS_HK', 'IGRINS2_HK', 'SPIROU', 'ESPRESSO', 'VIOLA_J', 'VIOLA_H', 'VIOLA_K')
    if instrument is not None and instrument.upper() not in supported:
        raise ValueError(f'Unsupported instrument {instrument!r}; choose from {supported}')
    if config_path is not None and planet is not None:
        raise ValueError('Provide an existing config OR --planet, not both')
    if config_path is None and planet is None:
        raise ValueError('Provide a config path or --planet for NASA config creation')
    if planet is not None:
        config_path = Path(data_dir) / planet.replace(' ', '_') / 'config'
        if fresh and config_path.exists():
            raise FileExistsError(f'Fresh NASA run refuses existing {config_path}; choose an empty --data-dir')
        if config_path.exists() and not refresh_config:
            if instrument is not None and str(Config.from_file(config_path).get('instrument', '')).upper() != instrument:
                raise ValueError(f'Existing config {config_path} has a different instrument. '
                                 'Use an instrument-specific config or a separate data root.')
            print(f'[CONFIG] Reusing {config_path} (NASA query skipped).')
            return config_path.resolve()
        from exoatmosim.nasa import find_parameters
        print('[1/5] NASA: retrieving planet parameters and creating config')
        find_parameters(planet, data_dir=data_dir, instrument=instrument or 'CRIRES_K',
                        site_config=site_config, print_params=False,
                        overwrite=refresh_config, preserve_existing=True)
        # Newly generated NASA config: select the intended chemistry explicitly.
        generated = Config.from_file(config_path)
        generated.set('chemistry', 'fastchem')
        generated.set('bulk_species', 'H2 He')
        generated.set('fastchem_species', 'H2 He Fe' if instrument == 'ESPRESSO' else 'H2 He H2O')
        generated.set('tellurics', 'False')
        generated.write(config_path)
        return config_path.resolve()
    config_path = Path(config_path).resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f'Config not found: {config_path}. '
                                'Use --planet to create one via NASA.')
    if refresh_config:
        raise ValueError('--refresh-config requires --planet')
    if fresh:
        raise ValueError('--fresh requires --planet')
    return config_path


def run(config_path: Path, *, snr: float | None = None, seed: int = 42,
        overwrite: bool = False, tellurics: bool | None = None):
    """Run with an optional seed-only telluric override.

    Keep the planet-level NASA template untouched. A temporary sibling config
    ensures legacy relative paths still resolve against the planet directory.
    The effective config is snapshotted into the seed folder by the runner.
    """
    config_path = Path(config_path).resolve()
    if tellurics is None:
        return _run_pipeline(config_path, snr=snr, seed=seed, overwrite=overwrite)

    cfg = Config.from_file(config_path)
    cfg.set('tellurics', bool(tellurics))
    fd, temporary = tempfile.mkstemp(prefix='.exoatmosim_effective_',
                                    suffix='.conf', dir=config_path.parent)
    os.close(fd)
    try:
        effective = Path(temporary)
        cfg.write(effective)
        print(f'[CONFIG] Tellurics override: {bool(tellurics)} (planet template unchanged)')
        return _run_pipeline(effective, snr=snr, seed=seed, overwrite=overwrite)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _run_pipeline(config_path: Path, *, snr: float | None = None, seed: int = 42,
                  overwrite: bool = False):
    import numpy as np
    from exoatmosim.poseidon import TransmissionRun
    from exoatmosim.spectrum_fits import write_spectrum_fits
    from exoatmosim.independent_crires import run as run_crires

    config_path = Path(config_path).resolve()
    cfg = Config.from_file(config_path)
    instrument = str(cfg.get('instrument', '')).upper()
    if instrument not in ('CRIRES_K', 'ANDES', 'ANDES_K', 'IGRINS_HK', 'IGRINS2_HK', 'SPIROU', 'ESPRESSO', 'VIOLA_J', 'VIOLA_H', 'VIOLA_K'):
        raise ValueError(f'Unsupported full-pipeline instrument: {instrument}')
    chemistry = str(cfg.get('chemistry', 'fastchem')).lower()
    if chemistry not in ('fastchem', 'kaltenegger'):
        raise ValueError('Full pipeline requires chemistry fastchem or kaltenegger')
    if instrument.startswith('VIOLA_'):
        if snr is not None:
            raise ValueError('VIOLA computes ETC electron noise: do not supply --snr')
    else:
        if snr is None: snr = 100.
        if not np.isfinite(snr) or snr <= 0:
            raise ValueError('--snr must be positive and finite')
    mock_path = config_path.parent / 'simulations' / instrument / f'seed_{seed:04d}'
    model_path = mock_path / f'{config_path.parent.name}.fits'
    snapshot_path = mock_path / 'config'
    incoming_config = config_path.read_bytes()
    # A seed is an immutable input snapshot unless overwrite is explicit.
    if snapshot_path.exists() and snapshot_path.read_bytes() != incoming_config and not overwrite:
        raise FileExistsError(
            f'Seed {seed} already has a different configuration: {snapshot_path}. '
            'Choose a new --seed or use --overwrite deliberately.')
    if not overwrite:
        existing = [p for p in [model_path, mock_path / 'order_0.fits'] if p.exists()]
        if existing:
            raise FileExistsError(f'Outputs already exist: {existing}. Use --overwrite to replace them.')

    diagnostic = mock_path / 'diagnostics'
    print('\n' + '=' * 64)
    print(f' ExoAtmoSim | {chemistry} / POSEIDON / PHOENIX / {instrument}')
    print('=' * 64)
    print(f'Planet: {config_path.parent.name}  |  Seed: {seed}  |  '
          + ('VIOLA ETC electron noise' if instrument.startswith('VIOLA_') else f'SNR: {snr:g}'))
    print(f'Output: {mock_path}')
    mock_path.mkdir(parents=True, exist_ok=True)
    if not snapshot_path.exists() or snapshot_path.read_bytes() != incoming_config:
        snapshot_path.write_bytes(incoming_config)
    print(f'Config snapshot: {snapshot_path}')
    diagnostic.mkdir(parents=True, exist_ok=True)
    print(f'[2/5] {chemistry}: chemistry and atmospheric structure')
    atmosphere_run = TransmissionRun.from_config(config_path)
    with publication_style():
        atmosphere_run.atmosphere.plot_tp(diagnostic / 'temperature_pressure.png')
        atmosphere_run.abundances.plot(diagnostic / f'{chemistry}_abundances.png',
                                       title=f'{chemistry} abundances')
    atmosphere_run.abundances.to_file(diagnostic / f'{chemistry}_abundances.dat')
    print('[3/5] POSEIDON: calculating planetary transmission')
    spectrum = atmosphere_run.compute()
    t = np.asarray(spectrum.transmission, float)
    if not np.isfinite(t).all() or not np.any(np.abs(t - 1.) > 1e-12):
        raise ValueError('Planet model has no finite, non-unity transmission signal')
    print(f'      Planet signal: min(T)={t.min():.8f}, max |1-T|={np.max(np.abs(1.-t)):.3e}')
    with publication_style():
        spectrum.plot(diagnostic / 'poseidon_transmission.png')
    atmosphere_run.save(spectrum, directory=diagnostic / 'poseidon_arrays')
    print('[4/5] Exporting planetary model to air-wavelength FITS')
    write_spectrum_fits(spectrum, model_path, overwrite=overwrite,
                              config_path=config_path)
    print(f'[5/5] Generating schedule, stellar spectrum and {instrument} exposures')
    written = run_crires(config_path, None, output=mock_path, snr=snr, seed=seed,
                         overwrite=overwrite, model_path=model_path)
    from exoatmosim.instruments import instrument_resource
    from exoatmosim.instruments import load_instrument_orders
    expected = len(load_instrument_orders(instrument))
    if len(written) != 2 * expected:
        raise RuntimeError(f'Expected {expected} order FITS and {expected} wavelength FITS; '
                           f'writer reported {len(written)} files')
    for i in range(expected):
        for prefix in ('order', 'wave'):
            f = mock_path / f'{prefix}_{i}.fits'
            if not f.is_file() or f.stat().st_size == 0:
                raise RuntimeError(f'Missing or empty output {f}')
    if not (mock_path / 'obs_times').is_file():
        raise RuntimeError('No obs_times file produced')
    print(f'[SUCCESS] {expected} {instrument} spectral orders and their wavelengths validated.')
    print(f'[SUCCESS] Observations: {mock_path}')
    print(f'[SUCCESS] Atmosphere diagnostics: {diagnostic}')
    return written


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', type=Path,
                        help='Existing planet config; omit when using --planet')
    parser.add_argument('--planet', help='NASA Exoplanet Archive name (create config if absent)')
    parser.add_argument('--site-config', type=Path, default=None,
                        help='Local backend paths (default: config/site_paths.conf or EXOATMOSIM_SITE_CONFIG)')
    parser.add_argument('--refresh-config', action='store_true',
                        help='Explicitly refresh NASA keys while preserving other config values')
    parser.add_argument('--instrument', type=str.upper, default=None,
                        help='CRIRES_K, ANDES (experimental), ANDES_K (provisional), IGRINS_HK, IGRINS2_HK, SPIROU, ESPRESSO (experimental)' )
    telluric_group = parser.add_mutually_exclusive_group()
    telluric_group.add_argument('--tellurics', dest='tellurics', action='store_true',
                               help='Enable TelFit for this seed, without editing planet config')
    telluric_group.add_argument('--no-tellurics', dest='tellurics', action='store_false',
                               help='Disable TelFit for this seed, without editing planet config')
    parser.set_defaults(tellurics=None)
    parser.add_argument('--snr', type=float, default=None,
                        help='Constant SNR for non-VIOLA; VIOLA derives ETC noise')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--data-dir', type=Path, default=Path('data'),
                        help='NASA generated files and outputs root (default: data)')
    parser.add_argument('--fresh', action='store_true',
                        help='Require a newly created NASA config; never reuse an existing one')
    parser.add_argument('--overwrite', action='store_true',
                        help='Explicitly replace existing seed-specific FITS products')
    args = parser.parse_args(argv)
    from exoatmosim.backend_paths import discover_site_config
    site_config = discover_site_config(args.site_config)
    if site_config is not None:
        from exoatmosim.nasa import load_site_config
        from exoatmosim.backend_paths import configure_backend_path
        site = load_site_config(site_config)
    config_path = resolve_config(args.config, planet=args.planet,
                                 site_config=site_config,
                                 refresh_config=args.refresh_config,
                                 data_dir=args.data_dir, fresh=args.fresh, instrument=args.instrument)
    run(config_path, snr=args.snr, seed=args.seed, overwrite=args.overwrite,
        tellurics=args.tellurics)


if __name__ == '__main__':
    main()
