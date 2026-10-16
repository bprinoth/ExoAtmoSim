"""Prepare legacy tayph.mock inputs from one ExoAtmoSim planet config.

The model FITS schema is based on a user-provided working tayph input.
This module stages files, exports FITS, and optionally calls existing tayph.mock.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

from .config import Config
from .backend_paths import configure_backend_path
from .instruments import instrument_resource, canonical_instrument


@dataclass(frozen=True)
class TayphPreparation:
    planet_dir: Path
    runfile: Path
    order_definition: Path
    blaze: Path | None
    expected_mock_fits: Path


def _stage_resource(source: Path, destination: Path) -> None:
    if destination.exists():
        if not destination.is_file() or destination.read_bytes() != source.read_bytes():
            raise FileExistsError(
                f"Refusing to overwrite existing instrument resource: {destination}"
            )
        return
    shutil.copy2(source, destination)


def prepare_tayph_mock(config_path: str | Path, *, project_root: str | Path = '.') -> TayphPreparation:
    """Prepare instrument files and a 4-key runfile for tayph on one planet.

    Planet config is never edited. The runfile is generated under the planet
    directory, and the planet path is expressed relative to project_root to
    match legacy tayph's expectation that the process runs from project root.
    """
    root = Path(project_root).resolve()
    config_path = Path(config_path).resolve()
    cfg = Config.from_file(config_path)
    planet_dir = config_path.parent
    try:
        relative = planet_dir.relative_to(root)
    except ValueError as exc:
        raise ValueError('Planet directory must lie inside project_root for tayph relative paths') from exc
    instrument = canonical_instrument(cfg.get('instrument', 'ESPRESSO'))
    # Missing blaze calibration is delegated to tayph.mock; none is invented here.
    # Optional user-supplied paths resolve relative to the planet config directory.
    def explicit_resource(key: str, kind: str) -> Path:
        if key in cfg:
            path = Path(cfg[key]).expanduser()
            path = (path if path.is_absolute() else planet_dir / path).resolve()
            if not path.is_file():
                raise FileNotFoundError(f'{key} does not exist: {path}')
            return path
        return instrument_resource(instrument, kind)

    source_orders = explicit_resource('instrument_order_file', 'orders')
    try:
        source_blaze = explicit_resource('instrument_blaze_file', 'blaze')
    except FileNotFoundError:
        if 'instrument_blaze_file' in cfg:
            raise  # An explicitly specified resource must exist.
        source_blaze = None
    # Validate supplied files without inventing a synthetic blaze.
    import numpy as np
    if source_orders.suffix != '.npy' or (source_blaze is not None and source_blaze.suffix != '.npy'):
        raise ValueError('tayph instrument_order_file and instrument_blaze_file must be .npy arrays')
    order_grid = np.load(source_orders, allow_pickle=False)
    if order_grid.ndim != 2 or order_grid.shape[0] < 1 or order_grid.shape[1] < 2:
        raise ValueError(f'Invalid {instrument} order definition shape {order_grid.shape}')
    if not np.any(np.isfinite(order_grid) & (order_grid > 0)):
        raise ValueError(f'{instrument} order definition has no finite positive wavelengths')
    if source_blaze is not None:
        blaze_grid = np.load(source_blaze, allow_pickle=False)
        if blaze_grid.shape != order_grid.shape:
            raise ValueError(f'{instrument} blaze shape {blaze_grid.shape} does not match order shape {order_grid.shape}')
        if not np.any(np.isfinite(blaze_grid)):
            raise ValueError(f'{instrument} blaze contains no finite values')
    destination_orders = planet_dir / f'order_def_{instrument}.npy'
    destination_blaze = planet_dir / f'blaze_{instrument}.npy'
    runfile = planet_dir / 'runfile.dat'
    planet_name = planet_dir.name
    relative_posix = relative.as_posix().rstrip('/') + '/'
    content = (f'datapath\t{relative_posix}\n'
               f'mockdatapath\t{relative_posix}\n'
               f'mockdata\t{planet_name}.fits\n'
               f'transit\tTrue\n')
    if runfile.exists() and runfile.read_text() != content:
        raise FileExistsError(f'Refusing to overwrite existing runfile: {runfile}')
    _stage_resource(source_orders, destination_orders)
    if source_blaze is not None:
        _stage_resource(source_blaze, destination_blaze)
    # Without a supplied calibration, let tayph.mock handle its own blaze.
    # Do not fabricate unity arrays or replace existing user calibrations.
    if source_blaze is None and destination_blaze.exists():
        existing = np.load(destination_blaze, allow_pickle=False)
        if existing.shape != order_grid.shape:
            raise ValueError(f'Existing {instrument} blaze shape {existing.shape} does not match orders {order_grid.shape}')
    staged_blaze = destination_blaze if destination_blaze.exists() else None
    if not runfile.exists():
        runfile.write_text(content)
    return TayphPreparation(planet_dir, runfile, destination_orders, staged_blaze,
                            planet_dir / f'{planet_name}.fits')


def write_tayph_spectrum_fits(spectrum, destination: str | Path, *, overwrite: bool = False,
                             config_path: str | Path | None = None) -> Path:
    """Write the validated three-HDU legacy tayph model FITS format.

    Schema established from the user's working TRAPPIST-1_b.fits example:
    PrimaryHDU([0.0]), ImageHDU(wavelength_nm), ImageHDU(transmission).
    Input wavelength: vacuum microns (POSEIDON). FITS wavelength: air nm
    (matches the supplied ESPRESSO order definitions). Conversion uses the
    installed tayph.operations.vactoair exactly once. The transmission array
    is already 1 - POSEIDON output and must not be inverted a second time.
    """
    import numpy as np
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError('Writing tayph-compatible FITS requires astropy') from exc

    vacuum_nm = np.asarray(spectrum.wavelength_um, dtype=np.float64) * 1000.0
    if vacuum_nm.ndim != 1 or vacuum_nm.size < 2 or not np.isfinite(vacuum_nm).all():
        raise ValueError('Input vacuum wavelength grid must be finite and one-dimensional')
    if not np.all(np.diff(vacuum_nm) > 0) or vacuum_nm[0] <= 0:
        raise ValueError('Input vacuum wavelength grid must be positive and strictly ascending')
    if config_path is not None:
        config_path = Path(config_path).resolve()
        cfg = Config.from_file(config_path)
        configure_backend_path(cfg, key='tayph_path', package='tayph',
                               config_dir=config_path.parent)
    try:
        from tayph.operations import vactoair
    except ImportError as exc:
        raise ImportError('Air-wavelength tayph FITS export needs tayph.operations.vactoair; '
                          'set tayph_path in the planet config or install tayph') from exc
    wavelength_nm = np.asarray(vactoair(vacuum_nm), dtype=np.float64)
    transmission = np.asarray(spectrum.transmission, dtype=np.float64)
    if wavelength_nm.ndim != 1 or wavelength_nm.size < 2 or transmission.shape != wavelength_nm.shape:
        raise ValueError('Expected matching 1-D wavelength and transmission arrays')
    if not np.isfinite(wavelength_nm).all() or not np.all(np.diff(wavelength_nm) > 0):
        raise ValueError('FITS wavelength grid must be finite, positive and strictly ascending')
    if wavelength_nm[0] <= 0 or not np.isfinite(transmission).all():
        raise ValueError('FITS transmission and wavelengths must be finite with positive wavelengths')
    destination = Path(destination)
    if destination.exists() and not overwrite:
        raise FileExistsError(f'Refusing to overwrite existing tayph input FITS: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    hdus = fits.HDUList([
        fits.PrimaryHDU(np.array([0.0], dtype=np.float64)),
        fits.ImageHDU(wavelength_nm),
        fits.ImageHDU(transmission),
    ])
    hdus.writeto(destination, overwrite=overwrite, output_verify='exception')
    return destination


def _legacy_rvtrace(cfg, prepared, bjds, *, project_root):
    """Use the scientist's integrated radvel-backed RVTraceEstimator and plotter.

    Its plotting code requires an on-disk planet directory; create a diagnostic
    copy rather than editing the scientifically correct circular omega=nan.
    """
    import numpy as np
    from .rvtrace import RVTraceEstimator as estimator_cls

    diagnostics = project_root / 'diagnostics' / prepared.planet_dir.name
    stage = diagnostics / 'rvtrace_input'
    stage.mkdir(parents=True, exist_ok=True)
    stage_cfg = Config.from_file(prepared.planet_dir / 'config')
    # The integrated estimator interprets omega=nan for e=0 in memory.
    stage_cfg.write(stage / 'config')
    if bjds is not None:
        np.savetxt(stage / 'obs_times', np.asarray(bjds, dtype=float) - 2400000.5, fmt='%.12f')
    else:
        source_times = prepared.planet_dir / 'obs_times'
        if source_times.is_file():
            shutil.copy2(source_times, stage / 'obs_times')
        else:
            # Do not reuse potentially stale observation times from an earlier run.
            (stage / 'obs_times').unlink(missing_ok=True)
    rv = estimator_cls(str(stage), cfg.get('observatory', 'Paranal'))
    # The original constructor computes these already; do not recalculate without need.
    b = np.asarray(rv.bjds, dtype=float)
    f = np.asarray(rv.true_anomaly, dtype=float)
    if b.ndim != 1 or f.shape != b.shape or b.size < 2:
        raise ValueError('Integrated RVTraceEstimator returned incompatible BJDs/true anomalies')
    if not np.isfinite(b).all() or not np.isfinite(f).all():
        raise ValueError('Integrated RVTraceEstimator returned non-finite orbital values')
    if bjds is not None and not np.allclose(b, bjds, atol=1e-8, rtol=0):
        raise ValueError('Integrated RVTraceEstimator BJDs differ from supplied observation times')
    # Reuse the integrated RVTraceEstimator figure (do not reinvent its diagnostics).
    import matplotlib.pyplot as plt
    rv.plot_doppler_trace(RV_ext=cfg.get_float('rv_plot_extent', 200.0))
    fig = plt.gcf()
    if not any(ax.lines or ax.collections or ax.images for ax in fig.axes):
        raise RuntimeError('Integrated RVTraceEstimator.plot_doppler_trace() did not draw a plot. '
                           'Check whether the plotting code is commented out in your local RVTraceEstimator.py.')
    plot_path = diagnostics / 'rv_doppler_trace.png'
    fig.savefig(plot_path, dpi=170, bbox_inches='tight')
    plt.close(fig)
    return b, f, plot_path


def run_tayph_mock(config_path: str | Path, *, project_root: str | Path = '.',
                   bjds=None, dry_run: bool = False, mock_module=None):
    """Execute original tayph.mock with integrated radvel-backed RVTraceEstimator.

    The integrated estimator provides BJDs, true anomaly and Doppler-trace plots
    for both circular and eccentric orbits. This follows the user's proven
    workflow, including eccentric_orbit=True even for a circular trajectory.
    """
    import numpy as np

    config_path = Path(config_path).resolve()
    root = Path(project_root).resolve()
    cfg = Config.from_file(config_path)
    prepared = prepare_tayph_mock(config_path, project_root=root)
    if not prepared.expected_mock_fits.is_file():
        raise FileNotFoundError(f'Export model FITS first: {prepared.expected_mock_fits}')
    if bjds is None and 'mock_bjds_file' in cfg:
        path = Path(cfg['mock_bjds_file']).expanduser()
        if not path.is_absolute():
            path = config_path.parent / path
        bjds = np.loadtxt(path, dtype=float, ndmin=1)
    if bjds is not None:
        bjds = np.asarray(bjds, dtype=float)
        if bjds.ndim != 1 or bjds.size < 2 or not np.isfinite(bjds).all() or np.any(np.diff(bjds) <= 0):
            raise ValueError('Observation BJDs must be finite and strictly increasing (at least two)')
    kwargs = dict(create_obs_file=True, mode=canonical_instrument(cfg.get('instrument', 'ESPRESSO')),
                  real_data=False, spec='depth', rot=False,
                  eccentric_orbit=True, tellurics=cfg.get_bool('tellurics', False))
    if dry_run:
        # No imports/plots/expensive integration during dry runs.
        if bjds is not None:
            kwargs['bjds'] = bjds
        return prepared, kwargs
    if Path.cwd().resolve() != root:
        raise RuntimeError(f'Run tayph from project root {root}, not {Path.cwd()}')
    orbital_bjds, true_anomaly, plot_path = _legacy_rvtrace(cfg, prepared, bjds, project_root=root)
    kwargs.update(phase_range=true_anomaly, bjds=orbital_bjds)
    if mock_module is None:
        configure_backend_path(cfg, key='tayph_path', package='tayph', config_dir=config_path.parent)
        try:
            import tayph.mock as mock_module
        except ImportError as exc:
            raise ImportError('tayph.mock required; install tayph or set tayph_path in config') from exc
    mock_module.create_mock_obs(str(prepared.runfile), **kwargs)
    return prepared, kwargs
