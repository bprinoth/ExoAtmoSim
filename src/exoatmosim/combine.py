"""Combine compatible Kp-Vsys CCF maps from multiple instruments and seeds."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np


def _run_snr_weight(run):
    """Median per-exposure VIOLA resolution-element S/N, squared.

    No S/N file means equal weighting; malformed existing data are errors.
    """
    from astropy.io import fits
    path = run / 'snr_res_0.fits'
    if not path.is_file():
        return 1.0, 'equal_weight_missing_snr'
    values = np.asarray(fits.getdata(path), float)
    if values.ndim != 2 or not values.size:
        raise ValueError(f'Invalid S/N shape: {path}')
    medians = []
    for exposure in values:
        valid = exposure[np.isfinite(exposure) & (exposure > 0)]
        if not valid.size:
            raise ValueError(f'No valid positive S/N for an exposure: {path}')
        medians.append(float(np.median(valid)))
    snr = float(np.median(medians))
    weight = snr ** 2
    if not np.isfinite(weight) or weight <= 0:
        raise ValueError(f'Invalid S/N weight: {path}')
    return weight, 'median_exposure_snr_squared'


def combine_maps(planet, *, instruments, seeds, species=('O2', 'O3', 'H2O'),
                 sysrem=5, data_dir='data', overwrite=False):
    """S/N-squared weighted mean of raw CCF maps.

    Every requested instrument/seed map must exist. Axes must agree exactly
    within floating-point tolerance. This is a sum of independent CCF
    contributions, not a calibrated detection significance estimate.
    """
    from astropy.io import fits
    from .cc import expected_planet_kp
    from POSEIDON.high_res import plot_CCF_Kp_Vsys
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if not planet.strip() or not instruments or not seeds or not species:
        raise ValueError('Specify planet, instruments, seeds, and species')
    if any(not str(x).strip() for x in instruments + species):
        raise ValueError('Blank instrument/species')
    if len(set(instruments)) != len(instruments) or len(set(seeds)) != len(seeds):
        raise ValueError('Duplicate instruments or seeds')
    if any(int(seed) < 0 for seed in seeds) or sysrem < 0:
        raise ValueError('Seeds and SYSREM must be nonnegative')
    root = Path(data_dir).resolve() / planet.replace(' ', '_')
    output = root / 'combined' / f'sysrem_{sysrem}'
    files = {}
    for molecule in species:
        selected = []
        for instrument in instruments:
            for seed in seeds:
                run = root / 'simulations' / instrument / f'seed_{int(seed):04d}'
                base = run / 'cc' / f'sysrem_{sysrem}' / molecule
                paths = [base / name for name in
                         ('kp_vsys.fits', 'kp_axis.fits', 'vsys_axis.fits')]
                for path in paths:
                    if not path.is_file():
                        raise FileNotFoundError(f'Missing CCF product: {path}')
                selected.append((run, paths))
        maps, weights, weight_info, reference_axes = [], [], {}, None
        for run, paths in selected:
            image, kp, vsys = (np.asarray(fits.getdata(p), float) for p in paths)
            if (kp.ndim != 1 or vsys.ndim != 1 or
                    image.shape != (len(kp), len(vsys)) or
                    not np.isfinite(image).all() or
                    not np.isfinite(kp).all() or not np.isfinite(vsys).all() or
                    np.any(np.diff(kp) <= 0) or np.any(np.diff(vsys) <= 0)):
                raise ValueError(f'Invalid CCF or velocity axes: {paths[0]}')
            if reference_axes is None:
                reference_axes = (kp, vsys)
            elif not (kp.shape == reference_axes[0].shape and
                      vsys.shape == reference_axes[1].shape and
                      np.allclose(kp, reference_axes[0], atol=1e-9, rtol=0) and
                      np.allclose(vsys, reference_axes[1], atol=1e-9, rtol=0)):
                raise ValueError(f'Incompatible Kp/Vsys axes: {paths[0]}')
            maps.append(image)
            weight, kind = _run_snr_weight(run)
            weights.append(weight)
            weight_info[f'{run.parent.name}/{run.name}'] = {
                'weight': float(weight), 'weight_kind': kind,
            }
        kp, vsys = reference_axes
        destination = output / molecule
        paths_out = [destination / name for name in ('kp_vsys.fits', 'kp_axis.fits',
                                                       'vsys_axis.fits', 'kp_vsys.png',
                                                       'combination.json')]
        if not overwrite and any(path.exists() for path in paths_out):
            raise FileExistsError(f'Combined output exists: {destination}; use --overwrite')
        combined = np.average(np.stack(maps), axis=0, weights=np.asarray(weights))
        kp_expected = expected_planet_kp(selected[0][0] / 'config')
        for run, _ in selected[1:]:
            other = expected_planet_kp(run / 'config')
            if not np.isclose(other, kp_expected, rtol=1e-5):
                raise ValueError(f'Inconsistent expected Kp across seeds: {run}')
        files[molecule] = (selected, combined, kp, vsys, kp_expected, destination, weight_info)
    # Validate every requested species before writing any output.
    for molecule, (selected, combined, kp, vsys, kp_expected, dest, weight_info) in files.items():
        dest.mkdir(parents=True, exist_ok=True)
        for name, array in [('kp_vsys.fits', combined), ('kp_axis.fits', kp),
                            ('vsys_axis.fits', vsys)]:
            fits.writeto(dest / name, array, overwrite=overwrite)
        # Feed the full combined raw map to POSEIDON's own significance
        # plotting routine, as in the single-run CC stage. Retain negative
        # Kp internally if POSEIDON needs it for noise normalization, but
        # never display that unphysical region.
        with plt.rc_context({'figure.constrained_layout.use': False,
                             'figure.autolayout': False}):
            fig, _ = plot_CCF_Kp_Vsys(
                kp, vsys, CCF_Kp_Vsys=combined, label=molecule,
                Kp_expected=kp_expected, Vsys_expected=0.,
                plot_slice=True, save_file_path=None)
            for ax in fig.axes:
                if ax.get_xlabel():
                    ax.set_xlabel('Vsys [km/s]')
                if ax.get_ylabel() and ('Kp' in ax.get_ylabel() or 'K_p' in ax.get_ylabel()):
                    ax.set_ylim(bottom=0)
                legend = ax.get_legend()
                if legend is not None:
                    legend.remove()
            # First panel is POSEIDON's Kp-Vsys image; ensure positive Kp only
            # even if the backend uses a different y-axis label.
            for ax in fig.axes:
                if ax.images:
                    ax.set_ylim(0, float(kp[-1]))
                    ax.set_xlabel('Vsys [km/s]')
                    ax.set_ylabel('Kp [km/s]')
            fig.savefig(dest / 'kp_vsys.png', dpi=180, bbox_inches='tight')
            plt.close(fig)
        (dest / 'combination.json').write_text(json.dumps({
            'planet': planet, 'molecule': molecule, 'sysrem': sysrem,
            'instruments': list(instruments), 'seeds': list(seeds),
            'n_maps': len(selected), 'expected_kp_kms': kp_expected,
            'aggregation': 'weighted_mean_of_raw_ccf_maps',
            'plot_normalization': 'POSEIDON.plot_CCF_Kp_Vsys',
            'plot_kp_min_kms': 0.0,
            'weight_kind': 'median_exposure_snr_squared_with_missing_file_equal_fallback',
            'weights': weight_info,
            'inputs': [str(paths[0]) for _, paths in selected],
        }, indent=2) + '\n')
    return output
