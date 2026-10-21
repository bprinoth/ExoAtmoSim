"""Combine compatible Kp-Vsys CCF maps from multiple instruments and seeds."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np


def combine_maps(planet, *, instruments, seeds, species=('O2', 'O3', 'H2O'),
                 sysrem=5, data_dir='data', overwrite=False):
    """Sum raw CCF maps (not individually normalized detection significances).

    Every requested instrument/seed map must exist. Axes must agree exactly
    within floating-point tolerance. This is a sum of independent CCF
    contributions, not a sigma-combined statistical significance estimate.
    """
    from astropy.io import fits
    from .cc import expected_planet_kp
    from .plotting import publication_style
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
        maps, reference_axes = [], None
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
            elif not (np.allclose(kp, reference_axes[0], atol=1e-9, rtol=0) and
                      np.allclose(vsys, reference_axes[1], atol=1e-9, rtol=0)):
                raise ValueError(f'Incompatible Kp/Vsys axes: {paths[0]}')
            maps.append(image)
        kp, vsys = reference_axes
        destination = output / molecule
        paths_out = [destination / name for name in ('kp_vsys.fits', 'kp_axis.fits',
                                                       'vsys_axis.fits', 'kp_vsys.png',
                                                       'combination.json')]
        if not overwrite and any(path.exists() for path in paths_out):
            raise FileExistsError(f'Combined output exists: {destination}; use --overwrite')
        combined = np.sum(np.stack(maps), axis=0)
        kp_expected = expected_planet_kp(selected[0][0] / 'config')
        for run, _ in selected[1:]:
            other = expected_planet_kp(run / 'config')
            if not np.isclose(other, kp_expected, rtol=1e-5):
                raise ValueError(f'Inconsistent expected Kp across seeds: {run}')
        files[molecule] = (selected, combined, kp, vsys, kp_expected, destination)
    # Validate every requested species before writing any output.
    for molecule, (selected, combined, kp, vsys, kp_expected, dest) in files.items():
        dest.mkdir(parents=True, exist_ok=True)
        for name, array in [('kp_vsys.fits', combined), ('kp_axis.fits', kp),
                            ('vsys_axis.fits', vsys)]:
            fits.writeto(dest / name, array, overwrite=overwrite)
        with publication_style():
            fig, ax = plt.subplots(figsize=(9, 6))
            extent = [vsys[0], vsys[-1], kp[0], kp[-1]]
            image = ax.imshow(combined, origin='lower', aspect='auto', extent=extent)
            ax.axhline(kp_expected, color='white', linestyle='--',
                       linewidth=1, label=f'Expected Kp {kp_expected:.1f} km/s')
            ax.axvline(0, color='white', linestyle=':', linewidth=1)
            ax.set(xlabel='Delta Vsys [km/s]', ylabel='Kp [km/s]',
                   title=f'{planet} — {molecule} — {len(selected)} CCF maps')
            ax.legend(loc='upper right')
            fig.colorbar(image, ax=ax, label='Sum of CCF values (arbitrary units)')
            fig.savefig(dest / 'kp_vsys.png', dpi=180, bbox_inches='tight')
            plt.close(fig)
        (dest / 'combination.json').write_text(json.dumps({
            'planet': planet, 'molecule': molecule, 'sysrem': sysrem,
            'instruments': list(instruments), 'seeds': list(seeds),
            'n_maps': len(selected), 'expected_kp_kms': kp_expected,
            'aggregation': 'sum_of_raw_ccf_maps',
            'inputs': [str(paths[0]) for _, paths in selected],
        }, indent=2) + '\n')
    return output
