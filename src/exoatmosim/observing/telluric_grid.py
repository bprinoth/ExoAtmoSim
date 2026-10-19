"""Read external TelFit HDF5 grid-builder v0.1.x products, in AIR nm.

Interpolation is bilinear in optical depth, not transmission. The PWV axis
is only as physically calibrated as the generating builder's mapping.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np

SITES = {'gemini-south': 'gemini_south', 'gemini_south': 'gemini_south',
         'la-palma': 'la_palma', 'lapalma': 'la_palma', 'roque': 'la_palma',
         'orm': 'la_palma', 'maunakea': 'mauna_kea', 'mauna_kea': 'mauna_kea',
         'paranal': 'paranal', 'armazones': 'armazones'}


def resolve_grid_file(site, *, config_path=None):
    """Find an external telluric library from private site settings at run time.

    No machine-specific paths are copied into planet configs or committed to Git.
    Search from the planet config's checkout as well as cwd/installed checkout;
    an explicit EXOATMOSIM_SITE_CONFIG overrides discovery.
    """
    from ..backend_paths import discover_site_config
    from ..config import Config

    path = discover_site_config(config_path=config_path)
    if path is None:
        raise FileNotFoundError(
            'Telluric HDF5 library not configured. Set telluric_grid_dir in '
            'config/site_paths.conf or EXOATMOSIM_SITE_CONFIG.')
    cfg = Config.from_file(path)
    source = cfg.get('telluric_grid_dir')
    if not source:
        raise FileNotFoundError(f'Set telluric_grid_dir in {path} (external HDF5 library)')
    base = Path(source).expanduser()
    if not base.is_absolute():
        base = path.parent / base
    base = base.resolve()
    key = str(site).strip().lower().replace(' ', '_')
    if key not in SITES:
        raise ValueError(f'No telluric library mapping for observatory {site!r}; supported: {sorted(SITES)}')
    result = base / f'telluric_{SITES[key]}.h5'
    if not result.is_file():
        raise FileNotFoundError(f'Atmospheric grid absent: {result}')
    return result


def _bounds(axis, value, label):
    axis = np.asarray(axis, float)
    if axis.ndim != 1 or not np.all(np.isfinite(axis)) or not np.all(np.diff(axis)>0):
        raise ValueError(f'Invalid telluric {label} axis')
    if not np.isfinite(value) or not axis[0] <= value <= axis[-1]:
        raise ValueError(f'{label}={value:g} outside supported grid [{axis[0]:g}, {axis[-1]:g}]')
    hi = int(np.searchsorted(axis, value, side='right'))
    lo = min(max(hi-1,0),len(axis)-1)
    hi = min(hi,len(axis)-1)
    weight = 0. if lo==hi else (value-axis[lo])/(axis[hi]-axis[lo])
    return lo,hi,float(weight)


def _select_chunk(g, wave, ix, ip):
    arr = np.asarray(g[f'transmission_a{ix:02d}_p{ip:02d}'][:], float)
    if arr.shape != wave.shape or not np.isfinite(arr).all() or np.any((arr<0)|(arr>1)):
        raise ValueError('Invalid or incomplete telluric grid cell')
    return -np.log(np.maximum(arr, 1e-10))


def grid_transmission(schedule, wavelength_air_nm, *, pwv_mm, grid_file,
                      allow_approximate_pwv=True):
    """Return (AIR-nm wavelengths, Nexp x Nwave transmission).

    Only loads chunks intersecting caller's wavelengths; never requires the
    entire multi-GB grid in memory. No extrapolation or missing cells allowed.
    """
    try:
        import h5py
    except ImportError as exc:
        raise ImportError('HDF5 grid tellurics need h5py; pip install h5py') from exc
    grid_file = Path(grid_file)
    wl = np.asarray(wavelength_air_nm,float)
    water = np.asarray(pwv_mm,float)
    airmass = np.asarray(schedule.airmass,float)
    if (wl.ndim != 1 or len(wl)<2 or not np.isfinite(wl).all() or
            np.any(np.diff(wl)<=0) or water.shape!=airmass.shape or not len(water)):
        raise ValueError('Invalid telluric target wavelengths or exposure atmospheric arrays')
    output = np.empty((len(water),len(wl)),dtype=np.float32)
    covered = np.zeros(len(wl),dtype=bool)
    with h5py.File(grid_file,'r') as f:
        cfg = json.loads(f.attrs['config_json'])
        if cfg.get('wavelength_medium')!='air':
            raise ValueError('Telluric HDF5 must use air wavelengths')
        if cfg.get('allow_approx_pwv') and not allow_approximate_pwv:
            raise ValueError('Grid uses uncalibrated approximate PWV-to-RH mapping')
        if cfg.get('allow_approx_pwv'):
            import warnings
            warnings.warn('Telluric grid PWV is approximately mapped to RH; not physically PWV-calibrated',RuntimeWarning)
        ax=np.asarray(f['airmass'][:],float)
        px=np.asarray(f['pwv_mm'][:],float)
        corners=[(*_bounds(ax,a,'airmass'),*_bounds(px,p,'PWV mm')) for a,p in zip(airmass,water)]
        keys=sorted((k for k in f if k.startswith('chunk_')))
        if not keys: raise ValueError('No HDF5 telluric chunks found')
        chunk_edges=[]
        for key in keys:
            item=f[key]['wavelength_air_nm']
            chunk_edges.append((key,float(item[0]),float(item[-1])))
        # A globally sampled log grid is stored in disjoint HDF5 chunks.
        # The interval between each pair of chunks requires interpolation
        # across its two boundary points, not artificial transmission=1.
        for key in keys:
            group=f[key]
            x=group['wavelength_air_nm'][:]
            if x.ndim!=1 or len(x)<2 or not np.all(np.diff(x)>0):
                raise ValueError(f'Invalid wavelengths in {key}')
            sel=np.flatnonzero((wl>=x[0])&(wl<=x[-1]))
            if len(sel)==0: continue
            # Global builder uses disjoint chunks. Wavelengths in gaps must fail.
            covered[sel]=True
            request=wl[sel]
            cell_cache={}
            for row,(a0,a1,wa,p0,p1,wp) in enumerate(corners):
                for ai,pi in ((a0,p0),(a0,p1),(a1,p0),(a1,p1)):
                    if (ai,pi) not in cell_cache:
                        cell_cache[ai,pi]=np.interp(request,x,_select_chunk(group,x,ai,pi))
                t00,t01=cell_cache[a0,p0],cell_cache[a0,p1]
                t10,t11=cell_cache[a1,p0],cell_cache[a1,p1]
                tau=(1-wa)*((1-wp)*t00+wp*t01)+wa*((1-wp)*t10+wp*t11)
                output[row,sel]=np.exp(-tau)
        for (kl,_,last),(kr,first,_) in zip(chunk_edges,chunk_edges[1:]):
            gap_idx=np.flatnonzero((wl>last)&(wl<first))
            if len(gap_idx)==0: continue
            left,right=f[kl],f[kr]
            xl=np.asarray(left['wavelength_air_nm'][-1:],float)
            xr=np.asarray(right['wavelength_air_nm'][:1],float)
            if not xr[0]>xl[0]:
                raise ValueError('Non-monotone boundaries in HDF5 chunks')
            fraction=(wl[gap_idx]-xl[0])/(xr[0]-xl[0])
            corner_cache={}
            for row,(a0,a1,wa,p0,p1,wp) in enumerate(corners):
                for ai,pi in ((a0,p0),(a0,p1),(a1,p0),(a1,p1)):
                    if (ai,pi) not in corner_cache:
                        dataset=f'transmission_a{ai:02d}_p{pi:02d}'
                        if dataset not in left or dataset not in right:
                            raise ValueError(f'Missing HDF5 boundary cell {dataset}')
                        v1=float(left[dataset][-1]);v2=float(right[dataset][0])
                        if not np.isfinite([v1,v2]).all() or min(v1,v2)<0 or max(v1,v2)>1:
                            raise ValueError('Invalid HDF5 chunk-boundary transmission')
                        corner_cache[ai,pi]=((1-fraction)*(-np.log(max(v1,1e-10))) +
                                              fraction*(-np.log(max(v2,1e-10))))
                tau=(1-wa)*((1-wp)*corner_cache[a0,p0]+wp*corner_cache[a0,p1])+wa*((1-wp)*corner_cache[a1,p0]+wp*corner_cache[a1,p1])
                output[row,gap_idx]=np.exp(-tau)
            covered[gap_idx]=True
        if not np.all(covered):
            i=np.flatnonzero(~covered)[0]
            raise ValueError(f'HDF5 telluric grid has no wavelength coverage at {wl[i]:.6f} air nm')
    return wl.copy(),output
