# Experimental VIOLA J/H/K quickstart

The native VIOLA backend uses PHOENIX stellar spectra (shifted once per
exposure), separate planet RVs, packaged VIOLA wavelength grids, R=300,000,
an **8 m reference telescope**, Vega J/H/K magnitudes from the planet config,
and ETC-inspired photon/detector noise. No external VIOLA ETC invocation.

VIOLA sky emission is **not modeled** in this version. The resulting noise
and S/N are optimistic; do not use them as calibrated sensitivity forecasts.

Use the same external HDF5 telluric transmission library for all instruments.
Place the grid files outside the repository and set the directory only in
your private `config/site_paths.conf`:

```text
telluric_grid_dir /absolute/path/to/telluric_grids
```

The site path is looked up at observation time and is not copied into planet
configs. If the working directory is separate from the source checkout, use
`EXOATMOSIM_SITE_CONFIG=/absolute/path/to/site_paths.conf` if necessary.

Planet config, using the automatically retrieved `Hmag` for H band:

```text
instrument VIOLA_H
resolution 300000
tellurics True
pwv 2.5
pwv_sigma 0.5
```

```bash
python -m pip install -e '.[tellurics]'
python -m unittest discover -s tests -q
exoatmosim spectrum data/TRAPPIST-1_b/config --seed 1 --overwrite
exoatmosim observe data/TRAPPIST-1_b/config --seed 1 --overwrite
```

For J or K, change `instrument` to `VIOLA_J` or `VIOLA_K` and provide the
corresponding magnitude. The grids must cover the entire band and the
sampled exposure airmass/PWV range. The grid builder's approximate PWV-to-RH
mapping needs scientific validation. No silent fallback to clear atmosphere.

VIOLA ETC origin: https://github.com/sirinrats/viola-etc . The ETC informs the
noise equations and configuration assumptions, not the telluric HDF5 data.
The 8 m telescope is an ExoAtmoSim reference choice, not the ETC default.

## Diagnostic plots

Every successful native VIOLA observation writes two plots under
`simulations/VIOLA_<band>/seed_<seed>/diagnostics/`:

- `planet_transmission.png` — the saved POSEIDON transmission model used for injection.
- `viola_spectral_comparison.png` — noise-free star+tellurics versus
  planet-injected electron expectations for the exposure with the strongest
  transit weight, plus their fractional difference.

They are regenerated when `observe --overwrite` runs. Plotting does not
recompute PHOENIX, POSEIDON, or the HDF5 telluric grid.
