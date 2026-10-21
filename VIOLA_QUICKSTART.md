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

## Stellar continuum without PHOENIX lines

To isolate the injected planetary transmission lines, set the instrument-specific
configuration to:

```text
stellar_source blackbody
```

The native VIOLA J/H/K simulator then uses a smooth Planck continuum at
`Teff`, calibrated to the corresponding J/H/K magnitude. The star's photon
counts and ETC noise remain present; this is **not** a no-star/no-noise mode.
The planet transmission, observer-frame RV shifts, instrumental LSF and any
enabled telluric absorption are retained. For a clean planetary-only line test,
also set `tellurics False` (and compare with the default `PHOENIX` setting).
Regenerate `exoatmosim observe ... --overwrite` for this setting to take effect;
the planetary spectrum/templates do not need to be regenerated. The config
snapshot stays in the seed directory for reproducibility.

## One-pass molecular templates and CCF recovery

For multiple species, avoid reloading POSEIDON's opacity tables for each
individual template:

```bash
exoatmosim spectrum data/TRAPPIST-1_b/config_viola_k --seed 0 --create-templates --overwrite
exoatmosim observe data/TRAPPIST-1_b/config_viola_k --seed 0 --overwrite
exoatmosim cc data/TRAPPIST-1_b/config_viola_k --seed 0 --species H2O --sysrem 3 --overwrite
```

The `--create-templates` option keeps the full chemical composition for the
pressure–radius/MMW calculation and generates every non-bulk gas's opacity
contribution in the same POSEIDON setup. All species must still have useful
line opacity to produce a detectable CCF; cold, depleted atmospheres can
produce nearly flat templates. The current implementation warns through a
raised validation error for peak-to-peak contribution amplitudes ≤ 1e-12.

Once individual CCF maps have been generated for the same planet, species and
velocity grids, combine instrument bands and/or seeds:

```bash
exoatmosim combine --planet "WASP-189 b" --instruments VIOLA_J VIOLA_H VIOLA_K \
  --seeds 1 2 3 4 --species O2 O3 H2O --sysrem 5 --overwrite
```

This uses the median per-exposure S/N squared for relative map weights when
`snr_res_0.fits` is available. The combined FITS map remains in raw CCF units;
the corresponding PNG uses POSEIDON's sigma-like normalization (not a
calibrated detection probability). For the full methodology see
[CC_QUICKSTART.md](CC_QUICKSTART.md).

## Per-exposure atmospheric and S/N diagnostics

The `observing_conditions.png` diagnostic shows PWV, airmass, and median
S/N per resolution element (derived from `snr_res_0.fits`). When
`tellurics False`, the PWV panel instead states "No tellurics considered".
