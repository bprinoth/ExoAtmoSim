# External telluric HDF5 integration — experimental

The builder v0.1.x writes one `telluric_<site>.h5` per observing site. Copy or download them **outside GitHub** into a local folder, then configure your private `config/site_paths.conf`:

```text
telluric_grid_dir /data/bibi/telluric_grids
```

The program reads this private path **at observation time**, including when
called from outside the source checkout. It is intentionally not copied into
planet-specific configs. Use `EXOATMOSIM_SITE_CONFIG=/absolute/path/to/site_paths.conf`
for installations where automatic discovery is not possible.

Planet config, e.g. WASP-189 b with ESPRESSO:

```text
instrument ESPRESSO
tellurics True
pwv 2.5
pwv_sigma 0.5
```

The grid must cover the entire instrument wavelength range, **all sampled PWV values**, and **all exposure airmasses**. Partial Paranal H-band pilot cannot be used for ESPRESSO. Out-of-range conditions and missing grid cells are errors, not invisible tellurics.

Run `exoatmosim spectrum data/WASP-189_b/config --seed 42 --overwrite` if no matching planetary FITS exists, then `exoatmosim observe data/WASP-189_b/config --seed 42 --overwrite`. Inspect `observing_conditions.ecsv` (including `pwv_mm`) and `diagnostics/telluric_transmission.png`.

**Scientific cautions:** Builder `allow_approx_pwv` uses an *uncalibrated* water-vapour mapping. This backend interpolates the stored spectra in optical depth. The existing constant-S/N observational noise model has not become a physical instrument ETC. The high-resolution telluric signal is evaluated at instrument sample centers; accurate convolution of stellar×telluric high-resolution products and transmission interpolation require further validation. There is no automatic download, no silent fallback to direct TelFit and no extrapolation. Existing observations with `tellurics False` are unaffected.


## Grid provenance and supported sites

The experimental `VIOLA_J`, `VIOLA_H`, and `VIOLA_K` wavelength configurations
and photon/noise calculations are informed by the independent
[VIOLA Exposure Time Calculator (viola-etc)](https://github.com/sirinrats/viola-etc)
by its upstream authors. ExoAtmoSim does **not** execute that ETC or redistribute
its sky-emission tables. The current VIOLA model adopts an 8 m reference
telescope (different from the ETC's default 2 m), uses the ExoAtmoSim PHOENIX
stellar pipeline, and omits sky-emission noise for now. Its sensitivity
estimates are therefore experimental, not validated VIOLA performance claims.

For `tellurics True`, ExoAtmoSim reads **externally generated, high-resolution
TelFit transmission grids** from HDF5 files; it does not call TelFit during
individual observation exposures. The observing schedule supplies each exposure's
airmass, while a seeded, time-correlated PWV series (mm) supplies the water-vapour
coordinate. Transmission is interpolated **in optical depth** across neighboring
airmass and PWV grid nodes, then used in the exposure's spectral synthesis.
The transmission grids are independent of VIOLA's ETC and contain
**absorption only**, not sky emission. Out-of-range conditions or missing grid
coverage are errors rather than an unannounced clear-sky fallback.

The telluric-grid builder supports these observing-site identifiers:

| Site | Grid identifier |
| --- | --- |
| ESO Paranal Observatory (Chile) | `paranal` |
| Cerro Armazones (Chile) | `armazones` |
| Mauna Kea (Hawaiʻi) | `mauna_kea` |
| Gemini South / Cerro Pachón (Chile) | `gemini_south` |
| Roque de los Muchachos / La Palma (Canary Islands) | `la_palma` |

These are **supported site definitions**, not a claim that all five production
grids are already complete or validated. In particular, the Gemini South
LBLRTM run has encountered failures. Install only verified grid files, kept
outside GitHub, and configure their directory in the private
`config/site_paths.conf` as `telluric_grid_dir /absolute/path/to/grids`.
See [TELLURIC_GRID_QUICKSTART.md](TELLURIC_GRID_QUICKSTART.md) for setup.
The builder's provisional conversion of PWV to TelFit relative humidity
requires scientific validation.
