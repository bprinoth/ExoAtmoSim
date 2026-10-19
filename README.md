# ExoAtmoSim

ExoAtmoSim is a modular Python framework for synthesizing high-resolution exoplanet transmission observations. Its current demonstrator runs **NASA Exoplanet Archive → FastChem → POSEIDON → PHOENIX → stellar/planetary velocities → instrument-sampled, constant-S/N time-series spectra**.

> **Scientific development version.** The current public baseline uses a user-specified constant S/N and should not be treated as a physically calibrated detector, throughput, or noise prediction.

## Install and discover

```bash
python -m pip install -e .
exoatmosim --help
exoatmosim simulate --help
exoatmosim atmosphere --help
python -m unittest discover -s tests -q
```

FastChem, POSEIDON, PHOENIX archival retrieval, and orbital dependencies must be configured on the user's system. Copy the checked-in `config/site_paths.example.conf` to **`config/site_paths.conf`**, set your local FastChem/POSEIDON/data paths there, and do not commit that machine-specific file (it is Git-ignored). The CLI detects this file automatically when run from the project root; alternatively set `EXOATMOSIM_SITE_CONFIG` or pass `--site-config`. No FastChem or POSEIDON installation is bundled, and `pip install -e .` does not install them.

## Modular commands and Python API

The public CLI supports five stages: `config` retrieves/reuses a planet configuration; `atmosphere` computes the T–P and chemistry diagnostics; `spectrum` computes a POSEIDON transmission spectrum and writes the seed-specific FITS model; `observe` uses that existing FITS model to synthesize observations without re-running radiative transfer; and `simulate` runs the end-to-end workflow. Use `exoatmosim <stage> --help` for arguments.

```bash
exoatmosim config --planet "TOI-1842 b" --instrument CRIRES_K
exoatmosim atmosphere data/TOI-1842_b/config --seed 42
exoatmosim spectrum data/TOI-1842_b/config --backend poseidon --seed 42
exoatmosim observe data/TOI-1842_b/config --seed 42 --snr 100
```

Equivalent Python functions are exposed from `exoatmosim`: `make_config`, `calculate_atmosphere`, `calculate_spectrum`, `simulate_observation`, and `simulate`. `observe` requires an existing matching planetary model; it will not calculate a missing one. For now, the supported radiative-transfer backend is POSEIDON; petitRADTRANS is planned.

Examples in `examples/` are short, commented API usage demonstrations, not competing implementations. See `examples/README.md`.

## Example complete simulation

```bash
exoatmosim simulate --planet 'TOI-1842 b' --instrument CRIRES_K --snr 100 --seed 42 --no-tellurics
```

The same command interface accepts ESPRESSO, IGRINS_HK, IGRINS2_HK, SPIROU, ANDES, and ANDES_K wavelength setups. The non-CRIRES modes are development instrument implementations requiring scientific validation.

| Instrument | Supplied wavelength setup | Status |
| --- | --- | --- |
| `CRIRES_K` | Seven K-band orders | End-to-end demonstrated on mantis |
| `ESPRESSO` | 170 sorted echelle orders | Constant-S/N development |
| `IGRINS_HK`, `IGRINS2_HK` | Supplied H/K sampling | Constant-S/N development |
| `SPIROU` | Supplied wavelength orders | Constant-S/N development |
| `ANDES` | UBV, RIZ, YJH ETC wavelengths, preserving overlaps and input order | **Experimental**; wavelength inputs only, never ETC noise/flux |
| `ANDES_K` | Provisional K grid split into computational chunks | **Experimental**, chunks are not calibrated echelle orders |

The simulator saves a reusable planet `config` under `data/<planet>/`, while each seed also receives its **own frozen `config`**, planetary FITS model, observing conditions, spectral-order FITS outputs, manifest and diagnostics under `data/<planet>/simulations/<instrument>/seed_<seed>/`.

## Examples and limitations

TOI-1842 b with CRIRES+ K demonstrates the complete FastChem/POSEIDON-to-observations chain, including eccentric orbital velocities. WASP-189 b and ESPRESSO are development examples for optical/Fe simulations; supported wavelength setups have been exercised with constant S/N.

**Under active development:** VIOLA J/H/K is an experimental native PHOENIX + HDF5 telluric grid + ETC-inspired detector-noise backend. Sky emission is not yet modeled, so VIOLA S/N is optimistic. Telluric and noise validation, emission spectroscopy, and petitRADTRANS remain under development.

The order simulator is independent of external legacy mock-observation software. See code docstrings for numerical and astrophysical assumptions.

### Observation-window duration

`observation_window_hours` is the **total** requested observing duration in hours,
centered on `Tobs` when specified or on the selected maximum-altitude epoch.
New NASA-generated planet configs default it to `T14 + 2` hours (one hour of
baseline on each side of the transit), when a positive catalog transit duration
`T14` is available. You can override the duration by editing
`observation_window_hours` in the planet's config. If the catalog has no usable
`T14`, choose the observing duration explicitly. The old ambiguous `t` field is
not used and should be removed from older configs. `T14` remains the physical
transit duration; `observation_window_hours` controls the simulated schedule.

### Exposure duration configuration

Use `exposure_time` (seconds) for the integration time of **each exposure**,
`readout` (seconds) for dead time between exposures, and
`observation_window_hours` (hours) for the complete observation window.
The old `texp` and `exp_time` configuration keys are no longer accepted; update
existing planet configs before simulating. The scheduler will explicitly report
a missing `exposure_time` instead of silently choosing a different parameter.

### Kaltenegger / Rugheimer Earth-analogue chemistry (experimental)

The third-party profile data are **not distributed with ExoAtmoSim**. Download
and extract them from [Zenodo record 14829587](https://zenodo.org/records/14829587).
Follow the archive's license and citation instructions. The expected layout is:

```
Kaltenegger_profiles/
    PT_profiles/      # AltPT_<stellar type>_<epoch>Ga.txt
    chem_profiles/    # Chem_<stellar type>_<epoch>Ga.txt
```

Set the absolute path to the **parent folder** in your machine-specific
`config/site_paths.conf` (copy `config/site_paths.example.conf` if needed):

```
kaltenegger_profiles_dir /path/to/Kaltenegger_profiles
```

Then use `chemistry kaltenegger` in a copy of the planet configuration. Neither
the original NASA `stellar_type` nor `Gyr_age` is modified. The nearest
available stellar template (F0V, F7V, G8V, K2V, K7V, M1A, M3A, M8A, Sun)
and complete profile age (0.0, 0.8, 2.0, 3.9 Ga lookback) are chosen only
for this chemistry backend. The age conversion uses an approximate Earth-
analogue clock `max(0, 4.6 - Gyr_age)`; this is not a measured geological age
of another planet.

PT files contain `Z,T,P` records (km, K, mbar). Chemistry files contain
`Molecule`, `<species>_MR` and `<species>_ALT` records, with multiline `^`
continuation. Missing N2 is calculated as the residual of the **complete**
chemistry profile, before unsupported POSEIDON opacity species are omitted;
the retained mixing ratios and N2 are never renormalized. Unsupported
absorption is reported as a modeling limitation. The FastChem workflow
remains unchanged. Related study: https://arxiv.org/abs/1912.11149 .

Run `exoatmosim atmosphere`, followed by `exoatmosim spectrum` and
`exoatmosim observe` with the same planet config and seed, as above.

### Contributions and integration requests

**Integration suggestions are welcome!** If you would like ExoAtmoSim to
support another atmospheric chemistry source, radiative-transfer backend,
instrument, or observing workflow, please contact the project maintainer via
the [author's website](https://bprinoth.github.io/) or open a GitHub issue
at https://github.com/bprinoth/ExoAtmoSim/issues . Contributions and
scientific use cases are very welcome.


## Experimental VIOLA J/H/K and HDF5 tellurics

VIOLA_J, VIOLA_H and VIOLA_K use the packaged wavelength arrays, PHOENIX
stellar models shifted per exposure, the planet's separate Doppler shift, an
8 m reference telescope, R=300,000, and VIOLA ETC-inspired electron noise.
The **external VIOLA ETC is not executed**. VIOLA sky emission is **not modeled**
in this release: its S/N is an optimistic upper limit.

With `tellurics True`, the same externally stored TelFit HDF5 transmission
grids are used for VIOLA and the other supported instruments. The atmospheric
transmission is selected at each exposure's airmass and sampled PWV, with
bilinear interpolation in optical depth. It is not derived from VIOLA's
sky-emission library. The grid-builder's approximate PWV calibration is
experimental and needs validation.

In a planet config set only the relevant atmospheric parameters:

```text
instrument VIOLA_H
resolution 300000
tellurics True
pwv 2.5
pwv_sigma 0.5
```

Keep the H-band magnitude `Hmag` (or Jmag/Kmag for those modes) in the
planet config. Put your **private, machine-local** HDF5 directory in
`config/site_paths.conf`:

```text
telluric_grid_dir /absolute/path/to/telluric_grids
```

The local path is resolved **at observation time**, including when running
from outside the checkout. It is deliberately not copied into planet-specific
configs or bundled into this repository. If needed, set
`EXOATMOSIM_SITE_CONFIG` to an absolute path to your private file.
See `TELLURIC_GRID_QUICKSTART.md` and `VIOLA_QUICKSTART.md`.

```bash
python -m pip install -e '.[tellurics]'
exoatmosim spectrum data/YOUR_PLANET/config --seed 42 --overwrite
exoatmosim observe data/YOUR_PLANET/config --seed 42 --overwrite
```

VIOLA uses ETC-inspired noise rather than a requested `--snr`; other
instrument modes use the existing constant-S/N approximation.

### Outputs and important limitations

`order_0.fits`, `wave_0.fits`, `stellar_flux_0.fits`, `model_flux_0.fits`,
`noise_0.fits`, `snr_res_0.fits`, `viola_native_manifest.json`.
Counts and uncertainties are **per resolution element evaluated at each
supplied wavelength center**, *not* per independent detector pixel. Oversampled
adjacent samples have correlated physical noise, but the first implementation
draws independent Gaussian errors; **do not interpret its cross-correlation
significance as validated**. The native wavelength grids are computational
segments, not measured echelle order boundaries. The ETC 2x2 detector-pixel
noise convention is retained. Magnitudes anchor PHOENIX at a single reference
wavelength (not full synthetic-filter photometry). Atmospheric transmission
uses TelFit on a finer internal grid, and combined star × planet × tellurics
are convolved to R=300,000 before sampling. No extra ETC stellar Doppler or
telluric multiplication is performed. Validation of sky emission, TelFit
wavelength convention, photometry and detector covariance remains necessary.

For test parameters: `viola_seeing_arcsec` (default 1),
`viola_slit_width_arcsec` (1.6), `viola_optics_throughput` (product of
ETC default factors), `viola_nx`/`viola_ny` (2 each),
`viola_read_noise_e_pixel` (20), `viola_dark_rate_e_s_pixel` (0.01),
`viola_nod_subtraction` (True), `viola_ambient_temperature_k` (275),
`viola_thermal_emissivity` (0.1), `viola_internal_sampling_r` (900000).


## Experimental external telluric HDF5 grids

See [TELLURIC_GRID_QUICKSTART.md](TELLURIC_GRID_QUICKSTART.md). Current planet configs use `tellurics`, `pwv` and `pwv_sigma`. No direct TelFit invocation in the grid-mode observation path; HDF5 grids live outside the repository.


### VIOLA ETC provenance and telluric sites

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


### External telluric grid distribution

The high-resolution HDF5 telluric models are **not** part of the Git
repository or the Python package. Obtain the separately distributed
`telluric_<site>.h5` data from the project maintainer (a permanent download
URL has not yet been published), and configure their containing directory
only in your private `config/site_paths.conf` using `telluric_grid_dir`.
A single grid library is shared by all supported instruments. The initial
PWV-to-RH mapping used in grid production is approximate and must be
validated before interpreting atmospheric retrieval or sensitivity results.

The native VIOLA backend uses its provided J/H/K sampling grids, PHOENIX
stellar flux anchored to Vega magnitudes, an 8 m reference telescope and
ETC-inspired detector variance. Its atmosphere is drawn from the same HDF5
library, not from the external VIOLA ETC. VIOLA sky emission is omitted, so
its present noise predictions are optimistic. The native model and numerical
limitations are documented in `VIOLA_QUICKSTART.md`.
