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

**Under active development, not validated in this no-tellurics release:** telluric modelling and sampling, VIOLA, non-constant/instrument-specific S/N, emission spectroscopy, and a petitRADTRANS spectral backend.

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

### Optional Kaltenegger/Rugheimer Earth-like chemistry

Download the third-party atmosphere profiles from
[Zenodo record 14829587](https://zenodo.org/records/14829587).
**They are not bundled or redistributed with ExoAtmoSim.** Extract
the files preserving this layout:

```text
Kaltenegger_profiles/
  PT_profiles/AltPT_F0V_0.0Ga.txt
  chem_profiles/Chem_F0V_0.0Ga.txt
  ... (all matching files)
```

Copy `config/site_paths.example.conf` to your private
`config/site_paths.conf` and set:

```text
kaltenegger_profiles_dir /absolute/path/to/Kaltenegger_profiles
```

Set `chemistry kaltenegger` in a copy of your planet config; other
settings and FastChem behavior remain unchanged. The original NASA
`stellar_type` and `Gyr_age` values are preserved. Only when
Kaltenegger chemistry is selected, ExoAtmoSim maps to the nearest
available F0V, F7V, G8V, K2V, K7V, M1A, M3A, M8A or Sun template
and 0.0, 0.8, 2.0 or 3.9 Ga geological lookback profile.

The lookback mapping, `max(0, 4.6 - Gyr_age)`, is an **approximate
Earth-analogue clock**, not the geological age of an exoplanet.
The PT format has `Z,T,P` columns (km, K, mbar); abundances have
`Molecule`, `<species>_MR` and `<species>_ALT` continuation records.
N2 is computed from all abundances *before* dropping species without
POSEIDON opacity; the remaining VMRs are not renormalized. Omitted
opacity species are reported as a limitation. Related work:
https://arxiv.org/abs/1912.11149 .

Run `exoatmosim atmosphere`, `exoatmosim spectrum`, then
`exoatmosim observe` using the same configuration and seed.

### Contributions and integration requests

**Integration wishes are welcome!** If you would like support for another
atmospheric chemistry model, radiative-transfer backend, instrument, or
observing workflow, please [get in touch](https://bprinoth.github.io/)
or [open a GitHub issue](https://github.com/bprinoth/ExoAtmoSim/issues).
Scientific collaborations and contributions are very welcome.
