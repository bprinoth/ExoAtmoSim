# ExoAtmoSim

**From planetary parameters to synthetic high-resolution transmission spectroscopy.**

ExoAtmoSim is an actively developed Python toolkit for reproducible exoplanet atmosphere and simulated high-resolution observation workflows. The initial public foundation demonstrates **NASA Exoplanet Archive → FastChem → POSEIDON → PHOENIX → orbital velocities and constant-S/N spectral orders**.

> **Development foundation, not a calibrated observatory simulator.** Instrument wavelength sampling and fixed user-specified S/N do not imply realistic detector throughput or noise.

## Instrument modes

| Instrument flag | Wavelength setup |
| --- | --- |
| `CRIRES_K` | Seven K-band spectral orders |
| `ESPRESSO` | 170 spectral orders, corrected wavelength order sorting |
| `IGRINS_HK` | H and K spectral orders |
| `IGRINS2_HK` | H and K spectral orders |
| `SPIROU` | Supplied wavelength orders |
| `ANDES` | UBV, RIZ and YJH ETC wavelengths only, including all overlaps; experimental |
| `ANDES_K` | Provisional K wavelength grid partitioned into computational segments, **not physical echelle orders**; experimental |

**All currently documented examples use constant S/N.** The ANDES ETC products provide wavelength samples only; their S/N or flux fields are **not** used.

## Quick start

```bash
python -m pip install -e .
exoatmosim simulate --help
exoatmosim simulate --planet "TOI-1842 b" --instrument CRIRES_K --snr 100 --seed 42 --no-tellurics
```

Full simulations require separately installed and configured scientific backends, including FastChem, POSEIDON, stellar-spectrum retrieval and orbital modelling dependencies. See `examples/site_paths.conf` for backend configuration when the implementation files are present.

The planet-level `data/<planet>/config` remains a reusable template. Each simulation stores its effective `config`, planet spectrum FITS, observation products, run manifest and diagnostics beneath `data/<planet>/simulations/<instrument>/seed_<seed>/`.

## In development

Telluric absorption and atmospheric sampling, VIOLA integration, non-constant/instrument-specific S/N, emission spectroscopy, and petitRADTRANS support are being developed. These are **not part of the minimal validated no-tellurics release**.

The simulated spectral-order workflow is independent of an external mock-observation runtime.
