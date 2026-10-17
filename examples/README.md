# Runnable ExoAtmoSim examples

The example scripts contain **only public API calls and explanatory comments**. They are not separate simulators. Run from the repository root after `python -m pip install -e .` and configure the external scientific backends.

| Script | Instrument | Status |
| --- | --- | --- |
| `crires_k.py` | CRIRES+ K | Demonstrated full workflow |
| `espresso.py` | ESPRESSO | Constant-S/N development |
| `igrins.py` | IGRINS H/K | Constant-S/N development |
| `igrins2.py` | IGRINS2 H/K | Constant-S/N development |
| `spirou.py` | SPIRou | Constant-S/N development |
| `andes.py` | ANDES UBV/RIZ/YJH | Experimental; ETC wavelengths only |
| `andes_k.py` | ANDES K | Experimental provisional wavelength chunks |

For example: `python examples/crires_k.py`, or use the equivalent terminal command:

```bash
exoatmosim simulate --planet "TOI-1842 b" --instrument CRIRES_K --snr 100 --seed 42 --no-tellurics
```

`--instrument` selects the wavelength setup; it **does not** automatically select atmospheric species, chemistry or opacities. For optical iron examples, edit the planet-level `config` and choose appropriate Fe species/opacity settings first. Simulations are constant S/N and do not model calibrated throughput or detailed detector noise.

For independent stages, run `exoatmosim config --help`, `exoatmosim atmosphere --help`, `exoatmosim spectrum --help` and `exoatmosim observe --help`.
