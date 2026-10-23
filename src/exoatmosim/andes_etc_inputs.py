"""Input-folder preparation for manually downloaded ANDES ETC curves.

The nine CSVs are supplied by the user; this module never fabricates S/N data.
"""
from pathlib import Path

AIRMASSES = ('1.0', '1.5', '3.0')
ARMS = ('UBV', 'RIZ', 'YJH')


def scaffold_andes_etc(planet_dir: str | Path) -> Path:
    """Prepare non-destructive per-airmass directories and print ETC instructions."""
    planet_dir = Path(planet_dir)
    root = planet_dir / 'andes_etc'
    for airmass in AIRMASSES:
        (root / f'airmass_{airmass}').mkdir(parents=True, exist_ok=True)

    instructions = """# ANDES ETC inputs for ExoAtmoSim

Generate **nine** official ANDES ETC CSV results manually:
one for each spectrograph arm (UBV, RIZ, YJH) at each of
airmasses 1.0, 1.5, and 3.0.

Save the exported CSVs at these paths:

    airmass_1.0/UBV.csv
    airmass_1.0/RIZ.csv
    airmass_1.0/YJH.csv
    airmass_1.5/UBV.csv
    airmass_1.5/RIZ.csv
    airmass_1.5/YJH.csv
    airmass_3.0/UBV.csv
    airmass_3.0/RIZ.csv
    airmass_3.0/YJH.csv

Use the existing `exposure_time` in the planet's config as the
**ETC integration time in seconds**, too. No separate reference exposure
keyword is needed. The `readout` config key is an observation overhead,
NOT ETC integration time.

Keep the stellar template, magnitude and band/system, seeing, transparency,
binning and NDIT identical across the nine ETC calculations. Set moon phase
to 0. Only airmass should change. Choose inputs appropriate to the actual
target; the earlier example G2V/V=12 is not a TRAPPIST-1 calibration.

The S/N output should be the ETC's **per resolution element** quantity.
Do not paste, rename, or convert columns within CSV files; rename only
the exported file itself to UBV.csv, RIZ.csv, or YJH.csv.

Generating folders does NOT mean ANDES ETC noise is implemented or calibrated.
The observe backend will need an explicit nine-file reader/interpolator
before these files can determine simulated noise.
"""
    guide = root / 'README.md'
    if not guide.exists():
        guide.write_text(instructions)
    print(f'[ANDES ETC] Input folders ready: {root}')
    print('[ANDES ETC] Download 9 CSVs from https://andes-etc.brera.inaf.it/')
    print('[ANDES ETC] Required: UBV.csv, RIZ.csv, YJH.csv in each of '
          'airmass_1.0, airmass_1.5 and airmass_3.0.')
    print('[ANDES ETC] Use config exposure_time as the ETC total integration '
          'time (seconds); moon phase = 0.')
    print(f'[ANDES ETC] Full instructions: {guide}')
    return root
