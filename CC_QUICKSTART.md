# High-resolution cross-correlation (experimental)

The `cc` stage operates on a **previously generated observation** and the
**already-generated species templates** for that same seed. It does not repeat
POSEIDON radiative transfer, recalculate templates, or alter observation FITS.

Example using WASP-189 b on Mantis:

```bash
conda activate exoatmosim
git fetch origin feature/poseidon-cc
git switch feature/poseidon-cc
pip install -e . --no-deps
pytest tests/test_cc.py -v

exoatmosim cc data/WASP-189_b/config --seed 0 --sysrem 5 \
  --species O2 O3 H2O --prepare-only

exoatmosim cc data/WASP-189_b/config --seed 0 --sysrem 5 \
  --species O2 O3 H2O --cc-only
```

Use `--overwrite` to regenerate existing prepared files or CCF products.
Without `--prepare-only` or `--cc-only`, both stages execute. Avoid
`--overwrite` when you want to preserve an existing CCF.

Expected input: `data/WASP-189_b/simulations/VIOLA_H/seed_0000/`
(or the instrument from the configuration) containing `config`,
`obs_times`, `diagnostics/observation_schedule.csv`,
`order_0.fits`, `wave_0.fits`, and `templates/O2.fits` etc.

Outputs: `<run>/cc/sysrem_5/prepared/order_0/` and
`<run>/cc/sysrem_5/O2/`, `O3/`, `H2O/`, containing CCF FITS arrays and figures.

**Dependencies:** editable POSEIDON, tayph, h5py, Matplotlib, Astropy.
These remain optional for other ExoAtmoSim commands.

**Orbital caution:** The saved `phi` in prepared HDF5 is *true anomaly in
radians*, obtained by solving Kepler's equation from the staged `Tc` and
orbital elements. The current calculation uses the edge-on approximation
`f(Tc)+omega=pi/2`; this may differ from the original orbital estimator
for inclined eccentric systems. Validate against `rvtrace_backend.log`
and the original RV model for accurate non-edge-on cases.

**POSEIDON compatibility:** This stage requires a POSEIDON `cross_correlate`
implementation whose Kp-Vsys velocity expression treats `phi` as true
anomaly. For eccentric orbits, the routine must either accept `ecc` and `omega_rad`
keywords or read the orbital elements from its `dp` configuration using a
correct true-anomaly velocity expression. The CC stage checks for these
interfaces; science validation still requires an injected-signal recovery. The repository cannot patch or
verify your separate editable POSEIDON installation.

**Reference marker:** the expected `Kp` is computed from the orbital masses,
period, inclination and eccentricity. The `Vsys_expected=0` marker denotes
the zero residual-velocity offset in the corrected observer-frame convention. The CCF frame and
velocity sign must be validated before interpreting any detection. Results
are exploratory, not calibrated detection significances.

## Combining VIOLA bands and seeds

Once the individual `cc` runs have produced maps for each requested band
and seed, combine them without re-running POSEIDON:

```bash
exoatmosim combine --planet "WASP-189 b" \
  --instruments VIOLA_J VIOLA_H VIOLA_K \
  --seeds 1 2 3 4 --sysrem 5 --species O2 O3 H2O
```

For every molecule, all twelve instrument–seed raw CCF maps are combined using an S/N-squared weighted mean (equal weights only when S/N FITS are absent). Outputs are
under `data/WASP-189_b/combined/sysrem_5/{species}/` and contain the raw weighted-mean map, Kp and Vsys axes, a POSEIDON-normalized significance figure restricted to nonnegative Kp and labeled Vsys, and `combination.json` provenance.
The command fails if an input is absent or the velocity grids differ. The FITS map remains in raw CCF units; the PNG uses POSEIDON's significance plotting normalization and is not by itself a calibrated false-alarm probability.
Use `--data-dir /data/bibi/VIOLA/data` for simulations stored outside the
repository. `--overwrite` regenerates the combined products.

For eccentric orbits, the Mantis POSEIDON fork accepts `dp` as its first
argument and reads `ecc` and `omega` from the seed config internally.
The CC stage recognizes this variant, checks that its source contains the
eccentric-velocity correction, and does not pass unsupported keyword arguments.
Always verify the model's sign and frame conventions using an injected signal.

## Species templates with fixed atmospheric structure

`exoatmosim template` constructs the full atmospheric composition in POSEIDON,
including the original mean molecular weight and pressure-radius structure.
It then uses `POSEIDON.contributions.spectral_contribution` to isolate the
requested species opacity (with `bulk_species=False` in the opacity
contribution calculation). Other species remain present in the atmospheric
structure but their opacity is not included in the isolated template. The
result is a species-only transit-depth FITS product, not a re-normalized
mixture or a spectrum generated from a reduced-composition atmosphere.

Regenerate any older molecular templates with `exoatmosim template ...
--species H2O --overwrite`, because their atmospheric structures may differ.
This feature requires a POSEIDON installation exposing the documented
`spectral_contribution` interface. Validate it with a known injection on the
installed POSEIDON fork; the Mantis POSEIDON variant has not been executed
within this PR update.

## One-pass full spectrum and molecular templates

Avoid repeatedly loading the same high-resolution opacity grids by running:

```bash
exoatmosim spectrum data/TRAPPIST-1_b/config_viola_k \
  --seed 0 --create-templates --overwrite
```

This constructs the complete atmospheric composition and mean molecular weight
once, reads the opacity grid once, and calls POSEIDON's
`spectral_contribution` once with all non-bulk species. It saves the normal
full-spectrum FITS plus `templates/<species>.fits` for each gas returned by
the contribution routine. This avoids repeating those expensive setup steps,
though each separate contribution still requires radiative-transfer work.

The current implementation rejects a contribution that is nonfinite, missing,
or effectively featureless (peak-to-peak depth ≤ 1e-12). This is a diagnostic
threshold, **not proof of a software error**: cold or otherwise depleted
atmospheres can have genuinely negligible gas-phase abundances. Check the
abundance profile before interpreting a flat contribution. The combined
workflow reduces repeated opacity setup but does not guarantee detectable
molecular features.
