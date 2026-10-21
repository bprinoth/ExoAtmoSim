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
anomaly. For eccentric orbits, the function must expose `ecc` and
`omega_rad` keyword parameters, or the CC stage raises an error rather
than silently assuming circular motion. The repository cannot patch or
verify your separate editable POSEIDON installation.

**Reference marker:** `Kp_expected=0` and `Vsys_expected=0` are purely
plotting placeholders. They are not physical estimates. The CCF frame and
velocity sign must be validated before interpreting any detection. Results
are exploratory, not calibrated detection significances.
