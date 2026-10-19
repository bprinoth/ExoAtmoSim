# VIOLA J/H/K and external HDF5 telluric-grid integration (draft)

This is a **staging PR**, not ready to merge until the tested release source has been uploaded to this feature branch.

Scope of the tested release to stage:
- Native VIOLA J/H/K wavelength sampling and PHOENIX-based stellar spectra with per-exposure Doppler shifts.
- External HDF5 telluric grid reader and interpolation driven by airmass and sampled PWV.
- External scientific-data paths resolved at runtime from private `config/site_paths.conf` without committing machine-specific paths.
- Updated configuration examples, diagnostics, and regression tests.
- Existing Kaltenegger chemistry preserved.

Important: generated grids are external data, not files to commit to GitHub. The source library must remain scientifically labeled as using an approximate PWV-to-relative-humidity mapping until calibrated.

Review checklist once code has been pushed:
- [ ] Inspect the full source diff against main; no private configs or simulation outputs.
- [ ] Run `python -m unittest discover -s tests -q`.
- [ ] Check ESPRESSO/WASP-189 b with HDF5 tellurics.
- [ ] Check VIOLA_K/TRAPPIST-1 b and convolution array-length handling.
- [ ] Verify site-path discovery both from checkout root and planet-specific config paths.

## Source attribution and upstream reference

The native VIOLA implementation adapts the instrument assumptions and photon/detector-noise formulas from the **VIOLA Exposure Time Calculator (VIOLA ETC)**:

- **Official upstream repository:** https://github.com/sirinrats/viola-etc
- Relevant upstream modules: [instrument configuration](https://github.com/sirinrats/viola-etc/blob/main/viola_etc/models.py), [signal and noise equations](https://github.com/sirinrats/viola-etc/blob/main/viola_etc/core/runner.py), and [SkyCalc grid handling](https://github.com/sirinrats/viola-etc/blob/main/viola_etc/observing/sky_io.py).
- ExoAtmoSim's reference telescope diameter is **8 m**, explicitly overriding the upstream ETC's 2 m code default; do not suggest that this reference value was sourced from the ETC.
- Any future VIOLA SkyCalc/other datasets should retain an explicit link to the upstream repo in their documentation, with applicable licence/provenance reviewed before redistribution.

**Attribution placement when the integration code is staged:** include the upstream URL in the project's README section for VIOLA, VIOLA backend and noise-model module docstrings, scientific assumptions/limitations documentation, and the PR description. Cite the specific upstream module next to adapted equations. Avoid implying formal endorsement by the VIOLA ETC authors.
