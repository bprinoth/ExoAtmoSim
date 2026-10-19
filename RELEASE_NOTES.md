# Release scope and provenance

This is a minimal constant-S/N foundation, derived from the foundation development source. `ANDES` uses UBV/RIZ/YJH ETC wavelength sampling with repeats and overlaps preserved; no ETC S/N or flux values feed the noise generator. `ANDES_K` is an experimental computational chunking of provisional wavelengths.

The public examples are tested primarily without tellurics. Optional telluric and VIOLA code exists internally but is not claimed as a validated supported workflow in this release.

Excluded from this source distribution: private machine config paths, historical benchmark notes, generated spectral outputs, user data, editable-install metadata, old legacy scripts and tests requiring obsolete benchmark helpers.

- Standardized planet-config exposure duration to `exposure_time` (seconds);
  removed `texp` and `exp_time` configuration aliases. Existing configs must
  rename one of those keys to `exposure_time`.

- Experimental `chemistry kaltenegger` workflow: retain NASA stellar type; map only when selecting 0.0 Ga template, load original Z,T,P / Molecule profiles, add residual N2 and warn about unsupported opacity species. Requires POSEIDON backend validation on mantis.

- v0.9.28.10 Kaltenegger: fill residual N2 from the entire original
  chemistry inventory BEFORE filtering unsupported POSEIDON species. No
  arbitrary omitted-opacity VMR cutoff; emit a warning and retain full
  inventory/omission metadata in the converter. This is an opacity-limited
  approximation and does not represent all input species spectroscopically.

- Kaltenegger profiles are downloaded by users from Zenodo, not redistributed.
  Set `kaltenegger_profiles_dir` in local `config/site_paths.conf`.

## 2026-10-19 diagnostic and NASA compatibility fix
- Blank NASA spectral types are represented as `stellar_type nan` instead of causing an invalid empty config entry.
- Observing-conditions diagnostic plots PWV (mm) when PWV is present, rather than displaying the unused legacy 50% RH series.
- HDF5 telluric grid already uses the PWV series saved to `observing_conditions.ecsv` for each exposure; no RH-to-PWV conversion is performed on this path.

## 2026-10-19 GitHub candidate: VIOLA J/H/K + HDF5 tellurics
- Reuse external high-resolution HDF5 atmospheric transmission for all instruments.
- Private `telluric_grid_dir` is discovered automatically from the planet checkout
  at observing time, or from `EXOATMOSIM_SITE_CONFIG`; no local path is copied
  into public planet configs.
- Preserve native PHOENIX stellar RV handling and VIOLA J/H/K grids.
- VIOLA sky emission is intentionally disabled in this release: S/N optimistic.
- HDF5 grid-builder's approximate PWV calibration remains unvalidated.
- The complete scientific workflow requires a mantis end-to-end validation.

- Documentation now credits https://github.com/sirinrats/viola-etc and describes
  optical-depth telluric interpolation, five site definitions and the
  incomplete/unvalidated Gemini South model generation.

## VIOLA diagnostics restored (October 2026)

Native VIOLA J/H/K observations again write diagnostic figures without
additional radiative transfer: `diagnostics/planet_transmission.png` shows
the input POSEIDON model, and `diagnostics/viola_spectral_comparison.png`
shows the selected in-transit star/tellurics versus planet-injected electron
spectra and their ratio. This release includes the K-band convolution-axis
fix (using `broaden_spectrum`'s returned wavelength axis).

The external HDF5 telluric libraries remain separate downloads referenced
with `telluric_grid_dir` in the private `config/site_paths.conf`. No 1.3 GB
resources are included. The provisional PWV-to-RH library calibration and
zero-sky VIOLA noise approximation are not science-validated.
