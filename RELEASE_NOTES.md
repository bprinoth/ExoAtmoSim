# Release scope and provenance

This is a minimal constant-S/N foundation, derived from the foundation development source. `ANDES` uses UBV/RIZ/YJH ETC wavelength sampling with repeats and overlaps preserved; no ETC S/N or flux values feed the noise generator. `ANDES_K` is an experimental computational chunking of provisional wavelengths.

The public examples are tested primarily without tellurics. Optional telluric and VIOLA code exists internally but is not claimed as a validated supported workflow in this release.

Excluded from this source distribution: private machine config paths, historical benchmark notes, generated spectral outputs, user data, editable-install metadata, old legacy scripts and tests requiring obsolete benchmark helpers.

- Standardized planet-config exposure duration to `exposure_time` (seconds);
  removed `texp` and `exp_time` configuration aliases. Existing configs must
  rename one of those keys to `exposure_time`.

## Development v0.9.28.11
- Optional Kaltenegger chemistry, external Zenodo profiles, no data redistribution.
- NASA spectral metadata preserved; N2 residual filled before unsupported-opacity filtering.
- FastChem unchanged; the approximate Earth-analogue geological clock is disclosed in README.
