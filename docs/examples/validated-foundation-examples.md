# Two foundation science examples

Both use the planet directory's authoritative configuration file named config. The commands and final paths will be confirmed after importing the foundation code onto this branch.

## TOI-1842 b — eccentric, H2O, CRIRES+ K

Physics: FastChem equilibrium chemistry → POSEIDON H2O transmission spectrum → independent CRIRES+ K mock.

Required checks: eccentricity nonzero; seven wavelength/order FITS with consistent shapes; observer-frame stellar and planetary RVs; observing times; atmospheric and spectral plots; reproducible noise seed. A frozen tayph comparison may be used only as a development benchmark.

## WASP-189 b — circular, Fe, ESPRESSO

Physics: FastChem equilibrium chemistry → POSEIDON Fe transmission → independent ESPRESSO mock with clearly specified constant-SNR assumptions.

Required checks: eccentricity zero, physically supported Fe opacity, wavelength coverage including orbital RV margins, correct observer-frame velocities and transit weights, same seed reproduces noise, instrument/seed-separated outputs.

**Caveat:** the current packaged examples/wasp189b_pipeline.py calls its ESPRESSO mock approximate. The GitHub release documentation must not imply detector-level verification until the full pipeline passes on mantis.

## Common output

data/<planet>/simulations/<instrument>/seed_<seed>/ with observation arrays, shared wavelengths, conditions, diagnostics and a run manifest. The underlying planet/star components are reusable when their physical input hashes match.