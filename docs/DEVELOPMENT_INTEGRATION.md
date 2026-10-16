# Foundation integration — development staging

This is a staging branch, not a release. The tested development source is not yet imported into GitHub.

## Delivery plan

1. Import the modular source, packaging, tests, and measured instrument resources from the v0.9.19 development snapshot, auditing dependencies, personal paths, and generated files.
2. Restore a complete eccentric TOI-1842 b / H2O / CRIRES+ K example: FastChem → POSEIDON → independent observation simulation.
3. Restore a complete circular WASP-189 b / Fe / ESPRESSO example with the same workflow. Its currently packaged mock labels itself approximate; do not advertise full detector fidelity before validation.
4. Run both examples on mantis; publish exact configs (each named config), dependency instructions, commands, and output validation.
5. Test clean installation and scientific invariants before preparing a release tag.

## Known integration blockers

- The v0.9.19 CRIRES+ full-pipeline runner checks outputs in independent_mock, but the updated independent runner defaults to simulations/<instrument>/seed_<seed>. Reconcile the output contract.
- Some tools and benchmarks still refer to tayph; eliminate mandatory tayph runtime and personal absolute paths from the release.
- Preserve reproducible instrument/seed outputs, manifest hashes, air/vacuum conventions, observer-frame RV, and diagnostic directories.
- VIOLA is intentionally out of scope for the initial validated-example integration; track its wavelength-coverage failure separately.

## Acceptance

- [ ] Actual foundation code reviewed and committed
- [ ] Full TOI-1842 b example run and validated on mantis
- [ ] Full WASP-189 b example run and validated on mantis
- [ ] Release requirements and clean-install smoke test pass
- [ ] No generated spectral FITS, caches, or intermediate ETC directories checked into Git