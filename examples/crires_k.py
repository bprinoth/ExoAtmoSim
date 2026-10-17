"""Example: simulate TOI-1842 b with CRIRES_K.

Seven K-band spectral orders; demonstrated end-to-end with FastChem/POSEIDON.

This calls the public Python API; scientific implementations live under
src/exoatmosim/. The simulation retrieves or reuses a NASA-generated planet
config, runs chemistry and radiative transfer, retrieves a stellar spectrum,
then synthesizes time-series spectral observations at constant S/N.

Outputs: data/<planet>/simulations/CRIRES_K/seed_0042/
Requires separately installed/configured FastChem, POSEIDON, orbital backends
and the PHOENIX stellar model. No tellurics in the release demonstration.
"""
from exoatmosim import simulate


def main():
    # Change these arguments to configure a different target or noise realization.
    simulate(
        planet='TOI-1842 b',
        instrument='CRIRES_K',
        snr=100.0,  # User-prescribed constant S/N; no ETC noise or throughput.
        seed=42,    # Each seed receives its own immutable config/model snapshot.
        tellurics=False,
    )


if __name__ == '__main__':
    main()
