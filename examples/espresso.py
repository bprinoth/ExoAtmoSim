"""Example: simulate WASP-189 b with ESPRESSO.

170 optical echelle orders. Configure Fe chemistry/opacity explicitly in the planet config.

This calls the public Python API; scientific implementations live under
src/exoatmosim/. The simulation retrieves or reuses a NASA-generated planet
config, runs chemistry and radiative transfer, retrieves a stellar spectrum,
then synthesizes time-series spectral observations at constant S/N.

Outputs: data/<planet>/simulations/ESPRESSO/seed_0042/
Requires separately installed/configured FastChem, POSEIDON, orbital backends
and the PHOENIX stellar model. No tellurics in the release demonstration.
"""
from exoatmosim import simulate


def main():
    # Change these arguments to configure a different target or noise realization.
    simulate(
        planet='WASP-189 b',
        instrument='ESPRESSO',
        snr=100.0,  # User-prescribed constant S/N; no ETC noise or throughput.
        seed=42,    # Each seed receives its own immutable config/model snapshot.
        tellurics=False,
    )


if __name__ == '__main__':
    main()
