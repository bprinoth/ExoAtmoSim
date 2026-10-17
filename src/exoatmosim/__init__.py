"""ExoAtmoSim: modular exoplanet spectroscopy simulation (foundation)."""
from .config import Config
from .objects import Planet, Star
from .atmosphere import Atmosphere
from .profiles import AbundanceProfile
from .poseidon import TransmissionRun, PoseidonTransmissionModel, TransmissionSpectrum

__all__ = ['Config', 'Planet', 'Star', 'Atmosphere', 'AbundanceProfile', 'TransmissionRun', 'PoseidonTransmissionModel', 'TransmissionSpectrum']

from .stellar_spectra import StellarSpectrum, PhoenixSpectrumProvider
__all__ += ["StellarSpectrum", "PhoenixSpectrumProvider"]

# User-facing pipeline functions; imports of scientific backends are deferred.
from .workflow import (make_config, calculate_atmosphere, calculate_spectrum,
                       simulate_observation, simulate)
__all__ += ['make_config', 'calculate_atmosphere', 'calculate_spectrum',
            'simulate_observation', 'simulate']
