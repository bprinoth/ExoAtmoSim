"""ExoAtmoSim: modular exoplanet spectroscopy simulation (foundation)."""
from .config import Config
from .objects import Planet, Star
from .atmosphere import Atmosphere
from .profiles import AbundanceProfile
from .poseidon import TransmissionRun, PoseidonTransmissionModel, TransmissionSpectrum

__all__ = ['Config', 'Planet', 'Star', 'Atmosphere', 'AbundanceProfile', 'TransmissionRun', 'PoseidonTransmissionModel', 'TransmissionSpectrum']

from .stellar_spectra import StellarSpectrum, PhoenixSpectrumProvider
__all__ += ["StellarSpectrum", "PhoenixSpectrumProvider"]
