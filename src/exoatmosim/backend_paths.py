"""Optional source-checkout locations for external backends.

This is an explicit compatibility mechanism for legacy research installations.
Editable pip installations are generally preferable. Python paths are process-wide.
"""
from __future__ import annotations

from pathlib import Path
import os
import sys

from .config import Config


def configure_backend_path(config: Config, *, key: str, package: str,
                           config_dir: str | Path | None = None) -> Path | None:
    """Expose a checkout root to Python when the config explicitly supplies one.

    E.g. tayph_path /home/user/GitHub/tayph points to the checkout root
    containing the `tayph/` Python package, not the package directory itself.
    No path changes occur if `key` is absent.
    """
    if key not in config:
        return None
    raw = config[key].strip()
    if not raw:
        raise ValueError(f'Config key {key!r} is empty')
    root = Path(raw).expanduser()
    if not root.is_absolute():
        root = (Path(config_dir) if config_dir is not None else Path.cwd()) / root
    root = root.resolve()
    if not (root / package / '__init__.py').is_file():
        raise FileNotFoundError(
            f'{key}={root} is not a checkout root containing {package}/__init__.py'
        )
    # Avoid silently using an already-imported, different checkout.
    loaded = sys.modules.get(package)
    if loaded is not None:
        location = getattr(loaded, '__file__', None)
        if location is not None and not Path(location).resolve().is_relative_to(root / package):
            raise RuntimeError(
                f'{package} is already imported from {location}; cannot switch to {root}. '
                'Start a fresh Python process for a different backend path.'
            )
    path = str(root)
    if path in sys.path:
        sys.path.remove(path)
    sys.path.insert(0, path)
    return root


def configure_poseidon_environment(config: Config) -> dict[str, str]:
    """Set optional POSEIDON data directories before importing POSEIDON.

    Explicit config entries override environment settings for this process;
    absent entries leave the caller's environment untouched. The paths must
    exist, so typos fail before an expensive opacity calculation starts.
    """
    keys = {
        'poseidon_input_data': 'POSEIDON_input_data',
        'pysyn_cdbs': 'PYSYN_CDBS',
    }
    configured = {}
    for config_key, env_key in keys.items():
        if config_key not in config:
            continue
        raw = config[config_key].strip()
        if not raw:
            raise ValueError(f'{config_key} must not be empty')
        path = Path(raw).expanduser().resolve()
        if not path.is_dir():
            raise FileNotFoundError(f'{config_key} directory does not exist: {path}')
        configured[env_key] = str(path)
    # Validate all configured paths before changing any environment values.
    os.environ.update(configured)
    return configured
