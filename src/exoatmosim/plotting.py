"""Consistent publication-style diagnostics; never modifies global rcParams.

Montserrat is optional. Supply ``font_path`` or set ``EXOATMOSIM_FONT_PATH``;
otherwise use a common sans-serif fallback. No font is bundled.
"""
from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
import warnings


@contextmanager
def publication_style(font_path: str | Path | None = None):
    """Use the project's thick-spine, inward-tick style within a plot block."""
    import matplotlib as mpl
    from matplotlib import font_manager

    filename = font_path if font_path is not None else os.getenv('EXOATMOSIM_FONT_PATH')
    if filename is None:
        default_montserrat = Path('/home/bibi/thesis/static/Montserrat-Bold.ttf')
        if default_montserrat.is_file():
            filename = default_montserrat
    family = 'DejaVu Sans'
    if filename:
        path = Path(filename).expanduser()
        if path.is_file():
            font_manager.fontManager.addfont(str(path))
            family = font_manager.FontProperties(fname=str(path)).get_name()
        else:
            warnings.warn(f'Plot font {path} not found; using {family}', stacklevel=2)
    settings = {
        'font.family': family,
        'font.size': 12,
        'axes.labelsize': 18,
        'xtick.labelsize': 15,
        'ytick.labelsize': 15,
        'legend.fontsize': 18,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 2,
        'ytick.major.width': 2,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.minor.width': 1.5,
        'ytick.minor.width': 1.5,
        'xtick.minor.size': 4,
        'ytick.minor.size': 4,
        'xtick.minor.visible': True,
        'ytick.minor.visible': True,
        'xtick.top': True,
        'xtick.bottom': True,
        'ytick.left': True,
        'ytick.right': True,
        'axes.linewidth': 3,
        'legend.facecolor': 'white',
        'legend.edgecolor': 'white',
        'legend.framealpha': 1,
    }
    with mpl.rc_context(settings):
        yield
