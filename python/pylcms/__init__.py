"""
PyLCMS - Python LC-MS Data Analysis Toolkit

A Python interface for LC-MS (Liquid Chromatography-Mass Spectrometry) data
analysis, providing NumPy integration and high-level workflows.

Example:
    >>> import pylcms
    >>> exp = pylcms.load_mzml("sample.mzML")
    >>> spec = exp.spectrum(0)
    >>> peaks = pylcms.pick_peaks(spec)
"""

__version__ = "1.0.0"
__author__ = "LCMS Toolkit Contributors"

from .spectrum import Spectrum, SpectrumType, Polarity
from .chromatogram import Chromatogram, ChromatogramType
from .peak import Peak, PeakList
from .feature import Feature, FeatureMap
from .experiment import MSExperiment
from .io import load_mzml, load_mzxml, save_mztab
from .algorithms import (
    pick_peaks,
    centroid_spectrum,
    smooth_spectrum,
    correct_baseline,
    estimate_noise,
)
from .visualization import (
    plot_spectrum,
    plot_chromatogram,
    plot_heatmap,
    plot_peaks,
)

__all__ = [
    # Version
    "__version__",
    # Data structures
    "Spectrum",
    "SpectrumType",
    "Polarity",
    "Chromatogram",
    "ChromatogramType",
    "Peak",
    "PeakList",
    "Feature",
    "FeatureMap",
    "MSExperiment",
    # I/O
    "load_mzml",
    "load_mzxml",
    "save_mztab",
    # Algorithms
    "pick_peaks",
    "centroid_spectrum",
    "smooth_spectrum",
    "correct_baseline",
    "estimate_noise",
    # Visualization
    "plot_spectrum",
    "plot_chromatogram",
    "plot_heatmap",
    "plot_peaks",
]


def version_info():
    """Return detailed version information."""
    import sys
    info = {
        "pylcms_version": __version__,
        "python_version": sys.version,
    }
    try:
        import numpy as np
        info["numpy_version"] = np.__version__
    except ImportError:
        pass
    try:
        import scipy
        info["scipy_version"] = scipy.__version__
    except ImportError:
        pass
    return info
