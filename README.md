# LCMS Toolkit

A comprehensive, cross-platform toolkit for LC-MS (Liquid Chromatography-Mass Spectrometry) data analysis. Provides implementations in C++, Python, and Java with a unified API for working with mass spectrometry data.

## Features

- **Multi-language support**: C++ core library, Python bindings, and Java implementation
- **Standard file formats**: Read mzML and mzXML files
- **Data structures**: Spectrum, Chromatogram, Peak, Feature, and MSExperiment containers
- **Signal processing**: Peak picking, smoothing, baseline correction, noise estimation
- **Visualization**: Plot spectra, chromatograms, heatmaps, and mirror plots
- **Workflows**: Automated peak picking and feature detection pipelines

## Installation

### Python (pylcms)

```bash
cd python
pip install -e .

# With visualization support
pip install -e ".[viz]"

# With all optional dependencies
pip install -e ".[full]"
```

### C++ (lcms_core)

```bash
cd core
mkdir build && cd build
cmake ..
make
```

### Java (lcms-java)

```bash
cd java
mvn package
```

## Quick Start

### Python

```python
import pylcms

# Load an mzML file
exp = pylcms.load_mzml("sample.mzML")

# Access spectra
spec = exp.spectrum(0)
print(f"Base peak: m/z {spec.base_peak_mz:.4f}")

# Pick peaks
peaks = pylcms.pick_peaks(spec, min_snr=5)
print(f"Found {len(peaks)} peaks")

# Generate chromatograms
tic = exp.generate_tic()
xic = exp.generate_xic(target_mz=500.0, tolerance=0.5)

# Visualization
from pylcms.visualization import plot_spectrum
plot_spectrum(spec, mz_range=(100, 500))
```

### C++

```cpp
#include <lcms/lcms.hpp>

// Create spectrum
std::vector<double> mz = {100.0, 200.0, 300.0};
std::vector<double> intensity = {1000.0, 5000.0, 2000.0};
lcms::Spectrum spec(mz, intensity);

// Access data
std::cout << "TIC: " << spec.tic() << std::endl;
std::cout << "Base peak: " << spec.basePeakMz() << std::endl;
```

### Java

```java
import org.lcms.core.*;

// Create spectrum
double[] mz = {100.0, 200.0, 300.0};
double[] intensity = {1000.0, 5000.0, 2000.0};
Spectrum spec = new Spectrum(mz, intensity);

System.out.println("TIC: " + spec.getTic());
System.out.println("Base peak: " + spec.getBasePeakMz());
```

## Project Structure

```
lcms-toolkit/
├── core/                   # C++ core library
│   ├── include/lcms/       # Header files
│   ├── src/                # Source files
│   ├── tests/              # Unit tests (Catch2)
│   └── CMakeLists.txt
├── python/                 # Python package
│   ├── pylcms/             # Python modules
│   └── setup.py
├── java/                   # Java implementation
│   ├── src/main/java/      # Java sources
│   └── pom.xml
└── examples/               # Example scripts
```

## API Overview

### Data Structures

| Class | Description |
|-------|-------------|
| `Spectrum` | Mass spectrum with m/z and intensity arrays |
| `Chromatogram` | Intensity vs retention time data |
| `Peak` | Detected peak with properties (m/z, RT, area, FWHM) |
| `PeakList` | Collection of peaks with filtering methods |
| `Feature` | 2D feature spanning multiple spectra |
| `FeatureMap` | Collection of features with spatial indexing |
| `MSExperiment` | Container for a complete LC-MS run |

### Algorithms

| Function | Description |
|----------|-------------|
| `pick_peaks()` | Detect peaks with SNR filtering |
| `smooth_spectrum()` | Gaussian, Savitzky-Golay, or moving average |
| `correct_baseline()` | SNIP, top-hat, or rolling ball methods |
| `centroid_spectrum()` | Convert profile to centroided |
| `estimate_noise()` | MAD, percentile, or STD methods |

### I/O Functions

| Function | Description |
|----------|-------------|
| `load_mzml()` | Load mzML files |
| `load_mzxml()` | Load mzXML files |
| `save_mztab()` | Export features to mzTab format |

## Requirements

### Python
- Python >= 3.8
- NumPy >= 1.20
- SciPy >= 1.7
- Matplotlib >= 3.4 (optional, for visualization)
- Pandas >= 1.3 (optional, for DataFrame support)

### C++
- C++17 compatible compiler
- CMake >= 3.14
- pugixml (optional, for XML parsing)
- zlib (optional, for compression)
- Eigen3 (optional, for linear algebra)
- Catch2 (for testing)

### Java
- Java 11+
- Maven 3.6+

## License

MIT License - see [LICENSE](LICENSE) for details.


