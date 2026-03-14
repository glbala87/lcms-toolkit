# LCMS Toolkit

A comprehensive, cross-platform toolkit for LC-MS (Liquid Chromatography-Mass Spectrometry) data analysis. Provides implementations in **C++**, **Python**, and **Java** with a unified API for mass spectrometry data processing, quantification, identification, and statistical analysis.

## Features

### Core Capabilities
- **Multi-language support**: C++ core library, Python (NumPy/SciPy), and Java implementation
- **Standard file formats**: mzML, mzXML, MGF, MSP, mzTab, imzML, mzIdentML, CSV/TSV
- **Data structures**: Spectrum, Chromatogram, Peak, Feature, and MSExperiment containers
- **Signal processing**: Peak picking, smoothing, baseline correction, centroiding, noise estimation

### Advanced Analysis
- **Isotope detection**: Averagine model, isotope pattern matching, charge state deconvolution
- **Spectral matching**: Cosine, modified cosine, and spectral entropy similarity scoring with library search (MGF/MSP)
- **Label-free quantification**: RT alignment (linear/LOESS), feature alignment, consensus maps, differential analysis
- **Isotope labeling**: SILAC (2/3-plex), TMT (6/10/11/16/18-plex), iTRAQ (4/8-plex), dimethyl labeling
- **Peptide identification**: Database search, fragment ion matching, target-decoy FDR, protein inference
- **RT prediction**: Ridge regression on peptide sequence features (hydrophobicity, composition, charge)
- **Spectrum annotation**: b/y/a/c/x/z ion series, neutral losses (H2O, NH3, phospho), immonium ions

### Statistics & Visualization
- **Multivariate statistics**: PCA (SVD-based), PLS-DA (NIPALS), hierarchical clustering, ANOVA with BH-FDR
- **Volcano plots**: Log2 fold-change vs. -log10 p-value computation
- **Static plots**: Matplotlib-based spectra, chromatograms, heatmaps, mirror plots
- **Interactive plots**: Plotly-based spectrum viewer, chromatogram overlay, feature maps, dashboard (Dash)
- **Report generation**: Styled HTML/PDF reports with embedded figures and tables

### Infrastructure
- **Parallel processing**: Python multiprocessing/threading, C++ ThreadPool, Java ExecutorService
- **Streaming access**: Indexed mzML/mzXML readers with random access by scan number or RT range
- **Memory-mapped arrays**: numpy memmap for large feature matrices and spectral libraries
- **Plugin architecture**: Register custom algorithms, readers, writers, and processing pipelines
- **Cloud/HPC integration**: Dask distributed backend, S3 streaming, Snakemake/Nextflow templates, SLURM/PBS job submission
- **Quality control**: TIC stability, peak capacity, mass accuracy, RT reproducibility metrics
- **CLI tool**: `pylcms` command with subcommands for peaks, convert, quantify, search, QC
- **Java GUI**: JavaFX spectrum/chromatogram viewer with zoom, annotation, and processing tools
- **pybind11 bindings**: C++ performance from Python for core algorithms

## Installation

### Python (pylcms)

```bash
cd python
pip install -e .

# With visualization
pip install -e ".[viz]"

# With interactive plots + dashboard
pip install -e ".[dashboard]"

# With cloud/HPC support
pip install -e ".[cloud]"

# With report generation (PDF)
pip install -e ".[report]"

# With everything
pip install -e ".[full]"

# For development
pip install -e ".[dev]"
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

# Chromatograms
tic = exp.generate_tic()
xic = exp.generate_xic(target_mz=500.0, tolerance=0.5)

# Isotope detection
patterns = pylcms.detect_isotope_patterns(spec)

# Spectral library search
lib = pylcms.SpectralLibrary()
lib.load_mgf("library.mgf")
matches = lib.search(query_spec, precursor_mz=500.0, top_k=10)

# Label-free quantification
from pylcms.quantification import FeatureAlignment, ConsensusMap
aligner = FeatureAlignment()
consensus = aligner.align(feature_maps)
consensus = pylcms.median_normalization(consensus.intensity_matrix)

# TMT reporter ion quantification
from pylcms.labeling import extract_reporter_ions, LabelingStrategy
quant = extract_reporter_ions(ms2_spec, LabelingStrategy.TMT10)
print(quant.channel_intensities)

# PCA
from pylcms.statistics import pca
result = pca(consensus, n_components=3)
print(f"PC1 explains {result.explained_variance_ratio[0]:.1%}")

# Spectrum annotation
from pylcms.annotation import annotate_spectrum, IonType
ann = annotate_spectrum(ms2_spec, "PEPTIDER", precursor_charge=2)
print(f"Coverage: {ann.coverage:.1%}, Matched: {ann.n_matched} peaks")

# Generate report
from pylcms.reporting import ReportBuilder, ReportConfig
builder = ReportBuilder(ReportConfig(title="My Analysis"))
builder.add_summary(n_spectra=5000, n_features=1200)
builder.add_qc_section(qc_metrics)
builder.save_html("report.html")

# Visualization
pylcms.plot_spectrum(spec, mz_range=(100, 500))
```

### CLI

```bash
# File info
pylcms info sample.mzML

# Peak picking
pylcms peaks sample.mzML -o peaks.csv --snr 3.0

# Convert formats
pylcms convert sample.mzML --format mztab -o output.mztab

# Run QC
pylcms qc sample.mzML -o qc_report.json

# Database search
pylcms search sample.mzML --fasta proteins.fasta -o results.csv
```

### C++

```cpp
#include <lcms/lcms.hpp>

// Create and analyze spectrum
lcms::Spectrum spec(mz, intensity);
std::cout << "TIC: " << spec.tic() << std::endl;

// Isotope detection
auto patterns = lcms::detectIsotopePatterns(mz, intensity);

// Spectral matching
double score = lcms::cosineSimilarity(mz1, int1, mz2, int2, 0.01);

// Parallel processing
auto results = lcms::parallelPeakPicking(spectra, params);
```

### Java

```java
import org.lcms.core.*;
import org.lcms.gui.LCMSViewer;

// Create spectrum
Spectrum spec = new Spectrum(mz, intensity);
System.out.println("TIC: " + spec.getTotalIonCurrent());

// Spectral library search
SpectralLibrary lib = new SpectralLibrary();
lib.loadMGF("library.mgf");
List<SpectralMatch> matches = lib.search(querySpec, 500.0, 10);

// Launch GUI viewer
LCMSViewer.launch(args);
```

### Plugin System

```python
from pylcms.plugins import PluginRegistry, register_as, ProcessingPipeline

# Register a custom algorithm
@register_as("algorithm", "my_peak_picker")
def my_peak_picker(spectrum, threshold=100):
    # Custom peak picking logic
    ...

# Build a processing pipeline
pipeline = ProcessingPipeline()
pipeline.add_step("smooth", window_size=5)
pipeline.add_step("my_peak_picker", threshold=50)
results = pipeline.run(spectrum)

# Discover plugins from installed packages
registry = PluginRegistry.instance()
registry.discover_plugins("pylcms_plugins")
```

### Cloud/HPC

```python
from pylcms.cloud import DaskBackend, S3FileHandler, generate_snakemake_workflow

# Distributed processing with Dask
backend = DaskBackend(n_workers=16)
results = backend.map(process_file, file_list)

# Stream from S3
s3 = S3FileHandler(bucket="my-lcms-data")
files = s3.list_files(prefix="raw/", suffix=".mzML")

# Generate workflow templates
generate_snakemake_workflow(output_path="Snakefile")
```

## Project Structure

```
lcms-toolkit/
├── core/                          # C++ core library
│   ├── include/lcms/              # Headers
│   │   ├── algorithms/            # Peak picking, isotope, matching, quantification
│   │   ├── io/                    # Indexed readers
│   │   └── parallel.hpp           # Thread pool
│   ├── src/                       # Implementations
│   ├── bindings/                  # pybind11 Python bindings
│   ├── tests/                     # Catch2 unit tests
│   └── CMakeLists.txt
├── python/                        # Python package
│   ├── pylcms/
│   │   ├── algorithms.py          # Signal processing
│   │   ├── isotope.py             # Isotope detection & deconvolution
│   │   ├── spectral_matching.py   # Library search & similarity
│   │   ├── quantification.py      # LFQ, alignment, normalization
│   │   ├── labeling.py            # SILAC, TMT, iTRAQ, dimethyl
│   │   ├── statistics.py          # PCA, PLS-DA, ANOVA, clustering
│   │   ├── identification.py      # Database search, FDR, protein inference
│   │   ├── rt_prediction.py       # ML-based RT prediction
│   │   ├── annotation.py          # Fragment ion annotation
│   │   ├── streaming.py           # Indexed file access
│   │   ├── parallel.py            # Batch processing
│   │   ├── qc.py                  # Quality control metrics
│   │   ├── formats.py             # MGF, imzML, mzIdentML, CSV I/O
│   │   ├── visualization.py       # Matplotlib plots
│   │   ├── interactive_viz.py     # Plotly/Dash interactive plots
│   │   ├── reporting.py           # HTML/PDF report generation
│   │   ├── memmap.py              # Memory-mapped large arrays
│   │   ├── plugins.py             # Plugin architecture
│   │   ├── cloud.py               # Dask, S3, Snakemake, HPC
│   │   └── cli.py                 # Command-line interface
│   ├── tests/                     # pytest unit tests
│   └── setup.py
├── java/                          # Java implementation
│   ├── src/main/java/org/lcms/
│   │   ├── core/                  # Spectrum, SpectralLibrary, Quantification
│   │   ├── gui/                   # JavaFX viewer
│   │   ├── io/                    # Indexed readers
│   │   └── tools/                 # Batch/parallel processing
│   ├── src/test/                  # JUnit 5 tests
│   └── pom.xml
└── examples/                      # Example scripts
```

## API Reference

### Data Structures

| Class | Description |
|-------|-------------|
| `Spectrum` | Mass spectrum with m/z and intensity arrays |
| `Chromatogram` | Intensity vs. retention time data |
| `Peak` / `PeakList` | Detected peak with m/z, RT, area, FWHM |
| `Feature` / `FeatureMap` | 2D feature spanning multiple spectra |
| `MSExperiment` | Container for a complete LC-MS run |
| `ConsensusMap` | Aligned feature intensity matrix across samples |
| `IsotopePattern` | Detected isotope envelope |
| `SpectralLibrary` / `SpectralMatch` | Library search infrastructure |
| `PeptideSpectrumMatch` / `ProteinGroup` | Identification results |
| `RTPrediction` | RT prediction result with confidence |

### Algorithms

| Function | Description |
|----------|-------------|
| `pick_peaks()` | Detect peaks with SNR filtering |
| `smooth_spectrum()` | Gaussian, Savitzky-Golay, or moving average |
| `correct_baseline()` | SNIP, top-hat, or rolling ball methods |
| `centroid_spectrum()` | Convert profile to centroided |
| `detect_isotope_patterns()` | Isotope envelope detection |
| `deconvolute_spectrum()` | Charge state deconvolution |
| `cosine_similarity()` | Spectral cosine similarity |
| `modified_cosine_similarity()` | Modified cosine with precursor shift |
| `spectral_entropy_similarity()` | Entropy-based similarity |
| `annotate_spectrum()` | Fragment ion annotation |

### Quantification & Statistics

| Function | Description |
|----------|-------------|
| `median_normalization()` | Column median normalization |
| `quantile_normalization()` | Quantile normalization |
| `tic_normalization()` | Total ion current normalization |
| `extract_reporter_ions()` | TMT/iTRAQ reporter extraction |
| `find_silac_pairs()` | SILAC pair matching |
| `pca()` | Principal Component Analysis |
| `plsda()` | PLS Discriminant Analysis |
| `anova()` | One-way ANOVA with BH-FDR |
| `volcano_data()` | Volcano plot data computation |
| `target_decoy_fdr()` | FDR control via target-decoy |

### I/O

| Function | Description |
|----------|-------------|
| `load_mzml()` / `load_mzxml()` | Load standard MS files |
| `load_mgf()` / `save_mgf()` | MGF spectral library format |
| `save_mztab()` | Export to mzTab |
| `save_mzidentml()` | Export identification results |
| `load_feature_table()` / `save_feature_table()` | CSV/TSV feature tables |
| `IndexedMzMLReader` | Random-access indexed reading |

## Requirements

### Python
- Python >= 3.8
- NumPy >= 1.20
- SciPy >= 1.7
- Optional: Matplotlib, Plotly, Dash, Pandas, WeasyPrint, Dask, boto3

### C++
- C++17 compatible compiler
- CMake >= 3.14
- pybind11 (for Python bindings)
- Catch2 (for testing)

### Java
- Java 11+
- Maven 3.6+
- JavaFX (for GUI viewer)

## Testing

```bash
# Python
cd python && pip install -e ".[dev]" && pytest tests/ -v

# C++
cd core/build && cmake .. -DBUILD_TESTS=ON && make && ctest

# Java
cd java && mvn test
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## Citation

If you use this toolkit in your research, please cite:

```bibtex
@software{lcms_toolkit,
  title = {LCMS Toolkit: A Cross-Platform LC-MS Data Analysis Library},
  year = {2024},
  url = {https://github.com/glbala87/lcms-toolkit}
}
```
