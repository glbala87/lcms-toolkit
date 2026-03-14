# LCMS Toolkit

A comprehensive, cross-platform toolkit for **LC-MS (Liquid Chromatography-Mass Spectrometry)** data analysis. Provides implementations in C++, Python, and Java with a unified API for mass spectrometry data processing, quantification, identification, and statistical analysis.

---

## Quick Start (No Data Files Needed)

```bash
git clone https://github.com/glbala87/lcms-toolkit.git
cd lcms-toolkit
pip install -e python/
python run_demo.py
```

This runs **13 interactive demos** using synthetic data — spectrum operations, peak picking, isotope detection, spectral matching, quantification, TMT/SILAC labeling, PCA/ANOVA statistics, peptide identification, RT prediction, spectrum annotation, HTML reporting, plugin pipelines, memory-mapped arrays, and cloud/HPC workflows.

Run a specific demo:

```bash
python run_demo.py spectrum        # Spectrum & peak picking
python run_demo.py matching        # Spectral similarity & library search
python run_demo.py statistics      # PCA, PLS-DA, ANOVA
python run_demo.py annotation      # b/y ion annotation
python run_demo.py --list          # See all 13 demos
```

Or use the one-step setup script:

```bash
./setup.sh
python run_demo.py
```

---

## Analyze Real Data

Use `analyze.py` for a complete analysis pipeline on real mzML/mzXML files:

```bash
# Basic analysis (peaks, TIC, QC metrics)
python analyze.py sample.mzML

# Full analysis with features, report, and custom SNR
python analyze.py sample.mzML --features --report --snr 5

# Extract ion chromatograms at specific m/z values
python analyze.py sample.mzML --xic 500.0 750.3 --xic-tolerance 0.5

# Search against a spectral library (MGF)
python analyze.py sample.mzML --search library.mgf --min-score 0.7

# Annotate MS2 spectra with a peptide sequence
python analyze.py sample.mzML --annotate PEPTIDER --charge 2

# Multi-file analysis with report
python analyze.py sample1.mzML sample2.mzML sample3.mzML --features --report

# Output to a custom directory in TSV format
python analyze.py sample.mzML -o my_results/ --format tsv

# Filter by MS level and RT range
python analyze.py sample.mzML --ms-level 1 --rt-range 60 600
```

**Output files** (saved to `lcms_results/` by default):

| File | Contents |
|------|----------|
| `*_peaks.csv` | All detected peaks (m/z, intensity, RT, SNR, FWHM, area) |
| `*_features.csv` | Linked features across scans (m/z, RT range, volume, quality) |
| `*_tic.csv` | Total ion chromatogram |
| `*_xic_*.csv` | Extracted ion chromatograms |
| `*_matches.csv` | Spectral library search results |
| `*_annotations.csv` | Fragment ion annotations |
| `*_qc.json` | Quality control metrics |
| `analysis_report.html` | Full HTML report (with `--report`) |

See all options: `python analyze.py --help`

---

## Features

### Core Capabilities
- **Multi-language support**: C++ core library, Python (NumPy/SciPy), and Java implementation
- **Standard file formats**: mzML, mzXML, MGF, MSP, mzTab, imzML, mzIdentML, CSV/TSV
- **Data structures**: Spectrum, Chromatogram, Peak, Feature, and MSExperiment containers
- **Signal processing**: Peak picking, smoothing, baseline correction, centroiding, noise estimation

### Advanced Analysis
- **Isotope detection**: Averagine model, isotope pattern matching, charge state deconvolution
- **Spectral matching**: Cosine, modified cosine, and spectral entropy similarity with library search (MGF/MSP)
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

---

## Installation

### Python (pylcms)

```bash
cd python
pip install -e .

# With visualization support
pip install -e ".[viz]"

# With interactive plots + dashboard
pip install -e ".[dashboard]"

# With cloud/HPC support (Dask, S3)
pip install -e ".[cloud]"

# With PDF report generation
pip install -e ".[report]"

# With everything
pip install -e ".[full]"

# For development (pytest, black, mypy, flake8)
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

---

## Usage Examples

### Python API

```python
import pylcms
import numpy as np

# Load data
exp = pylcms.load_mzml("sample.mzML")
spec = exp.spectrum(0)
print(f"Base peak: m/z {spec.base_peak_mz:.4f}, TIC: {spec.tic:.0f}")

# Peak picking pipeline
smoothed = pylcms.smooth_spectrum(spec, method="gaussian", window_size=5)
corrected = pylcms.correct_baseline(smoothed, method="snip")
peaks = pylcms.pick_peaks(corrected, min_snr=5)
print(f"Found {len(peaks)} peaks")

# Chromatograms
tic = exp.generate_tic()
xic = exp.generate_xic(target_mz=500.0, tolerance=0.5)

# Isotope detection
patterns = pylcms.detect_isotope_patterns(spec)

# Spectral library search
lib = pylcms.SpectralLibrary()
lib.load_mgf("library.mgf")
matches = lib.search(query_mz, query_int, query_precursor_mz=500.0, top_n=10)

# Label-free quantification
from pylcms.quantification import ConsensusMap, DifferentialAnalysis
consensus = ConsensusMap(intensity_matrix=matrix, feature_ids=ids, sample_names=names)
normalized = pylcms.median_normalization(consensus.intensity_matrix)

# TMT reporter ion quantification
from pylcms.labeling import extract_reporter_ions, LabelingStrategy
quant = extract_reporter_ions(ms2_spec, LabelingStrategy.TMT10)
print(quant.channel_intensities)

# PCA
from pylcms.statistics import pca
result = pca(consensus, n_components=3)
print(f"PC1 explains {result.explained_variance_ratio[0]:.1%}")

# Spectrum annotation
from pylcms.annotation import annotate_spectrum
ann = annotate_spectrum(ms2_spec, "PEPTIDER", precursor_charge=2)
print(f"Coverage: {ann.coverage:.0%}, Matched: {ann.n_matched} peaks")

# Peptide identification
from pylcms.identification import calculate_peptide_mass, generate_theoretical_fragments
mass = calculate_peptide_mass("PEPTIDER")
fragments = generate_theoretical_fragments("PEPTIDER", charge=1)

# RT prediction
predictor = pylcms.RTPredictor()
predictor.train(peptide_list, rt_list)
predicted_rt = predictor.predict_single("PEPTIDER")

# Generate HTML report
from pylcms.reporting import ReportBuilder, ReportConfig
builder = ReportBuilder(ReportConfig(title="My Analysis"))
builder.add_summary(n_spectra=5000, n_ms1=3000, n_ms2=2000, n_features=1200)
builder.save_html("report.html")

# Visualization
pylcms.plot_spectrum(spec, mz_range=(100, 500))
```

### CLI Tool

```bash
# File info
pylcms info sample.mzML

# Peak picking
pylcms peaks sample.mzML -o peaks.csv --snr 3.0

# Convert formats
pylcms convert sample.mzML --format mztab -o output.mztab

# Extract ion chromatogram
pylcms xic sample.mzML --mz 500.0 --tolerance 0.5

# Label-free quantification (multiple files)
pylcms quantify sample1.mzML sample2.mzML sample3.mzML -o consensus.csv

# Spectral library search
pylcms search sample.mzML --library library.mgf --min-score 0.7

# Quality control
pylcms qc sample1.mzML sample2.mzML -o qc_report.json
```

### Plugin System

```python
from pylcms.plugins import PluginRegistry, register_as, ProcessingPipeline

@register_as("processor", "my_filter")
def my_filter(data, threshold=100):
    return data[data > threshold]

pipeline = ProcessingPipeline()
pipeline.add_step("my_filter", threshold=50)
result = pipeline.run(data)

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

# Generate Snakemake/Nextflow workflow templates
generate_snakemake_workflow(output_path="Snakefile")

# Generate SLURM/PBS job scripts
from pylcms.cloud import HPCJobSubmitter
submitter = HPCJobSubmitter(scheduler="slurm")
script = submitter.generate_script("pylcms peaks sample.mzML -o peaks.csv",
                                   cpus=8, memory="16G", time="1:00:00")
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

// Launch JavaFX GUI viewer
LCMSViewer.launch(args);
```

---

## Project Structure

```
lcms-toolkit/
├── analyze.py                     # Analyze real mzML/mzXML files (full pipeline)
├── run_demo.py                    # One-command demo (13 features, no data needed)
├── setup.sh                       # One-step install script
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
│   ├── tests/                     # pytest unit tests (10 test files)
│   └── setup.py
├── java/                          # Java implementation
│   ├── src/main/java/org/lcms/
│   │   ├── core/                  # Spectrum, SpectralLibrary, Quantification
│   │   ├── gui/                   # JavaFX viewer (LCMSViewer)
│   │   ├── io/                    # Indexed mzML/mzXML readers
│   │   └── tools/                 # Batch/parallel processing
│   ├── src/test/                  # JUnit 5 tests
│   └── pom.xml
└── examples/                      # Example scripts
```

---

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
| `RTPrediction` / `RTPredictor` | RT prediction model and results |
| `ReportBuilder` | HTML/PDF report generator |
| `MemmapMatrix` / `MemmapSpectraStore` | Disk-backed large arrays |
| `PluginRegistry` / `ProcessingPipeline` | Plugin system |

### Algorithms

| Function | Description |
|----------|-------------|
| `pick_peaks()` | Detect peaks with SNR filtering |
| `smooth_spectrum()` | Gaussian, Savitzky-Golay, or moving average |
| `correct_baseline()` | SNIP, top-hat, or rolling ball methods |
| `centroid_spectrum()` | Convert profile to centroided |
| `estimate_noise()` | MAD, percentile, or STD noise estimation |
| `detect_isotope_patterns()` | Isotope envelope detection |
| `deconvolute_spectrum()` | Charge state deconvolution |
| `cosine_similarity()` | Spectral cosine similarity |
| `modified_cosine_similarity()` | Modified cosine with precursor shift |
| `spectral_entropy_similarity()` | Entropy-based similarity |
| `annotate_spectrum()` | Fragment ion annotation (b/y/a/c/x/z) |
| `compute_fragment_ions()` | Theoretical fragment m/z computation |

### Quantification & Statistics

| Function | Description |
|----------|-------------|
| `median_normalization()` | Column median normalization |
| `quantile_normalization()` | Quantile normalization |
| `tic_normalization()` | Total ion current normalization |
| `extract_reporter_ions()` | TMT/iTRAQ reporter extraction |
| `find_silac_pairs()` | SILAC pair matching |
| `normalize_reporter_intensities()` | Reporter ion normalization |
| `pca()` | Principal Component Analysis |
| `plsda()` | PLS Discriminant Analysis |
| `hierarchical_clustering()` | Hierarchical sample clustering |
| `anova()` | One-way ANOVA with BH-FDR |
| `volcano_data()` | Volcano plot data computation |
| `target_decoy_fdr()` | FDR control via target-decoy |
| `calculate_peptide_mass()` | Monoisotopic peptide mass |
| `generate_theoretical_fragments()` | Theoretical b/y ions |

### I/O

| Function | Description |
|----------|-------------|
| `load_mzml()` / `load_mzxml()` | Load standard MS files |
| `load_mgf()` / `save_mgf()` | MGF spectral library format |
| `save_mztab()` | Export to mzTab |
| `save_mzidentml()` | Export identification results |
| `load_feature_table()` / `save_feature_table()` | CSV/TSV feature tables |
| `IndexedMzMLReader` / `IndexedMzXMLReader` | Random-access indexed reading |

---

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

---

## Testing

```bash
# Python (10 test files covering all modules)
cd python && pip install -e ".[dev]" && pytest tests/ -v

# C++
cd core/build && cmake .. -DBUILD_TESTS=ON && make && ctest

# Java
cd java && mvn test
```

---

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
