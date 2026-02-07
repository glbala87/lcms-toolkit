#pragma once

#include "../ms_experiment.hpp"
#include <string>
#include <memory>
#include <functional>
#include <stdexcept>

namespace lcms {
namespace io {

/**
 * @brief Exception thrown when mzML parsing fails.
 */
class MzMLParseError : public std::runtime_error {
public:
    explicit MzMLParseError(const std::string& msg)
        : std::runtime_error("mzML parse error: " + msg) {}
};

/**
 * @brief Progress callback signature.
 *
 * @param current Current progress (e.g., spectrum count)
 * @param total Total expected items (-1 if unknown)
 * @return false to cancel loading, true to continue
 */
using ProgressCallback = std::function<bool(int current, int total)>;

/**
 * @brief Options for mzML reading.
 */
struct MzMLReaderOptions {
    /// Only load spectra at these MS levels (empty = all)
    std::vector<MSLevel> ms_levels;

    /// RT range filter (only load spectra in this range)
    std::optional<RTRange> rt_range;

    /// M/Z range filter (only keep data in this range)
    std::optional<MZRange> mz_range;

    /// Skip chromatogram loading
    bool skip_chromatograms = false;

    /// Skip spectrum loading (chromatograms only)
    bool skip_spectra = false;

    /// Maximum number of spectra to load (0 = unlimited)
    std::size_t max_spectra = 0;

    /// Load intensity data (false = load m/z only)
    bool load_intensity = true;

    /// Sort spectra by m/z after loading
    bool sort_by_mz = true;

    /// Intensity threshold (discard points below)
    Intensity intensity_threshold = 0.0;

    /// Progress callback (called periodically during loading)
    ProgressCallback progress_callback = nullptr;
};

/**
 * @brief Reader for mzML (Proteomics Standards Initiative) format files.
 *
 * The mzML format is an XML-based format for mass spectrometry data.
 * It supports both profile and centroid data, MS and MS/MS spectra,
 * and various compression schemes.
 *
 * This reader supports:
 * - mzML 1.0 and 1.1 formats
 * - Base64 encoded binary data
 * - zlib compression
 * - Both 32-bit and 64-bit precision
 *
 * Usage:
 * @code
 * MzMLReader reader;
 * MSExperiment exp = reader.read("sample.mzML");
 * @endcode
 */
class MzMLReader {
public:
    MzMLReader();
    ~MzMLReader();

    // Non-copyable
    MzMLReader(const MzMLReader&) = delete;
    MzMLReader& operator=(const MzMLReader&) = delete;

    // Movable
    MzMLReader(MzMLReader&&) noexcept;
    MzMLReader& operator=(MzMLReader&&) noexcept;

    /**
     * @brief Read an mzML file.
     *
     * @param filename Path to the mzML file
     * @return Parsed MSExperiment
     * @throws MzMLParseError if parsing fails
     */
    MSExperiment read(const std::string& filename);

    /**
     * @brief Read an mzML file with options.
     *
     * @param filename Path to the mzML file
     * @param options Reading options
     * @return Parsed MSExperiment
     * @throws MzMLParseError if parsing fails
     */
    MSExperiment read(const std::string& filename, const MzMLReaderOptions& options);

    /**
     * @brief Parse mzML content from a string.
     *
     * @param content mzML XML content
     * @return Parsed MSExperiment
     * @throws MzMLParseError if parsing fails
     */
    MSExperiment parseString(const std::string& content);

    /**
     * @brief Parse mzML content with options.
     *
     * @param content mzML XML content
     * @param options Reading options
     * @return Parsed MSExperiment
     * @throws MzMLParseError if parsing fails
     */
    MSExperiment parseString(const std::string& content,
                             const MzMLReaderOptions& options);

    /**
     * @brief Get the number of spectra in a file without loading data.
     *
     * @param filename Path to the mzML file
     * @return Number of spectra
     */
    std::size_t countSpectra(const std::string& filename);

    /**
     * @brief Get the number of chromatograms in a file.
     *
     * @param filename Path to the mzML file
     * @return Number of chromatograms
     */
    std::size_t countChromatograms(const std::string& filename);

    /**
     * @brief Check if a file appears to be valid mzML.
     *
     * @param filename Path to check
     * @return true if file appears to be mzML
     */
    static bool isValidMzML(const std::string& filename);

    /**
     * @brief Set default options for subsequent reads.
     *
     * @param options Default options
     */
    void setDefaultOptions(const MzMLReaderOptions& options) {
        default_options_ = options;
    }

    /**
     * @brief Get the last error message (if any).
     */
    [[nodiscard]] const std::string& lastError() const noexcept {
        return last_error_;
    }

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
    MzMLReaderOptions default_options_;
    std::string last_error_;
};

/**
 * @brief Convenience function to load an mzML file.
 *
 * @param filename Path to the mzML file
 * @return Parsed MSExperiment
 */
inline MSExperiment loadMzML(const std::string& filename) {
    MzMLReader reader;
    return reader.read(filename);
}

/**
 * @brief Convenience function to load an mzML file with options.
 *
 * @param filename Path to the mzML file
 * @param options Reading options
 * @return Parsed MSExperiment
 */
inline MSExperiment loadMzML(const std::string& filename,
                             const MzMLReaderOptions& options) {
    MzMLReader reader;
    return reader.read(filename, options);
}

} // namespace io
} // namespace lcms
