#pragma once

#include "../ms_experiment.hpp"
#include <string>
#include <memory>
#include <functional>
#include <stdexcept>

namespace lcms {
namespace io {

/**
 * @brief Exception thrown when mzXML parsing fails.
 */
class MzXMLParseError : public std::runtime_error {
public:
    explicit MzXMLParseError(const std::string& msg)
        : std::runtime_error("mzXML parse error: " + msg) {}
};

/**
 * @brief Progress callback for mzXML loading.
 */
using MzXMLProgressCallback = std::function<bool(int current, int total)>;

/**
 * @brief Options for mzXML reading.
 */
struct MzXMLReaderOptions {
    /// Only load spectra at these MS levels (empty = all)
    std::vector<MSLevel> ms_levels;

    /// RT range filter
    std::optional<RTRange> rt_range;

    /// M/Z range filter
    std::optional<MZRange> mz_range;

    /// Maximum number of spectra to load (0 = unlimited)
    std::size_t max_spectra = 0;

    /// Load intensity data
    bool load_intensity = true;

    /// Sort spectra by m/z after loading
    bool sort_by_mz = true;

    /// Intensity threshold
    Intensity intensity_threshold = 0.0;

    /// Progress callback
    MzXMLProgressCallback progress_callback = nullptr;
};

/**
 * @brief Reader for mzXML format files.
 *
 * mzXML is a legacy XML-based format for mass spectrometry data,
 * developed at the Institute for Systems Biology. While superseded
 * by mzML, many existing datasets still use this format.
 *
 * This reader supports:
 * - mzXML 1.x, 2.x, and 3.x formats
 * - Base64 encoded binary data
 * - zlib compression
 * - Network byte order (big endian)
 *
 * Usage:
 * @code
 * MzXMLReader reader;
 * MSExperiment exp = reader.read("sample.mzXML");
 * @endcode
 */
class MzXMLReader {
public:
    MzXMLReader();
    ~MzXMLReader();

    // Non-copyable
    MzXMLReader(const MzXMLReader&) = delete;
    MzXMLReader& operator=(const MzXMLReader&) = delete;

    // Movable
    MzXMLReader(MzXMLReader&&) noexcept;
    MzXMLReader& operator=(MzXMLReader&&) noexcept;

    /**
     * @brief Read an mzXML file.
     *
     * @param filename Path to the mzXML file
     * @return Parsed MSExperiment
     * @throws MzXMLParseError if parsing fails
     */
    MSExperiment read(const std::string& filename);

    /**
     * @brief Read an mzXML file with options.
     *
     * @param filename Path to the mzXML file
     * @param options Reading options
     * @return Parsed MSExperiment
     * @throws MzXMLParseError if parsing fails
     */
    MSExperiment read(const std::string& filename, const MzXMLReaderOptions& options);

    /**
     * @brief Parse mzXML content from a string.
     *
     * @param content mzXML XML content
     * @return Parsed MSExperiment
     */
    MSExperiment parseString(const std::string& content);

    /**
     * @brief Parse mzXML content with options.
     *
     * @param content mzXML XML content
     * @param options Reading options
     * @return Parsed MSExperiment
     */
    MSExperiment parseString(const std::string& content,
                             const MzXMLReaderOptions& options);

    /**
     * @brief Get the number of spectra in a file.
     *
     * @param filename Path to the mzXML file
     * @return Number of spectra
     */
    std::size_t countSpectra(const std::string& filename);

    /**
     * @brief Check if a file appears to be valid mzXML.
     *
     * @param filename Path to check
     * @return true if file appears to be mzXML
     */
    static bool isValidMzXML(const std::string& filename);

    /**
     * @brief Set default options.
     */
    void setDefaultOptions(const MzXMLReaderOptions& options) {
        default_options_ = options;
    }

    /**
     * @brief Get the last error message.
     */
    [[nodiscard]] const std::string& lastError() const noexcept {
        return last_error_;
    }

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
    MzXMLReaderOptions default_options_;
    std::string last_error_;
};

/**
 * @brief Convenience function to load an mzXML file.
 */
inline MSExperiment loadMzXML(const std::string& filename) {
    MzXMLReader reader;
    return reader.read(filename);
}

/**
 * @brief Convenience function to load an mzXML file with options.
 */
inline MSExperiment loadMzXML(const std::string& filename,
                              const MzXMLReaderOptions& options) {
    MzXMLReader reader;
    return reader.read(filename, options);
}

} // namespace io
} // namespace lcms
