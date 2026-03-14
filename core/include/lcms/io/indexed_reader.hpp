#pragma once

#include "../spectrum.hpp"
#include "../types.hpp"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <optional>
#include <memory>

namespace lcms {
namespace io {

/**
 * @brief Index entry for a spectrum in a file.
 */
struct SpectrumIndexEntry {
    Index scan_number = 0;
    std::streamoff offset = 0;
    std::streamoff length = 0;
    MSLevel ms_level = 1;
    RetentionTime rt = 0.0;
    MZ precursor_mz = 0.0;
};

/**
 * @brief File index for random access to spectra.
 */
class FileIndex {
public:
    FileIndex() = default;

    /// Number of indexed spectra
    [[nodiscard]] size_t size() const { return entries_.size(); }

    /// Add an entry
    void addEntry(const SpectrumIndexEntry& entry);

    /// Get entry by sequential index
    [[nodiscard]] const SpectrumIndexEntry* getByIndex(size_t index) const;

    /// Get entry by scan number
    [[nodiscard]] const SpectrumIndexEntry* getByScan(Index scan) const;

    /// Get entry closest to given RT
    [[nodiscard]] const SpectrumIndexEntry* getByRT(RetentionTime rt,
                                                      double tolerance = 0.01) const;

    /// Get RT range
    [[nodiscard]] RTRange rtRange() const;

    /// Get unique MS levels
    [[nodiscard]] std::vector<MSLevel> msLevels() const;

    /// Save index to file
    void save(const std::string& filepath) const;

    /// Load index from file
    static FileIndex load(const std::string& filepath);

    /// Build index from mzML file
    static FileIndex buildFromMzML(const std::string& filepath);

    /// Build index from mzXML file
    static FileIndex buildFromMzXML(const std::string& filepath);

private:
    std::vector<SpectrumIndexEntry> entries_;
    std::map<Index, size_t> scan_map_;
};

/**
 * @brief Indexed reader for random access to mzML files.
 */
class IndexedMzMLReader {
public:
    explicit IndexedMzMLReader(const std::string& filepath);
    ~IndexedMzMLReader();

    /// Open file and build/load index
    bool open();

    /// Close file
    void close();

    /// Number of spectra
    [[nodiscard]] size_t spectrumCount() const;

    /// Get spectrum by index
    [[nodiscard]] std::optional<Spectrum> getSpectrum(size_t index);

    /// Get spectrum by scan number
    [[nodiscard]] std::optional<Spectrum> getSpectrumByScan(Index scan);

    /// Get spectrum by RT
    [[nodiscard]] std::optional<Spectrum> getSpectrumByRT(RetentionTime rt,
                                                           double tolerance = 0.01);

    /// Save index for future fast loading
    void saveIndex(const std::string& filepath = "") const;

    /// Get the file index
    [[nodiscard]] const FileIndex& index() const { return index_; }

private:
    std::string filepath_;
    FileIndex index_;
    std::unique_ptr<std::ifstream> file_;
    bool is_open_ = false;

    std::optional<Spectrum> readSpectrumAt(const SpectrumIndexEntry& entry);
};

/**
 * @brief Indexed reader for random access to mzXML files.
 */
class IndexedMzXMLReader {
public:
    explicit IndexedMzXMLReader(const std::string& filepath);
    ~IndexedMzXMLReader();

    bool open();
    void close();

    [[nodiscard]] size_t spectrumCount() const;
    [[nodiscard]] std::optional<Spectrum> getSpectrum(size_t index);
    [[nodiscard]] std::optional<Spectrum> getSpectrumByScan(Index scan);
    [[nodiscard]] std::optional<Spectrum> getSpectrumByRT(RetentionTime rt,
                                                           double tolerance = 0.01);
    void saveIndex(const std::string& filepath = "") const;
    [[nodiscard]] const FileIndex& index() const { return index_; }

private:
    std::string filepath_;
    FileIndex index_;
    std::unique_ptr<std::ifstream> file_;
    bool is_open_ = false;

    std::optional<Spectrum> readSpectrumAt(const SpectrumIndexEntry& entry);
};

} // namespace io
} // namespace lcms
