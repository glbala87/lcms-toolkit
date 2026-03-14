#pragma once

#include "../spectrum.hpp"
#include "../types.hpp"
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>

namespace lcms {
namespace algorithms {

/**
 * @brief Result of a spectral library search.
 */
struct SpectralMatch {
    std::string name;
    double score = 0.0;
    int matched_peaks = 0;
    MZ precursor_mz = 0.0;
    MetaData metadata;
};

/**
 * @brief A spectrum entry in a spectral library.
 */
struct LibrarySpectrum {
    std::string name;
    MZ precursor_mz = 0.0;
    std::vector<MZ> mz;
    std::vector<Intensity> intensity;
    MetaData metadata;

    [[nodiscard]] size_t numPeaks() const { return mz.size(); }
};

/**
 * @brief Spectral library for MS/MS matching.
 */
class SpectralLibrary {
public:
    SpectralLibrary() = default;

    /// Add a spectrum entry
    void addSpectrum(const LibrarySpectrum& entry) {
        entries_.push_back(entry);
    }

    /// Number of entries
    [[nodiscard]] size_t size() const { return entries_.size(); }

    /// Access entry by index
    [[nodiscard]] const LibrarySpectrum& operator[](size_t idx) const {
        return entries_[idx];
    }

    /// Load from MGF file
    size_t loadMGF(const std::string& filepath);

    /// Load from MSP file
    size_t loadMSP(const std::string& filepath);

    /// Save to MGF file
    void saveMGF(const std::string& filepath) const;

    /// Search library
    std::vector<SpectralMatch> search(
        const std::vector<MZ>& query_mz,
        const std::vector<Intensity>& query_intensity,
        MZ query_precursor = 0.0,
        const std::string& method = "cosine",
        double tolerance = 0.01,
        double min_score = 0.0,
        int top_n = 10,
        double precursor_tolerance = 0.0) const;

private:
    std::vector<LibrarySpectrum> entries_;
};

/**
 * @brief Compute cosine similarity between two spectra.
 * @return Pair of (score, matched_peaks_count)
 */
std::pair<double, int> cosineSimilarity(
    const std::vector<MZ>& mz1, const std::vector<Intensity>& int1,
    const std::vector<MZ>& mz2, const std::vector<Intensity>& int2,
    double tolerance = 0.01);

/**
 * @brief Compute modified cosine similarity (shift-aware).
 */
std::pair<double, int> modifiedCosineSimilarity(
    const std::vector<MZ>& mz1, const std::vector<Intensity>& int1, MZ precursor1,
    const std::vector<MZ>& mz2, const std::vector<Intensity>& int2, MZ precursor2,
    double tolerance = 0.01);

/**
 * @brief Compute spectral entropy similarity.
 */
std::pair<double, int> spectralEntropySimilarity(
    const std::vector<MZ>& mz1, const std::vector<Intensity>& int1,
    const std::vector<MZ>& mz2, const std::vector<Intensity>& int2,
    double tolerance = 0.01);

} // namespace algorithms
} // namespace lcms
