#pragma once

#include "types.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace lcms {

/**
 * @brief Represents a mass spectrum with m/z and intensity data.
 *
 * A Spectrum stores paired arrays of m/z values and corresponding intensities,
 * along with metadata about the spectrum acquisition. It provides efficient
 * access to the raw data as well as common operations like finding peaks,
 * computing statistics, and extracting m/z ranges.
 */
class Spectrum {
public:
    /// Default constructor creates an empty spectrum
    Spectrum() = default;

    /// Construct from m/z and intensity vectors
    Spectrum(std::vector<MZ> mz, std::vector<Intensity> intensity)
        : mz_(std::move(mz)), intensity_(std::move(intensity)) {
        if (mz_.size() != intensity_.size()) {
            throw std::invalid_argument("m/z and intensity arrays must have same size");
        }
        updateRanges();
    }

    /// Move constructor
    Spectrum(Spectrum&&) noexcept = default;

    /// Move assignment
    Spectrum& operator=(Spectrum&&) noexcept = default;

    /// Copy constructor
    Spectrum(const Spectrum&) = default;

    /// Copy assignment
    Spectrum& operator=(const Spectrum&) = default;

    // =========================================================================
    // Data Access
    // =========================================================================

    /// Get number of data points
    [[nodiscard]] std::size_t size() const noexcept { return mz_.size(); }

    /// Check if spectrum is empty
    [[nodiscard]] bool empty() const noexcept { return mz_.empty(); }

    /// Get m/z array (const reference)
    [[nodiscard]] const std::vector<MZ>& mz() const noexcept { return mz_; }

    /// Get intensity array (const reference)
    [[nodiscard]] const std::vector<Intensity>& intensity() const noexcept {
        return intensity_;
    }

    /// Get m/z at index
    [[nodiscard]] MZ mzAt(Index i) const { return mz_.at(i); }

    /// Get intensity at index
    [[nodiscard]] Intensity intensityAt(Index i) const { return intensity_.at(i); }

    /// Get m/z array (mutable)
    std::vector<MZ>& mz() noexcept { return mz_; }

    /// Get intensity array (mutable)
    std::vector<Intensity>& intensity() noexcept { return intensity_; }

    // =========================================================================
    // Metadata
    // =========================================================================

    /// Get spectrum index in the run
    [[nodiscard]] Index index() const noexcept { return index_; }
    void setIndex(Index idx) noexcept { index_ = idx; }

    /// Get native ID (from source file)
    [[nodiscard]] const std::string& nativeId() const noexcept { return native_id_; }
    void setNativeId(std::string id) { native_id_ = std::move(id); }

    /// Get MS level (1 for MS1, 2 for MS/MS, etc.)
    [[nodiscard]] MSLevel msLevel() const noexcept { return ms_level_; }
    void setMsLevel(MSLevel level) noexcept { ms_level_ = level; }

    /// Get retention time in seconds
    [[nodiscard]] RetentionTime retentionTime() const noexcept { return rt_; }
    void setRetentionTime(RetentionTime rt) noexcept { rt_ = rt; }

    /// Get spectrum type (profile or centroid)
    [[nodiscard]] SpectrumType type() const noexcept { return type_; }
    void setType(SpectrumType t) noexcept { type_ = t; }

    /// Get polarity
    [[nodiscard]] Polarity polarity() const noexcept { return polarity_; }
    void setPolarity(Polarity p) noexcept { polarity_ = p; }

    /// Get precursor information (for MS/MS spectra)
    [[nodiscard]] const std::vector<Precursor>& precursors() const noexcept {
        return precursors_;
    }
    void addPrecursor(Precursor p) { precursors_.push_back(std::move(p)); }
    void clearPrecursors() { precursors_.clear(); }

    /// Get total ion current
    [[nodiscard]] Intensity tic() const noexcept { return tic_; }
    void setTic(Intensity tic) noexcept { tic_ = tic; }

    /// Get base peak intensity
    [[nodiscard]] Intensity basePeakIntensity() const noexcept {
        return base_peak_intensity_;
    }

    /// Get base peak m/z
    [[nodiscard]] MZ basePeakMz() const noexcept { return base_peak_mz_; }

    /// Get m/z range
    [[nodiscard]] const MZRange& mzRange() const noexcept { return mz_range_; }

    /// Get custom metadata
    [[nodiscard]] const MetaData& metadata() const noexcept { return metadata_; }
    MetaData& metadata() noexcept { return metadata_; }

    // =========================================================================
    // Data Operations
    // =========================================================================

    /// Set data from vectors
    void setData(std::vector<MZ> mz, std::vector<Intensity> intensity) {
        if (mz.size() != intensity.size()) {
            throw std::invalid_argument("m/z and intensity arrays must have same size");
        }
        mz_ = std::move(mz);
        intensity_ = std::move(intensity);
        updateRanges();
    }

    /// Reserve capacity for data
    void reserve(std::size_t n) {
        mz_.reserve(n);
        intensity_.reserve(n);
    }

    /// Clear all data
    void clear() {
        mz_.clear();
        intensity_.clear();
        mz_range_ = MZRange();
        base_peak_intensity_ = 0;
        base_peak_mz_ = 0;
        tic_ = 0;
    }

    /// Add a single data point
    void addPoint(MZ mz, Intensity intensity) {
        mz_.push_back(mz);
        intensity_.push_back(intensity);
        mz_range_.extend(mz);
        if (intensity > base_peak_intensity_) {
            base_peak_intensity_ = intensity;
            base_peak_mz_ = mz;
        }
        tic_ += intensity;
    }

    /// Sort by m/z (ascending)
    void sortByMz() {
        if (mz_.size() <= 1) return;

        // Create index array
        std::vector<Index> indices(mz_.size());
        std::iota(indices.begin(), indices.end(), 0);

        // Sort indices by m/z
        std::sort(indices.begin(), indices.end(),
            [this](Index a, Index b) { return mz_[a] < mz_[b]; });

        // Reorder arrays
        std::vector<MZ> new_mz(mz_.size());
        std::vector<Intensity> new_intensity(intensity_.size());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            new_mz[i] = mz_[indices[i]];
            new_intensity[i] = intensity_[indices[i]];
        }
        mz_ = std::move(new_mz);
        intensity_ = std::move(new_intensity);
    }

    /// Check if data is sorted by m/z
    [[nodiscard]] bool isSortedByMz() const {
        return std::is_sorted(mz_.begin(), mz_.end());
    }

    /// Find index of peak closest to given m/z (requires sorted data)
    [[nodiscard]] Index findNearestMz(MZ target) const {
        if (empty()) {
            throw std::runtime_error("Cannot find peak in empty spectrum");
        }

        auto it = std::lower_bound(mz_.begin(), mz_.end(), target);

        if (it == mz_.end()) {
            return mz_.size() - 1;
        }
        if (it == mz_.begin()) {
            return 0;
        }

        // Check which neighbor is closer
        auto prev = std::prev(it);
        if (std::abs(*it - target) < std::abs(*prev - target)) {
            return std::distance(mz_.begin(), it);
        }
        return std::distance(mz_.begin(), prev);
    }

    /// Extract peaks within m/z range
    [[nodiscard]] Spectrum extractRange(MZ low, MZ high) const {
        Spectrum result;
        result.setMsLevel(ms_level_);
        result.setRetentionTime(rt_);
        result.setType(type_);
        result.setPolarity(polarity_);

        for (std::size_t i = 0; i < mz_.size(); ++i) {
            if (mz_[i] >= low && mz_[i] <= high) {
                result.addPoint(mz_[i], intensity_[i]);
            }
        }
        return result;
    }

    /// Compute sum of intensities
    [[nodiscard]] Intensity sumIntensity() const {
        return std::accumulate(intensity_.begin(), intensity_.end(), 0.0);
    }

    /// Update cached statistics (TIC, base peak, etc.)
    void updateRanges() {
        mz_range_ = MZRange();
        base_peak_intensity_ = 0;
        base_peak_mz_ = 0;
        tic_ = 0;

        for (std::size_t i = 0; i < mz_.size(); ++i) {
            mz_range_.extend(mz_[i]);
            tic_ += intensity_[i];
            if (intensity_[i] > base_peak_intensity_) {
                base_peak_intensity_ = intensity_[i];
                base_peak_mz_ = mz_[i];
            }
        }
    }

private:
    // Raw data
    std::vector<MZ> mz_;
    std::vector<Intensity> intensity_;

    // Spectrum identification
    Index index_ = 0;
    std::string native_id_;

    // Acquisition parameters
    MSLevel ms_level_ = 1;
    RetentionTime rt_ = 0.0;
    SpectrumType type_ = SpectrumType::UNKNOWN;
    Polarity polarity_ = Polarity::UNKNOWN;

    // Precursor information (for MS/MS)
    std::vector<Precursor> precursors_;

    // Cached statistics
    MZRange mz_range_;
    Intensity tic_ = 0;
    Intensity base_peak_intensity_ = 0;
    MZ base_peak_mz_ = 0;

    // Custom metadata
    MetaData metadata_;
};

} // namespace lcms
