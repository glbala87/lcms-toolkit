#pragma once

#include "types.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>

namespace lcms {

/**
 * @brief Type of chromatogram
 */
enum class ChromatogramType : std::uint8_t {
    UNKNOWN = 0,
    TIC,            // Total Ion Current
    BPC,            // Base Peak Chromatogram
    XIC,            // Extracted Ion Chromatogram
    SRM,            // Selected Reaction Monitoring
    MRM,            // Multiple Reaction Monitoring
    SIM,            // Selected Ion Monitoring
    ABSORPTION,     // UV/Vis absorption
    EMISSION        // Fluorescence emission
};

/**
 * @brief Represents a chromatogram (intensity vs retention time).
 *
 * A Chromatogram stores paired arrays of retention times and intensities,
 * along with metadata describing the chromatogram type and acquisition
 * parameters.
 */
class Chromatogram {
public:
    /// Default constructor creates an empty chromatogram
    Chromatogram() = default;

    /// Construct from retention time and intensity vectors
    Chromatogram(std::vector<RetentionTime> rt, std::vector<Intensity> intensity)
        : rt_(std::move(rt)), intensity_(std::move(intensity)) {
        if (rt_.size() != intensity_.size()) {
            throw std::invalid_argument(
                "RT and intensity arrays must have same size");
        }
        updateRanges();
    }

    /// Move constructor
    Chromatogram(Chromatogram&&) noexcept = default;

    /// Move assignment
    Chromatogram& operator=(Chromatogram&&) noexcept = default;

    /// Copy constructor
    Chromatogram(const Chromatogram&) = default;

    /// Copy assignment
    Chromatogram& operator=(const Chromatogram&) = default;

    // =========================================================================
    // Data Access
    // =========================================================================

    /// Get number of data points
    [[nodiscard]] std::size_t size() const noexcept { return rt_.size(); }

    /// Check if chromatogram is empty
    [[nodiscard]] bool empty() const noexcept { return rt_.empty(); }

    /// Get retention time array (const reference)
    [[nodiscard]] const std::vector<RetentionTime>& rt() const noexcept {
        return rt_;
    }

    /// Get intensity array (const reference)
    [[nodiscard]] const std::vector<Intensity>& intensity() const noexcept {
        return intensity_;
    }

    /// Get retention time at index
    [[nodiscard]] RetentionTime rtAt(Index i) const { return rt_.at(i); }

    /// Get intensity at index
    [[nodiscard]] Intensity intensityAt(Index i) const {
        return intensity_.at(i);
    }

    /// Get RT array (mutable)
    std::vector<RetentionTime>& rt() noexcept { return rt_; }

    /// Get intensity array (mutable)
    std::vector<Intensity>& intensity() noexcept { return intensity_; }

    // =========================================================================
    // Metadata
    // =========================================================================

    /// Get chromatogram index
    [[nodiscard]] Index index() const noexcept { return index_; }
    void setIndex(Index idx) noexcept { index_ = idx; }

    /// Get native ID
    [[nodiscard]] const std::string& nativeId() const noexcept {
        return native_id_;
    }
    void setNativeId(std::string id) { native_id_ = std::move(id); }

    /// Get chromatogram type
    [[nodiscard]] ChromatogramType type() const noexcept { return type_; }
    void setType(ChromatogramType t) noexcept { type_ = t; }

    /// Get polarity
    [[nodiscard]] Polarity polarity() const noexcept { return polarity_; }
    void setPolarity(Polarity p) noexcept { polarity_ = p; }

    /// Get target m/z (for XIC, SIM, SRM)
    [[nodiscard]] MZ targetMz() const noexcept { return target_mz_; }
    void setTargetMz(MZ mz) noexcept { target_mz_ = mz; }

    /// Get m/z tolerance window
    [[nodiscard]] MZ mzTolerance() const noexcept { return mz_tolerance_; }
    void setMzTolerance(MZ tol) noexcept { mz_tolerance_ = tol; }

    /// Get precursor m/z (for SRM/MRM)
    [[nodiscard]] MZ precursorMz() const noexcept { return precursor_mz_; }
    void setPrecursorMz(MZ mz) noexcept { precursor_mz_ = mz; }

    /// Get product m/z (for SRM/MRM)
    [[nodiscard]] MZ productMz() const noexcept { return product_mz_; }
    void setProductMz(MZ mz) noexcept { product_mz_ = mz; }

    /// Get RT range
    [[nodiscard]] const RTRange& rtRange() const noexcept { return rt_range_; }

    /// Get intensity range
    [[nodiscard]] const IntensityRange& intensityRange() const noexcept {
        return intensity_range_;
    }

    /// Get maximum intensity
    [[nodiscard]] Intensity maxIntensity() const noexcept {
        return intensity_range_.max_value;
    }

    /// Get RT at maximum intensity
    [[nodiscard]] RetentionTime rtAtMaxIntensity() const noexcept {
        return rt_at_max_;
    }

    /// Get custom metadata
    [[nodiscard]] const MetaData& metadata() const noexcept { return metadata_; }
    MetaData& metadata() noexcept { return metadata_; }

    // =========================================================================
    // Data Operations
    // =========================================================================

    /// Set data from vectors
    void setData(std::vector<RetentionTime> rt, std::vector<Intensity> intensity) {
        if (rt.size() != intensity.size()) {
            throw std::invalid_argument(
                "RT and intensity arrays must have same size");
        }
        rt_ = std::move(rt);
        intensity_ = std::move(intensity);
        updateRanges();
    }

    /// Reserve capacity for data
    void reserve(std::size_t n) {
        rt_.reserve(n);
        intensity_.reserve(n);
    }

    /// Clear all data
    void clear() {
        rt_.clear();
        intensity_.clear();
        rt_range_ = RTRange();
        intensity_range_ = IntensityRange();
        rt_at_max_ = 0;
    }

    /// Add a single data point
    void addPoint(RetentionTime rt, Intensity intensity) {
        rt_.push_back(rt);
        intensity_.push_back(intensity);
        rt_range_.extend(rt);
        if (intensity > intensity_range_.max_value) {
            intensity_range_.max_value = intensity;
            rt_at_max_ = rt;
        }
        if (intensity < intensity_range_.min_value) {
            intensity_range_.min_value = intensity;
        }
    }

    /// Sort by retention time (ascending)
    void sortByRt() {
        if (rt_.size() <= 1) return;

        std::vector<Index> indices(rt_.size());
        std::iota(indices.begin(), indices.end(), 0);

        std::sort(indices.begin(), indices.end(),
            [this](Index a, Index b) { return rt_[a] < rt_[b]; });

        std::vector<RetentionTime> new_rt(rt_.size());
        std::vector<Intensity> new_intensity(intensity_.size());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            new_rt[i] = rt_[indices[i]];
            new_intensity[i] = intensity_[indices[i]];
        }
        rt_ = std::move(new_rt);
        intensity_ = std::move(new_intensity);
    }

    /// Check if data is sorted by RT
    [[nodiscard]] bool isSortedByRt() const {
        return std::is_sorted(rt_.begin(), rt_.end());
    }

    /// Find index of point closest to given RT
    [[nodiscard]] Index findNearestRt(RetentionTime target) const {
        if (empty()) {
            throw std::runtime_error("Cannot find point in empty chromatogram");
        }

        auto it = std::lower_bound(rt_.begin(), rt_.end(), target);

        if (it == rt_.end()) {
            return rt_.size() - 1;
        }
        if (it == rt_.begin()) {
            return 0;
        }

        auto prev = std::prev(it);
        if (std::abs(*it - target) < std::abs(*prev - target)) {
            return std::distance(rt_.begin(), it);
        }
        return std::distance(rt_.begin(), prev);
    }

    /// Extract points within RT range
    [[nodiscard]] Chromatogram extractRange(RetentionTime low,
                                             RetentionTime high) const {
        Chromatogram result;
        result.setType(type_);
        result.setPolarity(polarity_);
        result.setTargetMz(target_mz_);

        for (std::size_t i = 0; i < rt_.size(); ++i) {
            if (rt_[i] >= low && rt_[i] <= high) {
                result.addPoint(rt_[i], intensity_[i]);
            }
        }
        return result;
    }

    /// Compute area under the curve (trapezoidal integration)
    [[nodiscard]] double computeArea() const {
        if (size() < 2) return 0.0;

        double area = 0.0;
        for (std::size_t i = 1; i < rt_.size(); ++i) {
            double dt = rt_[i] - rt_[i-1];
            area += 0.5 * (intensity_[i] + intensity_[i-1]) * dt;
        }
        return area;
    }

    /// Compute area within RT range
    [[nodiscard]] double computeArea(RetentionTime low, RetentionTime high) const {
        return extractRange(low, high).computeArea();
    }

    /// Interpolate intensity at given RT
    [[nodiscard]] Intensity interpolateAt(RetentionTime target) const {
        if (empty()) return 0.0;
        if (rt_.size() == 1) return intensity_[0];

        // Find bracketing points
        auto it = std::lower_bound(rt_.begin(), rt_.end(), target);

        if (it == rt_.end()) {
            return intensity_.back();
        }
        if (it == rt_.begin()) {
            return intensity_.front();
        }

        std::size_t i = std::distance(rt_.begin(), it);

        // Linear interpolation
        double t = (target - rt_[i-1]) / (rt_[i] - rt_[i-1]);
        return intensity_[i-1] + t * (intensity_[i] - intensity_[i-1]);
    }

    /// Update cached statistics
    void updateRanges() {
        rt_range_ = RTRange();
        intensity_range_ = IntensityRange();
        rt_at_max_ = 0;

        for (std::size_t i = 0; i < rt_.size(); ++i) {
            rt_range_.extend(rt_[i]);
            intensity_range_.extend(intensity_[i]);
            if (intensity_[i] >= intensity_range_.max_value) {
                intensity_range_.max_value = intensity_[i];
                rt_at_max_ = rt_[i];
            }
        }
    }

private:
    // Raw data
    std::vector<RetentionTime> rt_;
    std::vector<Intensity> intensity_;

    // Identification
    Index index_ = 0;
    std::string native_id_;

    // Type and acquisition
    ChromatogramType type_ = ChromatogramType::UNKNOWN;
    Polarity polarity_ = Polarity::UNKNOWN;

    // Target parameters
    MZ target_mz_ = 0.0;
    MZ mz_tolerance_ = 0.0;
    MZ precursor_mz_ = 0.0;
    MZ product_mz_ = 0.0;

    // Cached statistics
    RTRange rt_range_;
    IntensityRange intensity_range_;
    RetentionTime rt_at_max_ = 0;

    // Custom metadata
    MetaData metadata_;
};

/// Convert chromatogram type to string
inline std::string toString(ChromatogramType t) {
    switch (t) {
        case ChromatogramType::TIC: return "TIC";
        case ChromatogramType::BPC: return "BPC";
        case ChromatogramType::XIC: return "XIC";
        case ChromatogramType::SRM: return "SRM";
        case ChromatogramType::MRM: return "MRM";
        case ChromatogramType::SIM: return "SIM";
        case ChromatogramType::ABSORPTION: return "absorption";
        case ChromatogramType::EMISSION: return "emission";
        default: return "unknown";
    }
}

} // namespace lcms
