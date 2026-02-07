#pragma once

#include "types.hpp"
#include "peak.hpp"
#include <vector>
#include <memory>
#include <algorithm>

namespace lcms {

/**
 * @brief Represents a feature (2D peak) in LC-MS data.
 *
 * A Feature represents a molecular entity detected across multiple spectra,
 * with a defined m/z and RT extent, intensity profile, and associated peaks.
 */
class Feature {
public:
    Feature() = default;

    /// Construct with basic parameters
    Feature(MZ mz, RetentionTime rt, Intensity intensity)
        : mz_(mz), rt_(rt), intensity_(intensity) {}

    // =========================================================================
    // Position and Intensity
    // =========================================================================

    /// Get m/z (centroid or monoisotopic)
    [[nodiscard]] MZ mz() const noexcept { return mz_; }
    void setMz(MZ mz) noexcept { mz_ = mz; }

    /// Get retention time (apex)
    [[nodiscard]] RetentionTime rt() const noexcept { return rt_; }
    void setRt(RetentionTime rt) noexcept { rt_ = rt; }

    /// Get apex intensity
    [[nodiscard]] Intensity intensity() const noexcept { return intensity_; }
    void setIntensity(Intensity i) noexcept { intensity_ = i; }

    /// Get integrated intensity (volume)
    [[nodiscard]] double volume() const noexcept { return volume_; }
    void setVolume(double v) noexcept { volume_ = v; }

    // =========================================================================
    // Boundaries
    // =========================================================================

    /// Get m/z range
    [[nodiscard]] const MZRange& mzRange() const noexcept { return mz_range_; }
    void setMzRange(MZRange r) noexcept { mz_range_ = r; }
    void setMzRange(MZ low, MZ high) noexcept { mz_range_ = MZRange(low, high); }

    /// Get RT range
    [[nodiscard]] const RTRange& rtRange() const noexcept { return rt_range_; }
    void setRtRange(RTRange r) noexcept { rt_range_ = r; }
    void setRtRange(RetentionTime low, RetentionTime high) noexcept {
        rt_range_ = RTRange(low, high);
    }

    /// Get m/z width
    [[nodiscard]] MZ mzWidth() const noexcept { return mz_range_.span(); }

    /// Get RT width
    [[nodiscard]] RetentionTime rtWidth() const noexcept {
        return rt_range_.span();
    }

    // =========================================================================
    // Charge and Isotopes
    // =========================================================================

    /// Get charge state
    [[nodiscard]] ChargeState charge() const noexcept { return charge_; }
    void setCharge(ChargeState z) noexcept { charge_ = z; }

    /// Check if charge state is known
    [[nodiscard]] bool hasCharge() const noexcept { return charge_ != 0; }

    /// Get monoisotopic m/z
    [[nodiscard]] MZ monoisotopicMz() const noexcept { return mono_mz_; }
    void setMonoisotopicMz(MZ mz) noexcept { mono_mz_ = mz; }

    /// Get number of isotope peaks
    [[nodiscard]] int isotopeCount() const noexcept { return isotope_count_; }
    void setIsotopeCount(int n) noexcept { isotope_count_ = n; }

    /// Compute neutral mass
    [[nodiscard]] double neutralMass() const {
        if (charge_ == 0) return 0.0;
        constexpr double proton_mass = 1.007276;
        MZ mz = (mono_mz_ > 0) ? mono_mz_ : mz_;
        return (mz - proton_mass) * std::abs(charge_);
    }

    // =========================================================================
    // Quality Metrics
    // =========================================================================

    /// Get overall quality score (0-1)
    [[nodiscard]] double quality() const noexcept { return quality_; }
    void setQuality(double q) noexcept { quality_ = q; }

    /// Get signal-to-noise ratio
    [[nodiscard]] double snr() const noexcept { return snr_; }
    void setSnr(double s) noexcept { snr_ = s; }

    /// Get isotope pattern fit score
    [[nodiscard]] double isotopeScore() const noexcept { return isotope_score_; }
    void setIsotopeScore(double s) noexcept { isotope_score_ = s; }

    // =========================================================================
    // Constituent Peaks
    // =========================================================================

    /// Get peaks contributing to this feature
    [[nodiscard]] const PeakList& peaks() const noexcept { return peaks_; }
    PeakList& peaks() noexcept { return peaks_; }

    /// Add a constituent peak
    void addPeak(Peak p) { peaks_.add(std::move(p)); }

    /// Get number of constituent peaks
    [[nodiscard]] std::size_t peakCount() const noexcept { return peaks_.size(); }

    // =========================================================================
    // Identification
    // =========================================================================

    /// Get feature ID
    [[nodiscard]] Index id() const noexcept { return id_; }
    void setId(Index id) noexcept { id_ = id; }

    /// Get custom metadata
    [[nodiscard]] const MetaData& metadata() const noexcept { return metadata_; }
    MetaData& metadata() noexcept { return metadata_; }

    // =========================================================================
    // Utility
    // =========================================================================

    /// Check if point is within feature boundaries
    [[nodiscard]] bool contains(MZ mz, RetentionTime rt) const noexcept {
        return mz_range_.contains(mz) && rt_range_.contains(rt);
    }

    /// Check if another feature overlaps
    [[nodiscard]] bool overlaps(const Feature& other) const noexcept {
        return mz_range_.overlaps(other.mz_range_) &&
               rt_range_.overlaps(other.rt_range_);
    }

    /// Compare by m/z
    [[nodiscard]] bool operator<(const Feature& other) const noexcept {
        return mz_ < other.mz_;
    }

private:
    Index id_ = 0;

    // Position
    MZ mz_ = 0.0;
    RetentionTime rt_ = 0.0;
    Intensity intensity_ = 0.0;
    double volume_ = 0.0;

    // Boundaries
    MZRange mz_range_;
    RTRange rt_range_;

    // Charge/isotope
    ChargeState charge_ = 0;
    MZ mono_mz_ = 0.0;
    int isotope_count_ = 0;

    // Quality
    double quality_ = 0.0;
    double snr_ = 0.0;
    double isotope_score_ = 0.0;

    // Constituent peaks
    PeakList peaks_;

    // Metadata
    MetaData metadata_;
};

/**
 * @brief Container for a collection of features with spatial indexing.
 */
class FeatureMap {
public:
    using iterator = std::vector<Feature>::iterator;
    using const_iterator = std::vector<Feature>::const_iterator;

    FeatureMap() = default;

    // =========================================================================
    // Container Operations
    // =========================================================================

    /// Get number of features
    [[nodiscard]] std::size_t size() const noexcept { return features_.size(); }

    /// Check if empty
    [[nodiscard]] bool empty() const noexcept { return features_.empty(); }

    /// Access feature by index
    [[nodiscard]] Feature& operator[](std::size_t i) { return features_[i]; }
    [[nodiscard]] const Feature& operator[](std::size_t i) const {
        return features_[i];
    }

    /// Iterator access
    iterator begin() noexcept { return features_.begin(); }
    iterator end() noexcept { return features_.end(); }
    const_iterator begin() const noexcept { return features_.begin(); }
    const_iterator end() const noexcept { return features_.end(); }
    const_iterator cbegin() const noexcept { return features_.cbegin(); }
    const_iterator cend() const noexcept { return features_.cend(); }

    /// Add a feature
    void add(Feature f) {
        f.setId(next_id_++);
        updateRanges(f);
        features_.push_back(std::move(f));
    }

    /// Reserve capacity
    void reserve(std::size_t n) { features_.reserve(n); }

    /// Clear all features
    void clear() {
        features_.clear();
        mz_range_ = MZRange();
        rt_range_ = RTRange();
        next_id_ = 0;
    }

    // =========================================================================
    // Sorting
    // =========================================================================

    /// Sort by m/z
    void sortByMz() {
        std::sort(features_.begin(), features_.end(),
            [](const Feature& a, const Feature& b) { return a.mz() < b.mz(); });
    }

    /// Sort by RT
    void sortByRt() {
        std::sort(features_.begin(), features_.end(),
            [](const Feature& a, const Feature& b) { return a.rt() < b.rt(); });
    }

    /// Sort by intensity (descending)
    void sortByIntensity() {
        std::sort(features_.begin(), features_.end(),
            [](const Feature& a, const Feature& b) {
                return a.intensity() > b.intensity();
            });
    }

    /// Sort by quality (descending)
    void sortByQuality() {
        std::sort(features_.begin(), features_.end(),
            [](const Feature& a, const Feature& b) {
                return a.quality() > b.quality();
            });
    }

    // =========================================================================
    // Search
    // =========================================================================

    /// Find features within m/z range
    [[nodiscard]] std::vector<const Feature*> findInMzRange(MZ low,
                                                             MZ high) const {
        std::vector<const Feature*> result;
        for (const auto& f : features_) {
            if (f.mz() >= low && f.mz() <= high) {
                result.push_back(&f);
            }
        }
        return result;
    }

    /// Find features within RT range
    [[nodiscard]] std::vector<const Feature*> findInRtRange(
        RetentionTime low, RetentionTime high) const {
        std::vector<const Feature*> result;
        for (const auto& f : features_) {
            if (f.rt() >= low && f.rt() <= high) {
                result.push_back(&f);
            }
        }
        return result;
    }

    /// Find features within both m/z and RT range
    [[nodiscard]] std::vector<const Feature*> findInRange(
        MZ mz_low, MZ mz_high, RetentionTime rt_low, RetentionTime rt_high) const {
        std::vector<const Feature*> result;
        for (const auto& f : features_) {
            if (f.mz() >= mz_low && f.mz() <= mz_high &&
                f.rt() >= rt_low && f.rt() <= rt_high) {
                result.push_back(&f);
            }
        }
        return result;
    }

    /// Find feature closest to m/z and RT
    [[nodiscard]] const Feature* findNearest(MZ mz, RetentionTime rt) const {
        if (empty()) return nullptr;

        const Feature* nearest = nullptr;
        double min_dist = std::numeric_limits<double>::max();

        for (const auto& f : features_) {
            // Normalize distances (approximate)
            double mz_dist = std::abs(f.mz() - mz) / (mz_range_.span() + 1e-10);
            double rt_dist = std::abs(f.rt() - rt) / (rt_range_.span() + 1e-10);
            double dist = std::sqrt(mz_dist * mz_dist + rt_dist * rt_dist);

            if (dist < min_dist) {
                min_dist = dist;
                nearest = &f;
            }
        }
        return nearest;
    }

    // =========================================================================
    // Ranges
    // =========================================================================

    /// Get overall m/z range
    [[nodiscard]] const MZRange& mzRange() const noexcept { return mz_range_; }

    /// Get overall RT range
    [[nodiscard]] const RTRange& rtRange() const noexcept { return rt_range_; }

    // =========================================================================
    // Metadata
    // =========================================================================

    /// Get source file path
    [[nodiscard]] const std::string& sourceFile() const noexcept {
        return source_file_;
    }
    void setSourceFile(std::string path) { source_file_ = std::move(path); }

    /// Get custom metadata
    [[nodiscard]] const MetaData& metadata() const noexcept { return metadata_; }
    MetaData& metadata() noexcept { return metadata_; }

    /// Get underlying vector
    [[nodiscard]] const std::vector<Feature>& features() const noexcept {
        return features_;
    }
    std::vector<Feature>& features() noexcept { return features_; }

private:
    void updateRanges(const Feature& f) {
        mz_range_.extend(f.mz());
        rt_range_.extend(f.rt());
    }

    std::vector<Feature> features_;
    MZRange mz_range_;
    RTRange rt_range_;
    Index next_id_ = 0;
    std::string source_file_;
    MetaData metadata_;
};

} // namespace lcms
