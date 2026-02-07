#pragma once

#include "types.hpp"
#include <vector>
#include <cmath>
#include <optional>

namespace lcms {

/**
 * @brief Represents a detected peak in a spectrum or chromatogram.
 *
 * A Peak stores information about a local maximum including its position
 * (m/z or RT), intensity, area, shape parameters, and quality metrics.
 */
class Peak {
public:
    /// Default constructor
    Peak() = default;

    /// Construct with basic parameters
    Peak(MZ mz, Intensity intensity, RetentionTime rt = 0.0)
        : mz_(mz), intensity_(intensity), rt_(rt) {}

    // =========================================================================
    // Position and Intensity
    // =========================================================================

    /// Get m/z position
    [[nodiscard]] MZ mz() const noexcept { return mz_; }
    void setMz(MZ mz) noexcept { mz_ = mz; }

    /// Get retention time
    [[nodiscard]] RetentionTime rt() const noexcept { return rt_; }
    void setRt(RetentionTime rt) noexcept { rt_ = rt; }

    /// Get peak intensity (apex height)
    [[nodiscard]] Intensity intensity() const noexcept { return intensity_; }
    void setIntensity(Intensity intensity) noexcept { intensity_ = intensity; }

    /// Get integrated peak area
    [[nodiscard]] double area() const noexcept { return area_; }
    void setArea(double area) noexcept { area_ = area; }

    // =========================================================================
    // Boundaries
    // =========================================================================

    /// Get left boundary (m/z)
    [[nodiscard]] MZ mzLeft() const noexcept { return mz_left_; }
    void setMzLeft(MZ mz) noexcept { mz_left_ = mz; }

    /// Get right boundary (m/z)
    [[nodiscard]] MZ mzRight() const noexcept { return mz_right_; }
    void setMzRight(MZ mz) noexcept { mz_right_ = mz; }

    /// Get left boundary (RT)
    [[nodiscard]] RetentionTime rtLeft() const noexcept { return rt_left_; }
    void setRtLeft(RetentionTime rt) noexcept { rt_left_ = rt; }

    /// Get right boundary (RT)
    [[nodiscard]] RetentionTime rtRight() const noexcept { return rt_right_; }
    void setRtRight(RetentionTime rt) noexcept { rt_right_ = rt; }

    /// Get m/z width (FWHM or total width)
    [[nodiscard]] MZ mzWidth() const noexcept { return mz_right_ - mz_left_; }

    /// Get RT width
    [[nodiscard]] RetentionTime rtWidth() const noexcept {
        return rt_right_ - rt_left_;
    }

    // =========================================================================
    // Shape Parameters
    // =========================================================================

    /// Get full width at half maximum (m/z)
    [[nodiscard]] MZ fwhmMz() const noexcept { return fwhm_mz_; }
    void setFwhmMz(MZ fwhm) noexcept { fwhm_mz_ = fwhm; }

    /// Get full width at half maximum (RT)
    [[nodiscard]] RetentionTime fwhmRt() const noexcept { return fwhm_rt_; }
    void setFwhmRt(RetentionTime fwhm) noexcept { fwhm_rt_ = fwhm; }

    /// Get peak shape type
    [[nodiscard]] PeakShape shape() const noexcept { return shape_; }
    void setShape(PeakShape shape) noexcept { shape_ = shape; }

    /// Get asymmetry factor (right width / left width at 10% height)
    [[nodiscard]] double asymmetry() const noexcept { return asymmetry_; }
    void setAsymmetry(double asym) noexcept { asymmetry_ = asym; }

    /// Get Gaussian sigma (if Gaussian fit performed)
    [[nodiscard]] std::optional<double> sigma() const noexcept { return sigma_; }
    void setSigma(double s) noexcept { sigma_ = s; }

    // =========================================================================
    // Charge and Isotopes
    // =========================================================================

    /// Get charge state
    [[nodiscard]] ChargeState charge() const noexcept { return charge_; }
    void setCharge(ChargeState z) noexcept { charge_ = z; }

    /// Check if charge state is known
    [[nodiscard]] bool hasCharge() const noexcept { return charge_ != 0; }

    /// Get monoisotopic m/z (if determined)
    [[nodiscard]] std::optional<MZ> monoisotopicMz() const noexcept {
        return mono_mz_;
    }
    void setMonoisotopicMz(MZ mz) noexcept { mono_mz_ = mz; }

    /// Get isotope index (0 = monoisotopic)
    [[nodiscard]] int isotopeIndex() const noexcept { return isotope_index_; }
    void setIsotopeIndex(int idx) noexcept { isotope_index_ = idx; }

    // =========================================================================
    // Quality Metrics
    // =========================================================================

    /// Get signal-to-noise ratio
    [[nodiscard]] double snr() const noexcept { return snr_; }
    void setSnr(double snr) noexcept { snr_ = snr; }

    /// Get overall quality score (0-1)
    [[nodiscard]] double quality() const noexcept { return quality_; }
    void setQuality(double q) noexcept { quality_ = q; }

    /// Get fit residual (if model fit performed)
    [[nodiscard]] std::optional<double> fitResidual() const noexcept {
        return fit_residual_;
    }
    void setFitResidual(double r) noexcept { fit_residual_ = r; }

    // =========================================================================
    // Metadata
    // =========================================================================

    /// Get spectrum index where peak was detected
    [[nodiscard]] Index spectrumIndex() const noexcept { return spectrum_index_; }
    void setSpectrumIndex(Index idx) noexcept { spectrum_index_ = idx; }

    /// Get custom metadata
    [[nodiscard]] const MetaData& metadata() const noexcept { return metadata_; }
    MetaData& metadata() noexcept { return metadata_; }

    // =========================================================================
    // Utility Functions
    // =========================================================================

    /// Compute neutral mass from m/z and charge
    [[nodiscard]] double neutralMass() const {
        if (charge_ == 0) return 0.0;
        constexpr double proton_mass = 1.007276;
        return (mz_ - proton_mass) * std::abs(charge_);
    }

    /// Check if m/z is within peak boundaries
    [[nodiscard]] bool containsMz(MZ mz) const noexcept {
        return mz >= mz_left_ && mz <= mz_right_;
    }

    /// Check if RT is within peak boundaries
    [[nodiscard]] bool containsRt(RetentionTime rt) const noexcept {
        return rt >= rt_left_ && rt <= rt_right_;
    }

    /// Check if point is within peak boundaries (2D)
    [[nodiscard]] bool contains(MZ mz, RetentionTime rt) const noexcept {
        return containsMz(mz) && containsRt(rt);
    }

    /// Compare peaks by m/z
    [[nodiscard]] bool operator<(const Peak& other) const noexcept {
        return mz_ < other.mz_;
    }

private:
    // Position
    MZ mz_ = 0.0;
    RetentionTime rt_ = 0.0;
    Intensity intensity_ = 0.0;
    double area_ = 0.0;

    // Boundaries
    MZ mz_left_ = 0.0;
    MZ mz_right_ = 0.0;
    RetentionTime rt_left_ = 0.0;
    RetentionTime rt_right_ = 0.0;

    // Shape
    MZ fwhm_mz_ = 0.0;
    RetentionTime fwhm_rt_ = 0.0;
    PeakShape shape_ = PeakShape::UNKNOWN;
    double asymmetry_ = 1.0;
    std::optional<double> sigma_;

    // Charge/isotope
    ChargeState charge_ = 0;
    std::optional<MZ> mono_mz_;
    int isotope_index_ = 0;

    // Quality
    double snr_ = 0.0;
    double quality_ = 0.0;
    std::optional<double> fit_residual_;

    // Metadata
    Index spectrum_index_ = 0;
    MetaData metadata_;
};

/**
 * @brief Collection of peaks with spatial indexing capabilities.
 */
class PeakList {
public:
    using iterator = std::vector<Peak>::iterator;
    using const_iterator = std::vector<Peak>::const_iterator;

    PeakList() = default;

    /// Get number of peaks
    [[nodiscard]] std::size_t size() const noexcept { return peaks_.size(); }

    /// Check if empty
    [[nodiscard]] bool empty() const noexcept { return peaks_.empty(); }

    /// Access peak by index
    [[nodiscard]] Peak& operator[](std::size_t i) { return peaks_[i]; }
    [[nodiscard]] const Peak& operator[](std::size_t i) const { return peaks_[i]; }

    /// Iterator access
    iterator begin() noexcept { return peaks_.begin(); }
    iterator end() noexcept { return peaks_.end(); }
    const_iterator begin() const noexcept { return peaks_.begin(); }
    const_iterator end() const noexcept { return peaks_.end(); }
    const_iterator cbegin() const noexcept { return peaks_.cbegin(); }
    const_iterator cend() const noexcept { return peaks_.cend(); }

    /// Add a peak
    void add(Peak peak) { peaks_.push_back(std::move(peak)); }

    /// Add peak with basic parameters
    void add(MZ mz, Intensity intensity, RetentionTime rt = 0.0) {
        peaks_.emplace_back(mz, intensity, rt);
    }

    /// Reserve capacity
    void reserve(std::size_t n) { peaks_.reserve(n); }

    /// Clear all peaks
    void clear() { peaks_.clear(); }

    /// Sort by m/z
    void sortByMz() {
        std::sort(peaks_.begin(), peaks_.end(),
            [](const Peak& a, const Peak& b) { return a.mz() < b.mz(); });
    }

    /// Sort by RT
    void sortByRt() {
        std::sort(peaks_.begin(), peaks_.end(),
            [](const Peak& a, const Peak& b) { return a.rt() < b.rt(); });
    }

    /// Sort by intensity (descending)
    void sortByIntensity() {
        std::sort(peaks_.begin(), peaks_.end(),
            [](const Peak& a, const Peak& b) {
                return a.intensity() > b.intensity();
            });
    }

    /// Find peaks within m/z range
    [[nodiscard]] PeakList findInMzRange(MZ low, MZ high) const {
        PeakList result;
        for (const auto& p : peaks_) {
            if (p.mz() >= low && p.mz() <= high) {
                result.add(p);
            }
        }
        return result;
    }

    /// Find peaks within RT range
    [[nodiscard]] PeakList findInRtRange(RetentionTime low,
                                          RetentionTime high) const {
        PeakList result;
        for (const auto& p : peaks_) {
            if (p.rt() >= low && p.rt() <= high) {
                result.add(p);
            }
        }
        return result;
    }

    /// Find peak closest to m/z
    [[nodiscard]] const Peak* findNearestMz(MZ target) const {
        if (empty()) return nullptr;

        const Peak* nearest = &peaks_[0];
        double min_dist = std::abs(peaks_[0].mz() - target);

        for (const auto& p : peaks_) {
            double dist = std::abs(p.mz() - target);
            if (dist < min_dist) {
                min_dist = dist;
                nearest = &p;
            }
        }
        return nearest;
    }

    /// Get underlying vector
    [[nodiscard]] const std::vector<Peak>& peaks() const noexcept {
        return peaks_;
    }
    std::vector<Peak>& peaks() noexcept { return peaks_; }

private:
    std::vector<Peak> peaks_;
};

} // namespace lcms
