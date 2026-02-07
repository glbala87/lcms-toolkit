#pragma once

#include "types.hpp"
#include "spectrum.hpp"
#include "chromatogram.hpp"
#include <vector>
#include <map>
#include <memory>
#include <algorithm>
#include <functional>

namespace lcms {

/**
 * @brief Container for a complete LC-MS experiment.
 *
 * MSExperiment holds all spectra and chromatograms from an LC-MS run,
 * along with instrument metadata and run information.
 */
class MSExperiment {
public:
    using SpectrumIterator = std::vector<Spectrum>::iterator;
    using SpectrumConstIterator = std::vector<Spectrum>::const_iterator;
    using ChromatogramIterator = std::vector<Chromatogram>::iterator;
    using ChromatogramConstIterator = std::vector<Chromatogram>::const_iterator;

    MSExperiment() = default;

    // Move/copy operations
    MSExperiment(MSExperiment&&) noexcept = default;
    MSExperiment& operator=(MSExperiment&&) noexcept = default;
    MSExperiment(const MSExperiment&) = default;
    MSExperiment& operator=(const MSExperiment&) = default;

    // =========================================================================
    // Spectra Access
    // =========================================================================

    /// Get number of spectra
    [[nodiscard]] std::size_t spectrumCount() const noexcept {
        return spectra_.size();
    }

    /// Check if experiment has spectra
    [[nodiscard]] bool hasSpectra() const noexcept { return !spectra_.empty(); }

    /// Access spectrum by index
    [[nodiscard]] Spectrum& spectrum(std::size_t i) { return spectra_.at(i); }
    [[nodiscard]] const Spectrum& spectrum(std::size_t i) const {
        return spectra_.at(i);
    }

    /// Operator access to spectra
    [[nodiscard]] Spectrum& operator[](std::size_t i) { return spectra_[i]; }
    [[nodiscard]] const Spectrum& operator[](std::size_t i) const {
        return spectra_[i];
    }

    /// Get all spectra
    [[nodiscard]] const std::vector<Spectrum>& spectra() const noexcept {
        return spectra_;
    }
    std::vector<Spectrum>& spectra() noexcept { return spectra_; }

    /// Spectrum iterators
    SpectrumIterator spectraBegin() noexcept { return spectra_.begin(); }
    SpectrumIterator spectraEnd() noexcept { return spectra_.end(); }
    SpectrumConstIterator spectraBegin() const noexcept {
        return spectra_.begin();
    }
    SpectrumConstIterator spectraEnd() const noexcept { return spectra_.end(); }

    /// Add a spectrum
    void addSpectrum(Spectrum s) {
        s.setIndex(spectra_.size());
        updateRangesFromSpectrum(s);
        spectra_.push_back(std::move(s));
    }

    /// Reserve capacity for spectra
    void reserveSpectra(std::size_t n) { spectra_.reserve(n); }

    // =========================================================================
    // Chromatograms Access
    // =========================================================================

    /// Get number of chromatograms
    [[nodiscard]] std::size_t chromatogramCount() const noexcept {
        return chromatograms_.size();
    }

    /// Check if experiment has chromatograms
    [[nodiscard]] bool hasChromatograms() const noexcept {
        return !chromatograms_.empty();
    }

    /// Access chromatogram by index
    [[nodiscard]] Chromatogram& chromatogram(std::size_t i) {
        return chromatograms_.at(i);
    }
    [[nodiscard]] const Chromatogram& chromatogram(std::size_t i) const {
        return chromatograms_.at(i);
    }

    /// Get all chromatograms
    [[nodiscard]] const std::vector<Chromatogram>& chromatograms() const noexcept {
        return chromatograms_;
    }
    std::vector<Chromatogram>& chromatograms() noexcept { return chromatograms_; }

    /// Chromatogram iterators
    ChromatogramIterator chromatogramsBegin() noexcept {
        return chromatograms_.begin();
    }
    ChromatogramIterator chromatogramsEnd() noexcept {
        return chromatograms_.end();
    }
    ChromatogramConstIterator chromatogramsBegin() const noexcept {
        return chromatograms_.begin();
    }
    ChromatogramConstIterator chromatogramsEnd() const noexcept {
        return chromatograms_.end();
    }

    /// Add a chromatogram
    void addChromatogram(Chromatogram c) {
        c.setIndex(chromatograms_.size());
        chromatograms_.push_back(std::move(c));
    }

    /// Reserve capacity for chromatograms
    void reserveChromatograms(std::size_t n) { chromatograms_.reserve(n); }

    // =========================================================================
    // Filtering and Selection
    // =========================================================================

    /// Get spectra by MS level
    [[nodiscard]] std::vector<const Spectrum*> getSpectraByLevel(
        MSLevel level) const {
        std::vector<const Spectrum*> result;
        for (const auto& s : spectra_) {
            if (s.msLevel() == level) {
                result.push_back(&s);
            }
        }
        return result;
    }

    /// Count spectra at given MS level
    [[nodiscard]] std::size_t countSpectraByLevel(MSLevel level) const {
        return std::count_if(spectra_.begin(), spectra_.end(),
            [level](const Spectrum& s) { return s.msLevel() == level; });
    }

    /// Get spectra within RT range
    [[nodiscard]] std::vector<const Spectrum*> getSpectraInRtRange(
        RetentionTime low, RetentionTime high) const {
        std::vector<const Spectrum*> result;
        for (const auto& s : spectra_) {
            if (s.retentionTime() >= low && s.retentionTime() <= high) {
                result.push_back(&s);
            }
        }
        return result;
    }

    /// Find spectrum by native ID
    [[nodiscard]] const Spectrum* findSpectrumByNativeId(
        const std::string& id) const {
        for (const auto& s : spectra_) {
            if (s.nativeId() == id) {
                return &s;
            }
        }
        return nullptr;
    }

    /// Find spectrum closest to given RT
    [[nodiscard]] const Spectrum* findSpectrumByRt(
        RetentionTime rt, MSLevel level = 0) const {
        const Spectrum* nearest = nullptr;
        double min_dist = std::numeric_limits<double>::max();

        for (const auto& s : spectra_) {
            if (level > 0 && s.msLevel() != level) continue;

            double dist = std::abs(s.retentionTime() - rt);
            if (dist < min_dist) {
                min_dist = dist;
                nearest = &s;
            }
        }
        return nearest;
    }

    /// Get TIC chromatogram
    [[nodiscard]] const Chromatogram* getTICChromatogram() const {
        for (const auto& c : chromatograms_) {
            if (c.type() == ChromatogramType::TIC) {
                return &c;
            }
        }
        return nullptr;
    }

    // =========================================================================
    // Extracted Chromatograms
    // =========================================================================

    /// Generate TIC from spectra
    [[nodiscard]] Chromatogram generateTIC(MSLevel level = 1) const {
        Chromatogram tic;
        tic.setType(ChromatogramType::TIC);
        tic.setNativeId("TIC");

        for (const auto& s : spectra_) {
            if (level == 0 || s.msLevel() == level) {
                tic.addPoint(s.retentionTime(), s.tic());
            }
        }

        tic.sortByRt();
        return tic;
    }

    /// Generate base peak chromatogram from spectra
    [[nodiscard]] Chromatogram generateBPC(MSLevel level = 1) const {
        Chromatogram bpc;
        bpc.setType(ChromatogramType::BPC);
        bpc.setNativeId("BPC");

        for (const auto& s : spectra_) {
            if (level == 0 || s.msLevel() == level) {
                bpc.addPoint(s.retentionTime(), s.basePeakIntensity());
            }
        }

        bpc.sortByRt();
        return bpc;
    }

    /// Generate extracted ion chromatogram (XIC)
    [[nodiscard]] Chromatogram generateXIC(
        MZ target_mz, MZTolerance tolerance, MSLevel level = 1) const {
        Chromatogram xic;
        xic.setType(ChromatogramType::XIC);
        xic.setTargetMz(target_mz);
        xic.setMzTolerance(tolerance.absoluteAt(target_mz));

        for (const auto& s : spectra_) {
            if (level > 0 && s.msLevel() != level) continue;
            if (!s.isSortedByMz()) continue;

            double tol = tolerance.absoluteAt(target_mz);
            MZ low = target_mz - tol;
            MZ high = target_mz + tol;

            // Sum intensities in range
            Intensity sum = 0;
            for (std::size_t i = 0; i < s.size(); ++i) {
                if (s.mzAt(i) >= low && s.mzAt(i) <= high) {
                    sum += s.intensityAt(i);
                } else if (s.mzAt(i) > high) {
                    break;  // Sorted, so can stop early
                }
            }

            xic.addPoint(s.retentionTime(), sum);
        }

        xic.sortByRt();
        return xic;
    }

    // =========================================================================
    // Ranges
    // =========================================================================

    /// Get overall m/z range
    [[nodiscard]] const MZRange& mzRange() const noexcept { return mz_range_; }

    /// Get overall RT range
    [[nodiscard]] const RTRange& rtRange() const noexcept { return rt_range_; }

    /// Update all ranges
    void updateRanges() {
        mz_range_ = MZRange();
        rt_range_ = RTRange();
        for (const auto& s : spectra_) {
            updateRangesFromSpectrum(s);
        }
    }

    // =========================================================================
    // Sorting
    // =========================================================================

    /// Sort spectra by retention time
    void sortSpectraByRt() {
        std::sort(spectra_.begin(), spectra_.end(),
            [](const Spectrum& a, const Spectrum& b) {
                return a.retentionTime() < b.retentionTime();
            });
        // Update indices
        for (std::size_t i = 0; i < spectra_.size(); ++i) {
            spectra_[i].setIndex(i);
        }
    }

    /// Sort spectra by MS level, then RT
    void sortSpectraByLevelAndRt() {
        std::sort(spectra_.begin(), spectra_.end(),
            [](const Spectrum& a, const Spectrum& b) {
                if (a.msLevel() != b.msLevel()) {
                    return a.msLevel() < b.msLevel();
                }
                return a.retentionTime() < b.retentionTime();
            });
        for (std::size_t i = 0; i < spectra_.size(); ++i) {
            spectra_[i].setIndex(i);
        }
    }

    // =========================================================================
    // Clear/Reset
    // =========================================================================

    /// Clear all data
    void clear() {
        spectra_.clear();
        chromatograms_.clear();
        mz_range_ = MZRange();
        rt_range_ = RTRange();
        metadata_.clear();
    }

    /// Clear only spectra
    void clearSpectra() {
        spectra_.clear();
        mz_range_ = MZRange();
        rt_range_ = RTRange();
    }

    /// Clear only chromatograms
    void clearChromatograms() {
        chromatograms_.clear();
    }

    // =========================================================================
    // Metadata
    // =========================================================================

    /// Get source file path
    [[nodiscard]] const std::string& sourceFile() const noexcept {
        return source_file_;
    }
    void setSourceFile(std::string path) { source_file_ = std::move(path); }

    /// Get experiment date/time
    [[nodiscard]] const std::string& dateTime() const noexcept {
        return date_time_;
    }
    void setDateTime(std::string dt) { date_time_ = std::move(dt); }

    /// Get instrument model
    [[nodiscard]] const std::string& instrumentModel() const noexcept {
        return instrument_model_;
    }
    void setInstrumentModel(std::string model) {
        instrument_model_ = std::move(model);
    }

    /// Get instrument serial number
    [[nodiscard]] const std::string& instrumentSerial() const noexcept {
        return instrument_serial_;
    }
    void setInstrumentSerial(std::string serial) {
        instrument_serial_ = std::move(serial);
    }

    /// Get software name
    [[nodiscard]] const std::string& software() const noexcept {
        return software_;
    }
    void setSoftware(std::string sw) { software_ = std::move(sw); }

    /// Get custom metadata
    [[nodiscard]] const MetaData& metadata() const noexcept { return metadata_; }
    MetaData& metadata() noexcept { return metadata_; }

    // =========================================================================
    // Statistics
    // =========================================================================

    /// Get total number of data points across all spectra
    [[nodiscard]] std::size_t totalDataPoints() const {
        std::size_t total = 0;
        for (const auto& s : spectra_) {
            total += s.size();
        }
        return total;
    }

    /// Get average spectrum size
    [[nodiscard]] double averageSpectrumSize() const {
        if (spectra_.empty()) return 0.0;
        return static_cast<double>(totalDataPoints()) / spectra_.size();
    }

    /// Get memory usage estimate (bytes)
    [[nodiscard]] std::size_t estimatedMemoryUsage() const {
        std::size_t bytes = 0;
        for (const auto& s : spectra_) {
            bytes += s.size() * (sizeof(MZ) + sizeof(Intensity));
        }
        for (const auto& c : chromatograms_) {
            bytes += c.size() * (sizeof(RetentionTime) + sizeof(Intensity));
        }
        return bytes;
    }

private:
    void updateRangesFromSpectrum(const Spectrum& s) {
        mz_range_.extend(s.mzRange());
        rt_range_.extend(s.retentionTime());
    }

    std::vector<Spectrum> spectra_;
    std::vector<Chromatogram> chromatograms_;

    MZRange mz_range_;
    RTRange rt_range_;

    std::string source_file_;
    std::string date_time_;
    std::string instrument_model_;
    std::string instrument_serial_;
    std::string software_;

    MetaData metadata_;
};

} // namespace lcms
