#pragma once

#include "../types.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace lcms {
namespace algorithms {

/**
 * @brief A consensus feature aligned across multiple samples.
 */
struct ConsensusFeature {
    std::string id;
    MZ mz = 0.0;
    RetentionTime rt = 0.0;
    std::vector<Intensity> intensities; // One per sample
};

/**
 * @brief Matrix of aligned features across samples.
 */
class ConsensusMap {
public:
    ConsensusMap() = default;

    ConsensusMap(const std::vector<std::string>& feature_ids,
                 const std::vector<std::string>& sample_names,
                 const std::vector<std::vector<double>>& matrix)
        : feature_ids_(feature_ids), sample_names_(sample_names), matrix_(matrix) {}

    [[nodiscard]] size_t numFeatures() const { return feature_ids_.size(); }
    [[nodiscard]] size_t numSamples() const { return sample_names_.size(); }

    [[nodiscard]] const std::vector<std::string>& featureIds() const { return feature_ids_; }
    [[nodiscard]] const std::vector<std::string>& sampleNames() const { return sample_names_; }

    /// Get intensity for feature i, sample j
    [[nodiscard]] double intensity(size_t i, size_t j) const {
        return matrix_[i][j];
    }

    /// Set intensity
    void setIntensity(size_t i, size_t j, double value) {
        matrix_[i][j] = value;
    }

    /// Get entire row (feature across samples)
    [[nodiscard]] const std::vector<double>& featureRow(size_t i) const {
        return matrix_[i];
    }

    /// Add a feature
    void addFeature(const std::string& id, const std::vector<double>& intensities) {
        feature_ids_.push_back(id);
        matrix_.push_back(intensities);
    }

    void setSampleNames(const std::vector<std::string>& names) {
        sample_names_ = names;
    }

    /// Filter by minimum presence fraction
    ConsensusMap filterByPresence(double min_fraction) const;

    /// Fill missing (zero) values with column minimum
    ConsensusMap fillMissing() const;

    /// Median normalization
    ConsensusMap medianNormalize() const;

private:
    std::vector<std::string> feature_ids_;
    std::vector<std::string> sample_names_;
    std::vector<std::vector<double>> matrix_; // features x samples
};

/**
 * @brief Align features across multiple runs.
 */
struct FeatureAlignmentOptions {
    double mz_tolerance = 0.01;
    double rt_tolerance = 30.0;
};

/**
 * @brief Median normalization of intensity matrix.
 */
std::vector<std::vector<double>> medianNormalization(
    const std::vector<std::vector<double>>& matrix);

/**
 * @brief Quantile normalization of intensity matrix.
 */
std::vector<std::vector<double>> quantileNormalization(
    const std::vector<std::vector<double>>& matrix);

} // namespace algorithms
} // namespace lcms
