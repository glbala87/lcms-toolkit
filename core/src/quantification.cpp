#include "lcms/algorithms/quantification.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace lcms {
namespace algorithms {

ConsensusMap ConsensusMap::filterByPresence(double min_fraction) const {
    ConsensusMap result;
    result.setSampleNames(sample_names_);

    size_t n_samples = numSamples();
    for (size_t i = 0; i < numFeatures(); ++i) {
        int present = 0;
        for (size_t j = 0; j < n_samples; ++j) {
            if (matrix_[i][j] > 0) present++;
        }
        double fraction = static_cast<double>(present) / n_samples;
        if (fraction >= min_fraction) {
            result.addFeature(feature_ids_[i], matrix_[i]);
        }
    }

    return result;
}

ConsensusMap ConsensusMap::fillMissing() const {
    ConsensusMap result;
    result.setSampleNames(sample_names_);

    auto filled_matrix = matrix_;
    size_t n_samples = numSamples();

    for (size_t j = 0; j < n_samples; ++j) {
        double col_min = std::numeric_limits<double>::max();
        for (size_t i = 0; i < numFeatures(); ++i) {
            if (filled_matrix[i][j] > 0 && filled_matrix[i][j] < col_min) {
                col_min = filled_matrix[i][j];
            }
        }
        if (col_min == std::numeric_limits<double>::max()) col_min = 1.0;

        for (size_t i = 0; i < numFeatures(); ++i) {
            if (filled_matrix[i][j] == 0) {
                filled_matrix[i][j] = col_min;
            }
        }
    }

    for (size_t i = 0; i < numFeatures(); ++i) {
        result.addFeature(feature_ids_[i], filled_matrix[i]);
    }

    return result;
}

ConsensusMap ConsensusMap::medianNormalize() const {
    auto normalized = medianNormalization(matrix_);
    ConsensusMap result(feature_ids_, sample_names_, normalized);
    return result;
}

std::vector<std::vector<double>> medianNormalization(
    const std::vector<std::vector<double>>& matrix) {

    if (matrix.empty()) return matrix;

    auto result = matrix;
    size_t n_features = matrix.size();
    size_t n_samples = matrix[0].size();

    std::vector<double> medians(n_samples);
    for (size_t j = 0; j < n_samples; ++j) {
        std::vector<double> nonzero;
        for (size_t i = 0; i < n_features; ++i) {
            if (matrix[i][j] > 0) nonzero.push_back(matrix[i][j]);
        }
        if (nonzero.empty()) {
            medians[j] = 1.0;
        } else {
            std::sort(nonzero.begin(), nonzero.end());
            size_t mid = nonzero.size() / 2;
            medians[j] = (nonzero.size() % 2 == 0)
                ? (nonzero[mid - 1] + nonzero[mid]) / 2.0
                : nonzero[mid];
        }
    }

    // Target median
    std::vector<double> valid_medians;
    for (double m : medians) {
        if (m > 0) valid_medians.push_back(m);
    }
    std::sort(valid_medians.begin(), valid_medians.end());
    double target = valid_medians.empty() ? 1.0 :
        valid_medians[valid_medians.size() / 2];

    for (size_t j = 0; j < n_samples; ++j) {
        if (medians[j] > 0) {
            double scale = target / medians[j];
            for (size_t i = 0; i < n_features; ++i) {
                result[i][j] *= scale;
            }
        }
    }

    return result;
}

std::vector<std::vector<double>> quantileNormalization(
    const std::vector<std::vector<double>>& matrix) {

    if (matrix.empty()) return matrix;

    size_t n_features = matrix.size();
    size_t n_samples = matrix[0].size();

    auto result = matrix;

    // Sort each column, compute row means
    std::vector<std::vector<double>> sorted_cols(n_samples);
    std::vector<std::vector<size_t>> ranks(n_samples);

    for (size_t j = 0; j < n_samples; ++j) {
        std::vector<size_t> order(n_features);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](size_t a, size_t b) { return matrix[a][j] < matrix[b][j]; });

        sorted_cols[j].resize(n_features);
        ranks[j].resize(n_features);
        for (size_t i = 0; i < n_features; ++i) {
            sorted_cols[j][i] = matrix[order[i]][j];
            ranks[j][order[i]] = i;
        }
    }

    // Row means of sorted values
    std::vector<double> row_means(n_features, 0.0);
    for (size_t i = 0; i < n_features; ++i) {
        for (size_t j = 0; j < n_samples; ++j) {
            row_means[i] += sorted_cols[j][i];
        }
        row_means[i] /= n_samples;
    }

    // Replace with row means
    for (size_t j = 0; j < n_samples; ++j) {
        for (size_t i = 0; i < n_features; ++i) {
            result[i][j] = row_means[ranks[j][i]];
        }
    }

    return result;
}

} // namespace algorithms
} // namespace lcms
