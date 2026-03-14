/**
 * Tests for isotope detection and deconvolution (C++).
 * Uses Catch2 testing framework.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>
#include <cmath>

#include "lcms/algorithms/isotope_detection.hpp"

using namespace lcms;
using Catch::Matchers::WithinAbs;

TEST_CASE("Averagine model generates valid distributions", "[isotope]") {
    SECTION("Small mass") {
        auto dist = averagineDistribution(1000.0, 5);
        REQUIRE(dist.size() == 5);
        double sum = 0;
        for (auto d : dist) sum += d;
        REQUIRE_THAT(sum, WithinAbs(1.0, 0.05));
        REQUIRE(dist[0] > dist[4]);
    }

    SECTION("Large mass") {
        auto dist = averagineDistribution(5000.0, 8);
        REQUIRE(dist.size() == 8);
        double sum = 0;
        for (auto d : dist) sum += d;
        REQUIRE(sum > 0.9);
    }
}

TEST_CASE("Isotope pattern detection", "[isotope]") {
    // Create synthetic isotope pattern at z=2
    double monoMz = 500.0;
    int charge = 2;
    double spacing = 1.003355 / charge;

    std::vector<double> mz, intensity;
    auto dist = averagineDistribution(monoMz * charge, 5);
    for (int i = 0; i < 5; i++) {
        mz.push_back(monoMz + i * spacing);
        intensity.push_back(dist[i] * 10000);
    }
    // Add noise
    for (int i = 0; i < 50; i++) {
        mz.push_back(400.0 + i * 4.0);
        intensity.push_back(50.0);
    }

    auto patterns = detectIsotopePatterns(mz, intensity);
    REQUIRE(patterns.size() >= 1);
}

TEST_CASE("Charge state assignment", "[isotope]") {
    std::vector<double> mz = {500.0, 500.5016, 501.0032};
    int charge = assignChargeState(mz);
    REQUIRE(charge == 2);
}
