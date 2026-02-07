#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "lcms/spectrum.hpp"

using namespace lcms;
using Catch::Approx;

TEST_CASE("Spectrum construction", "[spectrum]") {
    SECTION("Default construction") {
        Spectrum spec;
        REQUIRE(spec.size() == 0);
        REQUIRE(spec.empty());
        REQUIRE(spec.msLevel() == 1);
    }

    SECTION("Construction with data") {
        std::vector<MZ> mz = {100.0, 200.0, 300.0};
        std::vector<Intensity> intensity = {1000.0, 5000.0, 2000.0};

        Spectrum spec(mz, intensity);
        REQUIRE(spec.size() == 3);
        REQUIRE_FALSE(spec.empty());
        REQUIRE(spec.mzAt(0) == Approx(100.0));
        REQUIRE(spec.intensityAt(1) == Approx(5000.0));
    }

    SECTION("Mismatched array sizes throw") {
        std::vector<MZ> mz = {100.0, 200.0};
        std::vector<Intensity> intensity = {1000.0};

        REQUIRE_THROWS_AS(Spectrum(mz, intensity), std::invalid_argument);
    }
}

TEST_CASE("Spectrum statistics", "[spectrum]") {
    std::vector<MZ> mz = {100.0, 200.0, 300.0, 400.0};
    std::vector<Intensity> intensity = {1000.0, 5000.0, 2000.0, 500.0};
    Spectrum spec(mz, intensity);

    SECTION("TIC calculation") {
        REQUIRE(spec.tic() == Approx(8500.0));
    }

    SECTION("Base peak") {
        REQUIRE(spec.basePeakMz() == Approx(200.0));
        REQUIRE(spec.basePeakIntensity() == Approx(5000.0));
    }

    SECTION("m/z range") {
        auto range = spec.mzRange();
        REQUIRE(range.min_value == Approx(100.0));
        REQUIRE(range.max_value == Approx(400.0));
    }
}

TEST_CASE("Spectrum operations", "[spectrum]") {
    std::vector<MZ> mz = {300.0, 100.0, 200.0};
    std::vector<Intensity> intensity = {3000.0, 1000.0, 2000.0};
    Spectrum spec(mz, intensity);

    SECTION("Sort by m/z") {
        REQUIRE_FALSE(spec.isSortedByMz());
        spec.sortByMz();
        REQUIRE(spec.isSortedByMz());
        REQUIRE(spec.mzAt(0) == Approx(100.0));
        REQUIRE(spec.mzAt(1) == Approx(200.0));
        REQUIRE(spec.mzAt(2) == Approx(300.0));
        REQUIRE(spec.intensityAt(0) == Approx(1000.0));
    }

    SECTION("Find nearest m/z") {
        spec.sortByMz();
        REQUIRE(spec.findNearestMz(150.0) == 1);  // Closer to 100 than 200
        REQUIRE(spec.findNearestMz(199.0) == 1);  // Closer to 200
        REQUIRE(spec.findNearestMz(201.0) == 1);  // Closer to 200
    }

    SECTION("Extract range") {
        spec.sortByMz();
        auto extracted = spec.extractRange(150.0, 250.0);
        REQUIRE(extracted.size() == 1);
        REQUIRE(extracted.mzAt(0) == Approx(200.0));
    }
}

TEST_CASE("Spectrum metadata", "[spectrum]") {
    Spectrum spec;

    SECTION("MS level") {
        spec.setMsLevel(2);
        REQUIRE(spec.msLevel() == 2);
    }

    SECTION("Retention time") {
        spec.setRetentionTime(123.45);
        REQUIRE(spec.retentionTime() == Approx(123.45));
    }

    SECTION("Spectrum type") {
        spec.setType(SpectrumType::CENTROID);
        REQUIRE(spec.type() == SpectrumType::CENTROID);
    }

    SECTION("Polarity") {
        spec.setPolarity(Polarity::POSITIVE);
        REQUIRE(spec.polarity() == Polarity::POSITIVE);
    }

    SECTION("Precursors") {
        Precursor prec;
        prec.mz = 500.0;
        prec.charge = 2;
        spec.addPrecursor(prec);

        REQUIRE(spec.precursors().size() == 1);
        REQUIRE(spec.precursors()[0].mz == Approx(500.0));
        REQUIRE(spec.precursors()[0].charge == 2);
    }
}
