#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "lcms/io/base64.hpp"

using namespace lcms::io;
using Catch::Approx;

TEST_CASE("Base64 encoding", "[base64]") {
    SECTION("Encode empty") {
        std::vector<std::uint8_t> data;
        REQUIRE(Base64::encode(data).empty());
    }

    SECTION("Encode simple") {
        std::vector<std::uint8_t> data = {'H', 'e', 'l', 'l', 'o'};
        std::string encoded = Base64::encode(data);
        REQUIRE(encoded == "SGVsbG8=");
    }

    SECTION("Encode with padding") {
        std::vector<std::uint8_t> data = {'a'};
        std::string encoded = Base64::encode(data);
        REQUIRE(encoded == "YQ==");

        data = {'a', 'b'};
        encoded = Base64::encode(data);
        REQUIRE(encoded == "YWI=");

        data = {'a', 'b', 'c'};
        encoded = Base64::encode(data);
        REQUIRE(encoded == "YWJj");
    }
}

TEST_CASE("Base64 decoding", "[base64]") {
    SECTION("Decode empty") {
        auto decoded = Base64::decode("");
        REQUIRE(decoded.empty());
    }

    SECTION("Decode simple") {
        auto decoded = Base64::decode("SGVsbG8=");
        REQUIRE(decoded.size() == 5);
        REQUIRE(decoded[0] == 'H');
        REQUIRE(decoded[4] == 'o');
    }

    SECTION("Decode with whitespace") {
        auto decoded = Base64::decode("SGVs\nbG8=");
        REQUIRE(decoded.size() == 5);
    }
}

TEST_CASE("Base64 roundtrip", "[base64]") {
    SECTION("Binary data roundtrip") {
        std::vector<std::uint8_t> original = {0, 1, 127, 128, 255};
        std::string encoded = Base64::encode(original);
        auto decoded = Base64::decode(encoded);

        REQUIRE(decoded == original);
    }
}

TEST_CASE("Base64 float decoding", "[base64]") {
    SECTION("Decode float64 little endian") {
        // 1.0 as little-endian double: 0x3FF0000000000000
        // Base64 encoded
        std::string encoded = "AAAAAAAA8D8=";
        auto values = Base64::decodeFloat64(encoded, true);
        REQUIRE(values.size() == 1);
        REQUIRE(values[0] == Approx(1.0));
    }

    SECTION("Decode multiple float64") {
        // 1.0 and 2.0 as little-endian doubles
        std::string encoded = "AAAAAAAA8D8AAAAAAAAAQAAc";
        auto values = Base64::decodeFloat64(encoded, true);
        REQUIRE(values.size() == 2);
        REQUIRE(values[0] == Approx(1.0));
        REQUIRE(values[1] == Approx(2.0));
    }
}

TEST_CASE("Base64 validation", "[base64]") {
    SECTION("Valid Base64") {
        REQUIRE(Base64::isValid("SGVsbG8="));
        REQUIRE(Base64::isValid("YWJj"));
        REQUIRE(Base64::isValid(""));
    }

    SECTION("Invalid Base64") {
        REQUIRE_FALSE(Base64::isValid("Hello!"));
        REQUIRE_FALSE(Base64::isValid("Test#$%"));
    }
}

TEST_CASE("Base64 decoded size", "[base64]") {
    SECTION("Calculate decoded size") {
        REQUIRE(Base64::decodedSize("YQ==") == 1);
        REQUIRE(Base64::decodedSize("YWI=") == 2);
        REQUIRE(Base64::decodedSize("YWJj") == 3);
        REQUIRE(Base64::decodedSize("YWJjZA==") == 4);
    }
}
