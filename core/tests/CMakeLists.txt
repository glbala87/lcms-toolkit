include(FetchContent)

# Fetch Catch2 for testing
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(Catch2)

# Test executable
add_executable(lcms_tests
    test_spectrum.cpp
    test_chromatogram.cpp
    test_peak.cpp
    test_base64.cpp
)

target_link_libraries(lcms_tests
    PRIVATE
        lcms_core
        Catch2::Catch2WithMain
)

# Register tests with CTest
include(Catch)
catch_discover_tests(lcms_tests)
