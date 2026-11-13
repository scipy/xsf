#include "../testing_utils.h"
#include <cmath>
#include <limits>
#include <tuple>
#include <xsf/cpu/specfun.h>
#include <xsf/hyp2f1.h>

TEST_CASE("hyp2f1(-1, b, c, z) / gh-4446", "[hyp2f1][xsf_tests]") {
    constexpr double rtol = 1e-14;
    constexpr double atol = 1e-14; // absolute tolerance
    constexpr int size = 41;
    constexpr double start = -10.0;
    constexpr double end = 10.0;

    const std::vector<double> test_values = linspace(start, end, size);

    for (double b : test_values) {
        for (double c : test_values) {
            if (c == 0.0)
                continue;

            for (double z : test_values) {
                double expected = 1 - (b / c) * z;
                double result = xsf::hyp2f1(-1.0, b, c, z);

                double error;
                if (std::abs(expected) < atol) {
                    error = std::abs(result - expected);
                } else {
                    error = xsf::extended_relative_error(result, expected);
                }

                CAPTURE(b, c, z, result, expected, rtol, atol, error);
                REQUIRE(error <= std::max(rtol * std::abs(expected), atol));
            }
        }
    }
}

TEST_CASE("hyp2f1(nan, nan, nan, complex('nan')) / gh-94", "[hyp2f1][xsf_tests]") {
    double a = std::numeric_limits<double>::quiet_NaN();
    std::complex<double> z = {std::numeric_limits<double>::quiet_NaN(), 0.0};

    std::complex<double> result = xsf::hyp2f1(a, a, a, z);
    REQUIRE(std::isnan(result.real()));
    REQUIRE(std::isnan(result.imag()));
}

TEST_CASE("hyp2f1(a, b, c, x) very negative c - a - b, real path / gh-371", "[hyp2f1][xsf_tests]") {
    /*
    import mpmath

    mpmath.mp.dps = 50

    for a, b, c, x in [(1, 200.5, 101, 0.5), (1, 201, 101, 0.5), (1, 100.25, 51, 0.875)]:
        print(float(mpmath.hyp2f1(a, b, c, x)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double, double>;
    auto [a, b, c, x, desired] = GENERATE(
        test_case{1, 200.5, 101, 0.5, 18.236244010510156}, test_case{1, 201, 101, 0.5, 18.746707942830703},
        test_case{1, 100.25, 51, 0.875, 3.172339232831432e+19}
    );

    const double out = xsf::hyp2f1(a, b, c, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-13;
    CAPTURE(a, b, c, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("hyp2f1(a, b, c, z) very negative c - a - b, complex path / gh-371", "[hyp2f1][xsf_tests]") {
    /*
    import mpmath

    mpmath.mp.dps = 50

    for a, b, c, z in [(1, 200, 101, 0.5), (1, 150, 51, 0.8)]:
        print(complex(mpmath.hyp2f1(a, b, c, z)))
    */
    SET_FP_FORMAT()
    using cd = std::complex<double>;
    using test_case = std::tuple<double, double, double, cd, cd>;
    auto [a, b, c, z, desired] = GENERATE(
        test_case{1, 200, 101, {0.5, 0}, {17.746707942830703, 0}},
        test_case{1, 150, 51, {0.8, 0}, {4.1188640665521755e+34, 0}}
    );

    const cd out = xsf::hyp2f1(a, b, c, z);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-13;
    CAPTURE(a, b, c, z, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("hyp1f1(nan, nan, complex('nan')) / gh-94", "[hyp1f1][xsf_tests]") {
    double a = std::numeric_limits<double>::quiet_NaN();
    std::complex<double> z = {std::numeric_limits<double>::quiet_NaN(), 0.0};

    std::complex<double> result = xsf::hyp1f1(a, a, z);
    REQUIRE(std::isnan(result.real()));
    REQUIRE(std::isnan(result.imag()));
}

TEST_CASE("hyp1f1(a, b, z) large |z| with ph z < 0 after Kummer", "[hyp1f1][xsf_tests]") {
    // Regression test for gh-398: the asymptotic expansion used e^{+pi i a}
    // for Im z < 0 (after Kummer's transformation).
    /*
    import mpmath

    mpmath.mp.dps = 50

    for a, b, z in [(0.5, 1.5, 5 - 25j), (0.5, 1.0, -5 + 25j), (2.5, 4.0, 1 - 40j)]:
        print(complex(mpmath.hyp1f1(a, b, z)))
    */
    SET_FP_FORMAT()
    using cd = std::complex<double>;
    using test_case = std::tuple<double, double, cd, cd>;
    auto [a, b, z, desired] = GENERATE(
        test_case{0.5, 1.5, {5, -25}, {0.24299529873403347, 2.780418451436882}},
        test_case{0.5, 1.0, {-5, 25}, {0.08623955921134369, 0.07114671464042453}},
        test_case{2.5, 4.0, {1, -40}, {0.047946430070203216, -0.001373279980288495}}
    );

    const cd out = xsf::hyp1f1(a, b, z);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-8; // limited by the 8-term asymptotic series
    CAPTURE(a, b, z, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("hyp1f1(a, b, z) tiny a and b", "[hyp1f1][xsf_tests]") {
    // Regression test for gh-454: the series rounded `a + j - 1.0` and `b + j - 1.0`.
    /*
    import mpmath

    mpmath.mp.dps = 50

    for a, b, z in [(1e-17, 2e-17, 1 + 1j), (0.5, 1e-10, 2 + 1j), (-2.0, 1e-17, 1 + 1j)]:
        print(complex(mpmath.hyp1f1(a, b, z)))
    for a, b, x in [(1e-17, 2e-17, 1.0), (0.5, 1e-10, 2.0)]:
        print(float(mpmath.hyp1f1(a, b, x)))
    */
    SET_FP_FORMAT()
    using cd = std::complex<double>;
    using test_case = std::tuple<double, double, cd, cd>;
    auto [a, b, z, desired] = GENERATE(
        test_case{1e-17, 2e-17, {1, 1}, {1.2343469699579426, 1.1436776435894211}},
        test_case{0.5, 1e-10, {2, 1}, {14042459676.44257, 52867238789.94098}},
        test_case{-2.0, 1e-17, {1, 1}, {-2e+17, -2.0}}, test_case{1e-17, 2e-17, {1, 0}, {1.8591409142295225, 0.0}},
        test_case{0.5, 1e-10, {2, 0}, {49777855912.16767, 0.0}}
    );

    const cd out = xsf::hyp1f1(a, b, z);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-14;
    CAPTURE(a, b, z, out, desired, error, rtol);
    REQUIRE(error <= rtol);
    if (z.imag() == 0.0) {
        const double out_real = xsf::hyp1f1(a, b, z.real());
        const double error_real = xsf::extended_relative_error(out_real, desired.real());
        CAPTURE(out_real, error_real);
        REQUIRE(error_real <= rtol);
    }
}
