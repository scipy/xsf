#include "../testing_utils.h"

#include <xsf/log_exp.h>

#include <cmath>
#include <limits>
#include <tuple>

TEST_CASE("logdiffexp accuracy", "[logdiffexp][xsf_tests]") {
    SET_FP_FORMAT()

    // References rounded to double from 80-digit mpmath evaluations
    // using the exact binary floating-point inputs.
    using test_case = std::tuple<double, double, double>;
    auto [a, b, desired] = GENERATE(
        test_case{1000.0, 999.0, 999.5413248546129}, test_case{1000.0, 999.999999999, 979.276723587097},
        test_case{0.0, -1e-12, -27.631021115929048}, test_case{0.0, -50.0, -1.9287498479639178e-22}
    );

    const double result = xsf::logdiffexp(a, b);
    const double error = xsf::extended_relative_error(result, desired);
    CAPTURE(a, b, result, desired, error);
    REQUIRE(error <= 2e-15);
}

TEST_CASE("logdiffexp special values", "[logdiffexp][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    REQUIRE(xsf::logdiffexp(1.0, 1.0) == -inf);
    REQUIRE(xsf::logdiffexp(-inf, -inf) == -inf);
    REQUIRE(xsf::logdiffexp(1.0, -inf) == 1.0);
    REQUIRE(xsf::logdiffexp(inf, 1.0) == inf);
    REQUIRE(xsf::logdiffexp(inf, -inf) == inf);

    REQUIRE(std::isnan(xsf::logdiffexp(inf, inf)));
    REQUIRE(std::isnan(xsf::logdiffexp(1.0, 2.0)));
    REQUIRE(std::isnan(xsf::logdiffexp(-inf, 1.0)));
    REQUIRE(std::isnan(xsf::logdiffexp(nan, 1.0)));
    REQUIRE(std::isnan(xsf::logdiffexp(1.0, nan)));
    REQUIRE(std::isnan(xsf::logdiffexp(nan, nan)));
}

TEST_CASE("logdiffexp accuracy near a zero result", "[logdiffexp][xsf_tests]") {
    SET_FP_FORMAT()

    using test_case = std::tuple<double, double>;
    auto [a, desired] = GENERATE(
        test_case{0.6931471805599452, -2.6842554120195732e-16}, test_case{0.6931471805599453, -4.6380936276925993e-17},
        test_case{0.6931471805599454, 1.7566366864810531e-16}
    );

    const double result = xsf::logdiffexp(a, 0.0);
    const double error = xsf::extended_relative_error(result, desired);
    CAPTURE(a, result, desired, error);
    REQUIRE(error <= 2e-15);
}

TEST_CASE("logdiffexp subnormal results", "[logdiffexp][xsf_tests]") {
    SET_FP_FORMAT()

    using test_case = std::tuple<double, double, double>;
    auto [a, b, desired] = GENERATE(
        test_case{0x1.9788af0d610e5p-1012, -701.0, -0x0.00000000002e9p-1022},
        test_case{0x1.9788af0d610e6p-1012, -701.0, 0x0.0000000000117p-1022},
        test_case{0x1.9788af0d610e7p-1012, -701.0, 0x0.0000000000517p-1022},
        test_case{0x1.78694fe9f73cbp-1009, -699.0, -0x0.0000000002829p-1022},
        test_case{0x1.78694fe9f73ccp-1009, -699.0, -0x0.0000000000829p-1022},
        test_case{0x1.78694fe9f73cdp-1009, -699.0, 0x0.00000000017d7p-1022}
    );

    const double result = xsf::logdiffexp(a, b);
    CAPTURE(a, b, result, desired);

    // Allow one subnormal spacing instead of a relative tolerance.
    REQUIRE(std::abs(result - desired) <= std::numeric_limits<double>::denorm_min());
}

TEST_CASE("logdiffexp float accuracy", "[logdiffexp][xsf_tests]") {
    SET_FP_FORMAT()

    using test_case = std::tuple<float, float, float>;
    auto [a, b, desired] = GENERATE(
        test_case{1000.0f, 999.0f, 0x1.f3c54ap+9f}, test_case{0.0f, -0x1.197998p-40f, -0x1.ba18aap+4f},
        test_case{0.0f, -50.0f, -0x1.d257d6p-73f}, test_case{0x1.62e430p-1f, 0.0f, 0x1.05c610p-28f}
    );

    const float result = xsf::logdiffexp(a, b);
    const float error = std::abs((result - desired) / desired);
    CAPTURE(a, b, result, desired, error);
    REQUIRE(error <= 4 * std::numeric_limits<float>::epsilon());
}
