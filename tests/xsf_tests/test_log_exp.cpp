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

TEST_CASE("logdiffexp branch boundaries", "[logdiffexp][xsf_tests]") {
    SET_FP_FORMAT()

    // references rounded from 500-digit mpmath evaluations of exact inputs
    using test_case = std::tuple<double, double, double>;
    auto [a, b, desired] = GENERATE(
        // neighboring inputs around the Taylor-series boundary at a = 0.5
        test_case{0x1.fffffffffffffp-2, -0x1.bb235fcf0fcccp-2, -0x1.cbd13a873a7e2p-54},
        test_case{0x1.0000000000000p-1, -0x1.bb235fcf0fcccp-2, -0x1.2df50d2ad0a2ep-57},
        test_case{0x1.0000000000001p-1, -0x1.bb235fcf0fcccp-2, 0x1.9333480f335f9p-53},

        // inputs around the positive fallback threshold for a = 1
        test_case{1.0, 0x1.137280c3c0cbfp-4, 0x1.0000000000001p-1},
        test_case{1.0, 0x1.137280c3c0ccfp-4, 0x1.0000000000000p-1},
        test_case{1.0, 0x1.137280c3c0cdfp-4, 0x1.ffffffffffffdp-2},

        // inputs around the negative fallback threshold for a = 1
        test_case{1.0, 0x1.7eba9e8033f95p-1, -0x1.ffffffffffffcp-2},
        test_case{1.0, 0x1.7eba9e8033f96p-1, -0x1.0000000000001p-1},
        test_case{1.0, 0x1.7eba9e8033f97p-1, -0x1.0000000000005p-1},

        // inputs around the positive fallback threshold for a < 1
        test_case{0.25, -0x1.e42b9ea8ca691p+0, 0x1.0000000000001p-3},
        test_case{0.25, -0x1.e42b9ea8ca690p+0, 0x1.fffffffffffffp-4},
        test_case{0.25, -0x1.e42b9ea8ca68fp+0, 0x1.ffffffffffffdp-4},

        // inputs around the negative fallback threshold for a < 1
        test_case{0.25, -0x1.d330265c61456p-1, -0x1.ffffffffffffdp-4},
        test_case{0.25, -0x1.d330265c61455p-1, -0x1.0000000000000p-3},
        test_case{0.25, -0x1.d330265c61454p-1, -0x1.0000000000002p-3}
    );

    const double result = xsf::logdiffexp(a, b);
    const double error = xsf::extended_relative_error(result, desired);
    CAPTURE(a, b, result, desired, error);
    REQUIRE(error <= 2e-15);
}

TEST_CASE("logdiffexp scaled-tail boundary", "[logdiffexp][xsf_tests]") {
    SET_FP_FORMAT()

    // neighboring inputs around b - a = -700
    using test_case = std::tuple<double, double, double>;
    auto [a, b, desired] = GENERATE(
        test_case{0x1.14f2b0fb92e55p-1010, -0x1.5e00000000001p+9, -0x0.0000000000256p-1022},
        test_case{0x1.14f2b0fb9307fp-1010, -0x1.5e00000000000p+9, -0x0.00000000000acp-1022},
        test_case{0x1.14f2b0fb932a9p-1010, -0x1.5dfffffffffffp+9, 0x0.00000000000fep-1022}
    );

    const double result = xsf::logdiffexp(a, b);
    CAPTURE(a, b, result, desired);
    REQUIRE(std::abs(result - desired) <= std::numeric_limits<double>::denorm_min());
}
