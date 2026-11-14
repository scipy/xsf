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
