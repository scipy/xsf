/* C++ translation of tests from scipy/special/tests/test_boxcox.py
 * https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_boxcox.py
 */

#include "../testing_utils.h"

#include <xsf/boxcox.h>

TEST_CASE("boxcox basic", "[boxcox][xsf_tests]") {
    constexpr double atol = 1.5e-7;
    constexpr double rtol = 0.0;
    std::vector<double> xs = {0.5, 1.0, 2.0, 4.0};

    for (const auto x : xs) {
        // lambda = 0  =>  y = log(x)
        auto output = xsf::boxcox(x, 0.0);
        auto expected = std::log(x);
        auto abs_error = std::abs(output - expected);
        auto tol = atol + rtol * std::abs(expected);
        CAPTURE(x, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);

        // lambda = 1  =>  y = x - 1
        output = xsf::boxcox(x, 1.0);
        expected = x - 1.0;
        abs_error = std::abs(output - expected);
        tol = atol + rtol * std::abs(expected);
        CAPTURE(x, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);

        // lambda = 2  =>  y = 0.5*(x**2 - 1)
        output = xsf::boxcox(x, 2.0);
        expected = 0.5 * (x * x - 1.0);
        abs_error = std::abs(output - expected);
        tol = atol + rtol * std::abs(expected);
        CAPTURE(x, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);
    }

    // x = 0 and lambda > 0  =>  y = -1 / lambda
    for (const auto lmbda : {0.5, 1.0, 2.0}) {
        const auto output = xsf::boxcox(0.0, lmbda);
        const auto expected = -1.0 / lmbda;
        const auto abs_error = std::abs(output - expected);
        const auto tol = atol + rtol * std::abs(expected);
        CAPTURE(lmbda, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);
    }
}

TEST_CASE("boxcox underflow", "[boxcox][xsf_tests]") {
    constexpr double x = 1.0 + 1e-15;
    constexpr double lmbda = 1e-306;
    constexpr double rtol = 1e-14;
    const auto output = xsf::boxcox(x, lmbda);
    const auto expected = std::log(x);
    const auto rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(x, lmbda, output, expected, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("boxcox nonfinite", "[boxcox][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    auto [x, lmbda] = GENERATE(test_case{-1.0, 0.5}, test_case{-1.0, 2.0}, test_case{-0.5, -1.5});
    CAPTURE(x, lmbda);
    REQUIRE(std::isnan(xsf::boxcox(x, lmbda)));

    REQUIRE(xsf::boxcox(0.0, -2.5) == -std::numeric_limits<double>::infinity());
    REQUIRE(xsf::boxcox(0.0, 0.0) == -std::numeric_limits<double>::infinity());
}

TEST_CASE("boxcox1p basic", "[boxcox][xsf_tests]") {
    constexpr double atol = 1.5e-7;
    constexpr double rtol = 0.0;
    std::vector<double> xs = {-0.25, -1e-20, 0.0, 1e-20, 0.25, 1.0, 3.0};

    for (const auto x : xs) {
        // lambda = 0  =>  y = log(1+x)
        auto output = xsf::boxcox1p(x, 0.0);
        auto expected = std::log1p(x);
        auto abs_error = std::abs(output - expected);
        auto tol = atol + rtol * std::abs(expected);
        CAPTURE(x, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);

        // lambda = 1  =>  y = x
        output = xsf::boxcox1p(x, 1.0);
        expected = x;
        abs_error = std::abs(output - expected);
        tol = atol + rtol * std::abs(expected);
        CAPTURE(x, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);

        // lambda = 2  =>  y = 0.5*((1+x)**2 - 1) = 0.5*x*(2 + x)
        output = xsf::boxcox1p(x, 2.0);
        expected = 0.5 * x * (2.0 + x);
        abs_error = std::abs(output - expected);
        tol = atol + rtol * std::abs(expected);
        CAPTURE(x, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);
    }

    // x = -1 and lambda > 0  =>  y = -1 / lambda
    for (const auto lmbda : {0.5, 1.0, 2.0}) {
        const auto output = xsf::boxcox1p(-1.0, lmbda);
        const auto expected = -1.0 / lmbda;
        const auto abs_error = std::abs(output - expected);
        const auto tol = atol + rtol * std::abs(expected);
        CAPTURE(lmbda, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);
    }
}

TEST_CASE("boxcox1p underflow", "[boxcox][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    auto [x, lmbda] = GENERATE(test_case{1e-15, 1e-306}, test_case{1e-306, 1e-18});
    const double rtol = 1e-14;
    const auto output = xsf::boxcox1p(x, lmbda);
    const auto expected = std::log1p(x);
    const auto rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(x, lmbda, output, expected, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("boxcox1p nonfinite", "[boxcox][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    auto [x, lmbda] = GENERATE(test_case{-2.0, 0.5}, test_case{-2.0, 2.0}, test_case{-1.5, -1.5});
    CAPTURE(x, lmbda);
    REQUIRE(std::isnan(xsf::boxcox1p(x, lmbda)));

    REQUIRE(xsf::boxcox1p(-1.0, -2.5) == -std::numeric_limits<double>::infinity());
    REQUIRE(xsf::boxcox1p(-1.0, 0.0) == -std::numeric_limits<double>::infinity());
}

TEST_CASE("boxcox inverse", "[boxcox][xsf_tests]") {
    constexpr double atol = 1.5e-7;
    constexpr double rtol = 0.0;
    std::vector<double> xs = {0.0, 1.0, 2.0};
    std::vector<double> lmbdas = {0.0, 1.0, 2.0};

    for (std::size_t i = 0; i < xs.size(); ++i) {
        const double x = xs[i];
        const double lmbda = lmbdas[i];
        const double y = xsf::boxcox(x, lmbda);
        auto output = xsf::inv_boxcox(y, lmbda);
        auto expected = x;
        auto abs_error = std::abs(output - expected);
        auto tol = atol + rtol * std::abs(expected);
        CAPTURE(x, lmbda, y, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);

        const double y1p = xsf::boxcox1p(x, lmbda);
        output = xsf::inv_boxcox1p(y1p, lmbda);
        expected = x;
        abs_error = std::abs(output - expected);
        tol = atol + rtol * std::abs(expected);
        CAPTURE(x, lmbda, y1p, output, expected, atol, rtol, abs_error, tol);
        REQUIRE(abs_error <= tol);
    }
}

TEST_CASE("inv_boxcox1p underflow", "[boxcox][xsf_tests]") {
    constexpr double x = 1e-15;
    constexpr double lmbda = 1e-306;
    constexpr double rtol = 1e-14;
    const auto output = xsf::inv_boxcox1p(x, lmbda);
    const auto rel_error = xsf::extended_relative_error(output, x);
    CAPTURE(x, lmbda, output, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("inv_boxcox premature overflow", "[boxcox][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    auto [x, lmbda] = GENERATE(test_case{100.0, 155.0}, test_case{0.01, -155.0});

    const double y = xsf::boxcox(x, lmbda);
    CAPTURE(x, lmbda, y);
    REQUIRE(std::isfinite(y));
    constexpr double rtol = 1e-14;
    auto output = xsf::inv_boxcox(y, lmbda);
    auto rel_error = xsf::extended_relative_error(output, x);
    CAPTURE(x, lmbda, y, output, rtol, rel_error);
    REQUIRE(rel_error <= rtol);

    const double y1p = xsf::boxcox1p(x - 1.0, lmbda);
    CAPTURE(y1p);
    REQUIRE(std::isfinite(y1p));
    const auto expected = x - 1.0;
    output = xsf::inv_boxcox1p(y1p, lmbda);
    rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(x, lmbda, y1p, output, expected, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
