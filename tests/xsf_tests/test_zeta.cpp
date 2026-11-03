#include "../testing_utils.h"
#include <xsf/zeta.h>

TEST_CASE("zeta(x, q) for x near 1 and tiny q", "[zeta][xsf_tests]") {
    // Regression test for gh-345: the direct sum stopped early, dropping the
    // tail of about (q + n)^(1 - x) / (x - 1).
    // Reference values: float(mpmath.zeta(x, q)) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, q, desired] = GENERATE(
        test_case{1.0000000000000002, 1e-15, 5503599627370504.0},
        test_case{1.000000000001, 1e-17, 1.0000099991502206e+17}
    );

    const double out = xsf::zeta(x, q);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-16;
    CAPTURE(x, q, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("zeta(x, q=1) == riemann_zeta(x) for float, double and complex", "[zeta][xsf_tests]") {
    SECTION("double") {
        for (double x : linspace(-20.0, 20.0, 81)) {
            REQUIRE(xsf::zeta(x, 1.0) == xsf::riemann_zeta(x));
        }
    }

    SECTION("float") {
        for (float x : linspace(-20.0f, 20.0f, 81)) {
            REQUIRE(xsf::zeta(x, 1.0f) == xsf::riemann_zeta(x));
        }
    }

    SECTION("complex") {
        using std::complex;
        for (double re : linspace(-5.0, 5.0, 21)) {
            for (double im : linspace(-10.0, 10.0, 41)) {
                complex<double> z(re, im);
                REQUIRE(xsf::zeta(z, 1.0) == xsf::riemann_zeta(z));
            }
        }
    }
}

TEST_CASE("zeta(x, q) for q > 1e8", "[zeta][xsf_tests]") {
    // Regression test for gh-345: the asymptotic expansion omitted the B_2 term.
    // Reference values: float(mpmath.zeta(x, q)) with mpmath.mp.dps = 50.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, q, desired] = GENERATE(
        test_case{20.0, 1.0000001e8, 5.263148394745859e-154}, test_case{30.0, 1.0000001e8, 3.4482663620824903e-234},
        test_case{35.0, 1.2e8, 5.975287124820555e-277}
    );

    const double out = xsf::zeta(x, q);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-16;
    CAPTURE(x, q, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
