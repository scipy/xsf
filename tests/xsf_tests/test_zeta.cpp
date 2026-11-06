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

TEST_CASE("Hurwitz zeta at negative shifts and large orders", "[zeta][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    const double inf = std::numeric_limits<double>::infinity();
    // Values from 150-digit finite sums and positive-argument Hurwitz zeta tails.
    auto [x, q, desired] = GENERATE_COPY(
        test_case{2.0, -0.75, 18.975106932284888}, test_case{3.0, -0.75, 62.29349959839809},
        test_case{31.0, -7.5, 1.5931355392154464e-29}, test_case{2.0, -10000000000.5, 9.869604400989358},
        test_case{3.0, -10000000000.5, 4.999999999e-21}, test_case{31.0, -0.49999999999999994, -1.1305221409216952e-05},
        test_case{1031.0, -0.49999999999999994, -5.267741079373558e+297},
        test_case{31.0, -0.5000000000000001, 3.304063430367369e-05}, test_case{3.0, -1e-15, -9.999999999999998e+44},
        test_case{1031.0, -0.5, 2.8178113466323453e-182}, test_case{1e300, 1.5, 0.0}, test_case{inf, 1.5, 0.0},
        test_case{inf, 0.5, inf}, test_case{inf, 1.0, 1.0}
    );
    const double out = xsf::zeta(x, q);
    CAPTURE(x, q, out, desired);
    REQUIRE((out == desired || xsf::extended_relative_error(out, desired) < 2e-13));
}
