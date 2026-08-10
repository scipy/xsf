#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/tetragamma.h>

// Apery's constant zeta(3), used throughout via the A&S 6.4.4 closed form
// psi''(1) = -2*zeta(3) and the A&S 6.4.6 recurrence psi''(z+1) = psi''(z) + 2/z^3.
static constexpr double ZETA3 = 1.2020569031595942854;

TEST_CASE("tetragamma real", "[tetragamma][xsf_tests]") {
    double rtol = 1e-14;
    using test_case = std::tuple<double, double>;

    auto [z, ref_tetragamma] = GENERATE(
        // Closed-form reference values, A&S 6.4.4 (psi''(1) = -2*zeta(3)) and
        // A&S 6.4.6 (psi''(1/2) = -14*zeta(3)), extended by the integer recurrence
        // psi''(n+1) = psi''(n) + 2/n^3.
        test_case{0.5, -14.0 * ZETA3},
        test_case{1.0, -2.0 * ZETA3},
        test_case{2.0, -2.0 * ZETA3 + 2.0},
        test_case{3.0, -2.0 * ZETA3 + 9.0 / 4.0},
        test_case{4.0, -2.0 * ZETA3 + 251.0 / 108.0},
        test_case{5.0, -2.0 * ZETA3 + 2035.0 / 864.0},
        test_case{10.0, -2.0 * ZETA3 + 2.0 * (19148110939.0 / 16003008000.0)},
        // Reference values from test_trigamma computed using mpmath
        test_case{-7.5, -0.015564511795903923324}, test_case{-3.2, 246.9761555554604978}
    );

    double result_tetragamma = xsf::tetragamma(z);
    double rel_err_tetragamma = xsf::extended_relative_error(result_tetragamma, ref_tetragamma);

    CAPTURE(z, result_tetragamma, ref_tetragamma, rel_err_tetragamma, rtol);
    REQUIRE(rel_err_tetragamma <= rtol);
}

TEST_CASE("tetragamma real poles", "[tetragamma][xsf_tests]") {
    auto z = GENERATE(-1.0, 0.0);
    double result_tetragamma = xsf::tetragamma(z);
    CAPTURE(z, result_tetragamma);
    REQUIRE(std::isinf(result_tetragamma));
}

TEST_CASE("tetragamma complex", "[tetragamma][xsf_tests]") {
    double rtol = 1e-14;
    using test_case = std::tuple<std::complex<double>, std::complex<double>>;

    auto [z, ref_tetragamma] = GENERATE(
        // Test points used in test_trigamma, taken from
        // https://github.com/JuliaMath/SpecialFunctions.jl/blob/c9ff36085c9008b9bdaa457e745f459a1827e52e/test/gamma.jl#L261-L294
        // computed using mpmath (mp.dps=60).
        test_case{
            {8.0, 0.0},
            {-0.0176995691957677739092916777362138785173389598024749258780319, 0.0}
        },
        test_case{
            {0.0, 8.0},
            {0.0155022836781734350054183962027714397271268582016030387460863,
             -0.00195312500000000001834374260890084586707709090054569497956919}
        },
        test_case{
            {-3.2, 0.1},
            {29.285628461389160743289806565288260993728667020276678436109,
             177.776489607267829881651358046772576380967686212268333777227}
        },
        test_case{
            {-10.0, 1.0e100},
            {1.0e-200, -2.0999999999999998e-299}
        }
    );

    std::complex<double> result_tetragamma = xsf::tetragamma(z);
    double rel_err_tetragamma = xsf::extended_relative_error(result_tetragamma, ref_tetragamma);

    CAPTURE(z, result_tetragamma, ref_tetragamma, rel_err_tetragamma, rtol);
    REQUIRE(rel_err_tetragamma <= rtol);
}

TEST_CASE("tetragamma complex poles", "[tetragamma][xsf_tests]") {
    auto z = GENERATE(std::complex<double>{0.0, 0.0}, std::complex<double>{-1.0, 0.0},
                       std::complex<double>{-2.0, 0.0});
    std::complex<double> result_tetragamma = xsf::tetragamma(z);
    CAPTURE(z, result_tetragamma);
    REQUIRE(std::isinf(result_tetragamma.real()));
    REQUIRE(result_tetragamma.imag() == 0.0);
}
