#include "../testing_utils.h"
#include <tuple>
#include <xsf/cpu/bessel.h>

TEST_CASE("iv tiny inputs", "[bessel][xsf_tests]") {
    using test_case = std::tuple<long, double, double, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref, rtol] = GENERATE(
        test_case{1, 1e-50, 5e-51, 1e-14}, test_case{100, 0.1, 8.45293498689205e-289, 1e-10},
        test_case{10, 5.2250558491838786e-27, 4.0817497902418315e-273, 1e-14},
        test_case{3, 1e-50, 2.0833333333333333e-152, 1e-14}, test_case{10, 1e-29, 2.6911444554673707e-300, 1e-8},
        test_case{20, 1e-14, 3.919904349624791e-305, 1e-8}, test_case{200, 5, 5.065429051896385e-296, 1e-10}
    );

    double result = xsf::cyl_bessel_i(n, z); // calls cephes::iv
    double rel_err = xsf::extended_relative_error(result, ref);

    CAPTURE(n, z, result, ref, rel_err, rtol);
    REQUIRE(rel_err <= rtol);
}

TEST_CASE("jv power series overflow keeps the sign", "[bessel][xsf_tests]") {
    // Regression test for gh-341: the first case returned +inf.
    // mpmath.besselj(v, x) is -4.2372e+431 and 1.6148e+428, respectively.
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE(xsf::cyl_bessel_j(-123.57388919180612, 0.028605192190944706) == -inf);
    REQUIRE(xsf::cyl_bessel_j(-124.68311766689978, 0.032989216721374115) == inf);
}

TEST_CASE("real jv and jve retain tiny normal values", "[bessel][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    // Reference values computed with mpmath at 120 digits.
    auto [v, x, expected] = GENERATE(
        test_case{5.0, 1e-60, 2.604166666666666e-304}, test_case{2.5, 1e-122, 5.319230405352436e-307},
        test_case{5.0, -1e-60, -2.604166666666666e-304}, test_case{-5.0, 1e-60, -2.604166666666666e-304},
        test_case{-5.0, -1e-60, 2.604166666666666e-304}, test_case{6.0, -1e-50, 2.170138888888889e-305},
        test_case{-6.0, 1e-50, 2.170138888888889e-305},
        test_case{0.1, std::numeric_limits<double>::denorm_min(), 4.580731034319154e-33},
        test_case{0.5, std::numeric_limits<double>::denorm_min(), 1.7735048886036273e-162}
    );
    const double j = xsf::cyl_bessel_j(v, x);
    const double je = xsf::cyl_bessel_je(v, x);
    CAPTURE(v, x, j, je, expected);
    CHECK(xsf::extended_relative_error(j, expected) <= 2e-13);
    CHECK(xsf::extended_relative_error(je, expected) <= 2e-13);
}

TEST_CASE("real jv and jve zero arguments and underflow", "[bessel][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    auto [v, x, expected] = GENERATE(test_case{0.0, 0.0, 1.0}, test_case{5.0, 0.0, 0.0}, test_case{5.0, 1e-100, 0.0});
    CHECK(xsf::cyl_bessel_j(v, x) == expected);
    CHECK(xsf::cyl_bessel_je(v, x) == expected);
}
