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

TEST_CASE("real iv and ive retain tiny normal values", "[bessel][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double>;
    // Reference values computed with mpmath at 120 digits.
    auto [v, x, expected_i, expected_ie] = GENERATE(
        test_case{46.5, 1e-5, 8.445815628803095e-306, 8.445731171069096e-306},
        test_case{49.5, 2.5e-5, 4.618931182184556e-307, 4.618815710348406e-307},
        test_case{47.0, -1.2e-5, -1.4469089975685643e-305, -1.4468916347647704e-305},
        test_case{-47.0, 1.2e-5, 1.4469089975685643e-305, 1.4468916347647704e-305},
        test_case{0.1, std::numeric_limits<double>::denorm_min(), 4.580731034319154e-33, 4.580731034319154e-33},
        test_case{0.5, std::numeric_limits<double>::denorm_min(), 1.7735048886036274e-162, 1.7735048886036274e-162}
    );
    const double i = xsf::cyl_bessel_i(v, x);
    const double ie = xsf::cyl_bessel_ie(v, x);
    CAPTURE(v, x, i, ie, expected_i, expected_ie);
    CHECK(xsf::extended_relative_error(i, expected_i) <= 2e-13);
    CHECK(xsf::extended_relative_error(ie, expected_ie) <= 2e-13);
}

TEST_CASE("iv finite results near negative integer orders", "[bessel][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    auto [v, x, expected] = GENERATE(
        test_case{-48.000001, 1e-5, 7.279714030645103e+307}, test_case{-47.999999, 1e-5, -7.279480110682128e+307},
        test_case{-6.999999999999544, 9.291987883093336e-45, 7.021010206327728e+300},
        test_case{-0.9999999999999831, 4.77905e-319, 7.0622397754471685e+304},
        test_case{-50.00000000000041, 1.9561109286550917e-5, 7.602367782332035e+300}
    );
    const double actual = xsf::cyl_bessel_i(v, x);
    CAPTURE(v, x, actual, expected);
    CHECK(xsf::extended_relative_error(actual, expected) <= 5e-13);
}

TEST_CASE("iv overflow keeps the sign at large negative order", "[bessel][xsf_tests]") {
    const double actual = xsf::cyl_bessel_i(-2205871995.2311125, 0.21626314137829938);
    CHECK(actual == -std::numeric_limits<double>::infinity());
}

TEST_CASE("jv power series overflow keeps the sign", "[bessel][xsf_tests]") {
    // Regression test for gh-341: the first case returned +inf.
    // mpmath.besselj(v, x) is -4.2372e+431 and 1.6148e+428, respectively.
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE(xsf::cyl_bessel_j(-123.57388919180612, 0.028605192190944706) == -inf);
    REQUIRE(xsf::cyl_bessel_j(-124.68311766689978, 0.032989216721374115) == inf);
}
