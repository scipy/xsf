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

TEST_CASE("iv large v where u_1(t) vanishes", "[bessel][xsf_tests]") {
    // Regression test for gh-468: x = |v| * sqrt(2/3), where u_1(t) = 0.
    // from mpmath import mp, mpf, besseli
    // mp.dps = 50
    // float(besseli(v, mpf(x)))
    using test_case = std::tuple<double, double, double>;
    auto [v, x, ref] = GENERATE(
        test_case{100, 81.6496580927726, 6392289838.26386}, test_case{-100.5, 82.05790638323647, 7258961449.952213}
    );

    double result = xsf::cyl_bessel_i(v, x);
    double rel_err = xsf::extended_relative_error(result, ref);

    SET_FP_FORMAT()
    CAPTURE(v, x, result, ref, rel_err);
    REQUIRE(rel_err <= 1e-13);
}

TEST_CASE("jv and yv for negative v on overflow", "[bessel][xsf_tests]") {
    // Regression test for gh-383: the first four returned inf, NaN, -inf, and NaN,
    // and the last one must keep using the finite cephes::jv fallback.
    // from mpmath import mp, mpf, besselj, bessely
    // mp.dps = 50
    // bessely(10.25, mpf(1e-100))  # -2.4778e+1033
    // besselj(-10.25, mpf(1e-100))  # 1.7521e+1033
    // bessely(-200.5, 1)  # 3.9331e-437
    // besselj(-1611.5, 200)  # -3.8575e+1246
    // float(besselj(mpf(-10.001), mpf(1e-30)))
    const double inf = std::numeric_limits<double>::infinity();
    CHECK(xsf::cyl_bessel_ye(10.25, 1e-100) == -inf);
    CHECK(xsf::cyl_bessel_je(-10.25, 1e-100) == inf);
    CHECK(xsf::cyl_bessel_y(-200.5, 1.0) == 0);
    CHECK(xsf::cyl_bessel_j(-1611.5, 200.0) == -inf);

    double rel_err = xsf::extended_relative_error(xsf::cyl_bessel_j(-10.001, 1e-30), 3.993385638016181e+305);
    CHECK(rel_err <= 1e-13);
}

TEST_CASE("ive for negative non-integer v and large Im(z)", "[bessel][xsf_tests]") {
    // Regression test for gh-517: the relative error was up to eps * Im(z).
    // from mpmath import mp, mpc, besseli
    // mp.dps = 50
    // besseli(-7.25, mpc(0, 10**9))   # 1.7732342063040383e-06+4.2809660701231e-06j
    // besseli(-0.5, mpc(0, 10**9))    # 1.494894726886609e-05-1.494894726886609e-05j
    // besseli(-7.25, mpc(0, 10**8))   # 3.052711541179276e-05+7.369897604727882e-05j
    using test_case = std::tuple<double, double, std::complex<double>>;
    auto [v, y, ref] = GENERATE(
        test_case{-7.25, 1e9, {1.7732342063040383e-06, 4.2809660701231e-06}},
        test_case{-0.5, 1e9, {1.494894726886609e-05, -1.494894726886609e-05}},
        test_case{-7.25, 1e8, {3.052711541179276e-05, 7.369897604727882e-05}}
    );
    SET_FP_FORMAT();

    std::complex<double> result = xsf::cyl_bessel_ie(v, std::complex<double>(0.0, y));
    double rel_err = xsf::extended_relative_error(result, ref);

    CAPTURE(v, y, result, ref, rel_err);
    REQUIRE(rel_err <= 1e-12);
}
