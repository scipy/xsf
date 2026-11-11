#include "../testing_utils.h"
#include <complex>
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

TEST_CASE(
    "cyl_hankel_1e, cyl_hankel_2e and cyl_bessel_ke near the axes",
    "[cyl_hankel_1e][cyl_hankel_2e][cyl_bessel_ke][xsf_tests]"
) {
    // Regression test for gh-516: the phase exp(-2i*Im(z)) was applied with an error of ~eps*abs(z).
    /*
    import mpmath as mp

    mp.mp.dps = 40
    h1e = mp.hankel1(0, mp.mpc(1e15, -1)) * mp.exp(-1j * mp.mpc(1e15, -1))
    h2e = mp.hankel2(0.25, mp.mpc(-1e15, 0)) * mp.exp(1j * mp.mpc(-1e15, 0))
    kve = mp.besselk(0, mp.mpc(-1, 1e9)) * mp.exp(mp.mpc(-1, 1e9))
    */
    using std::complex;
    const double rtol = 1e-13;
    {
        const complex<double> z{1e15, -1.0}, ref{1.784124116152772e-08, -1.7841241161527706e-08};
        const complex<double> w = xsf::cyl_hankel_1e(0.0, z);
        const double err = xsf::extended_relative_error(w, ref);
        CAPTURE(z, w, ref, err);
        CHECK(err <= rtol);
    }
    {
        const complex<double> z{-1e15, 0.0}, ref{-1.2192465936490278e-08, -1.3228340439591988e-08};
        const complex<double> w = xsf::cyl_hankel_2e(0.25, z);
        const double err = xsf::extended_relative_error(w, ref);
        CAPTURE(z, w, ref, err);
        CHECK(err <= rtol);
    }
    {
        const complex<double> z{-1.0, 1e9}, ref{2.8024956071480284e-05, -2.8024956092499003e-05};
        const complex<double> w = xsf::cyl_bessel_ke(0.0, z);
        const double err = xsf::extended_relative_error(w, ref);
        CAPTURE(z, w, ref, err);
        CHECK(err <= rtol);
    }
}
