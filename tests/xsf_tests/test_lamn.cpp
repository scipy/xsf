#include "../testing_utils.h"

#include <xsf/cpu/specfun/specfun.h>

// Regression test for lamn with tiny x: |x| < 1e-100 returned Lambda_k(x) = 0 for k >= 1 and
// Lambda_1'(x) = 0.5, but Lambda_k(0) = 1 and Lambda_k'(x) = -x / (2 (k + 1)) Lambda_{k + 1}(x).

/*
import mpmath

mpmath.mp.dps = 30

x = mpmath.mpf("1e-101")
for k in [0, 1, 3]:
    lam = lambda v: mpmath.gamma(v + 1) * (2 / x) ** v * mpmath.besselj(v, x)
    dl = -x / (2 * (k + 1)) * lam(k + 1)
    print(f"lamn: {k}, {float(lam(k))!r}, {float(dl)!r}")
*/

TEST_CASE("specfun::lamn tiny x", "[lamn][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<int, double, double>;
    auto [k, bl_k, dl_k] =
        GENERATE(test_case{0, 1.0, -5e-102}, test_case{1, 1.0, -2.5e-102}, test_case{3, 1.0, -1.25e-102});

    int nm = 0;
    double bl[4] = {0.0, 0.0, 0.0, 0.0}, dl[4] = {0.0, 0.0, 0.0, 0.0};
    xsf::specfun::lamn(3, 1e-101, &nm, bl, dl);

    const double error = xsf::extended_relative_error(dl[k], dl_k);
    const double rtol = 5e-16;
    CAPTURE(k, bl[k], dl[k], dl_k, error, rtol);
    REQUIRE(nm == 3);
    REQUIRE(bl[k] == bl_k);
    REQUIRE(error <= rtol);
}

// Regression test for lamn with x < -12: the power series for |x| <= 12 was used, which cancels.

/*
import mpmath

mpmath.mp.dps = 30

x = mpmath.mpf(-60)
for k in [0, 1, 2]:
    lam = lambda v: mpmath.gamma(v + 1) * (2 / x) ** v * mpmath.besselj(v, x)
    dl = -x / (2 * (k + 1)) * lam(k + 1)
    print(f"lamn: {k}, {float(lam(k))!r}, {float(dl)!r}")
*/

TEST_CASE("specfun::lamn x < -12", "[lamn][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<int, double, double>;
    auto [k, bl_k, dl_k] = GENERATE(
        test_case{0, -0.09147180408906187, 0.046598383758166315},
        test_case{1, 0.001553279458605544, 0.0031008361182555803},
        test_case{2, 0.00020672240788370536, -8.977047004812257e-05}
    );

    int nm = 0;
    double bl[3] = {0.0, 0.0, 0.0}, dl[3] = {0.0, 0.0, 0.0};
    xsf::specfun::lamn(2, -60.0, &nm, bl, dl);

    const double rtol = 1e-13;
    const double bl_error = xsf::extended_relative_error(bl[k], bl_k);
    const double dl_error = xsf::extended_relative_error(dl[k], dl_k);
    CAPTURE(k, bl[k], bl_k, bl_error, dl[k], dl_k, dl_error, rtol);
    REQUIRE(nm == 2);
    REQUIRE(bl_error <= rtol);
    REQUIRE(dl_error <= rtol);
}

TEST_CASE("specfun::lamn parity", "[lamn][xsf_tests]") {
    // Lambda_n(-x) = Lambda_n(x) and Lambda_n'(-x) = -Lambda_n'(x)
    SET_FP_FORMAT()
    const int n = GENERATE(0, 2, 10);
    const double x = GENERATE(5.0, 20.0, 60.0);

    // lamn writes up to index max(n, 1)
    int nm_pos = 0, nm_neg = 0;
    double bl_pos[16] = {}, dl_pos[16] = {}, bl_neg[16] = {}, dl_neg[16] = {};
    xsf::specfun::lamn(n, x, &nm_pos, bl_pos, dl_pos);
    xsf::specfun::lamn(n, -x, &nm_neg, bl_neg, dl_neg);

    CAPTURE(n, x, nm_pos, nm_neg);
    REQUIRE(nm_neg == nm_pos);
    for (int k = 0; k <= nm_pos; k++) {
        CAPTURE(k, bl_pos[k], bl_neg[k], dl_pos[k], dl_neg[k]);
        REQUIRE(bl_neg[k] == bl_pos[k]);
        REQUIRE(dl_neg[k] == -dl_pos[k]);
    }
}

// Regression test for lamn with large x: the backward recurrence overflowed from about x = 4476, giving nan.
// The second case checks the stored values that are rescaled too.

/*
import mpmath

mpmath.mp.dps = 40

lam = lambda v, x: mpmath.gamma(v + 1) * (2 / x) ** v * mpmath.besselj(v, x)
for n, x, k in [(2, 5000, 2), (1480, 800, 856)]:
    x = mpmath.mpf(x)
    print(f"lamn: {n}, {float(x)!r}, {k}, {float(lam(k, x))!r}, {float(-x / (2 * (k + 1)) * lam(k + 1, x))!r}")
*/

TEST_CASE("specfun::lamn large x", "[lamn][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<int, double, int, double, double>;
    auto [n, x, k, bl_k, dl_k] = GENERATE(
        test_case{2, 5000.0, 2, 2.1265079325321248e-09, -2.919271034712797e-09},
        test_case{1480, 800.0, 856, 1.7983969600283228e-95, -1.2321840028709226e-95}
    );

    int nm = 0;
    std::vector<double> bl(n + 1, 0.0), dl(n + 1, 0.0);
    xsf::specfun::lamn(n, x, &nm, bl.data(), dl.data());

    const double rtol = 1e-12;
    const double bl_error = xsf::extended_relative_error(bl[k], bl_k);
    const double dl_error = xsf::extended_relative_error(dl[k], dl_k);
    CAPTURE(n, x, k, bl[k], bl_k, bl_error, dl[k], dl_k, dl_error, rtol);
    REQUIRE(nm == n);
    REQUIRE(bl_error <= rtol);
    REQUIRE(dl_error <= rtol);
}
