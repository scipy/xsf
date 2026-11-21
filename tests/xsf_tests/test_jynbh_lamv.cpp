#include "../testing_utils.h"

#include <xsf/cpu/specfun/specfun.h>

// Regression tests for gh-351: (int)(0.9 * x) overflowed for x > 2**31 / 0.9.
// J and Y are ill-conditioned here (relative condition number about x), hence the loose rtol.

/*
import mpmath

mpmath.mp.dps = 30

for x in [2.4e9, 3e9]:
    j0, j1, y0, y1 = [f(n, x) for f in (mpmath.besselj, mpmath.bessely) for n in (0, 1)]
    print(f"jynbh: {x!r}, {float(j0)!r}, {float(j1)!r}, {float(y0)!r}, {float(y1)!r}")
    v = 2.5
    lam = mpmath.gamma(v + 1) * (2 / mpmath.mpf(x)) ** v * mpmath.besselj(v, x)
    print(f"lamv: {x!r}, {float(lam)!r}")
*/

TEST_CASE("specfun::jynbh x > 2**31 / 0.9 gh-351", "[jynbh][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double, double>;
    auto [x, j0, j1, y0, y1] = GENERATE(
        test_case{2.4e9, -4.597762215918059e-06, 1.5624302258418982e-05, 1.5624302259376848e-05, 4.597762219173122e-06},
        test_case{3e9, 8.511574125077861e-06, 1.1822000538658425e-05, 1.182200053723983e-05, -8.511574123107527e-06}
    );

    int nm = 0;
    double bj[2] = {0.0, 0.0}, by[2] = {0.0, 0.0};
    xsf::specfun::jynbh(1, 0, x, &nm, bj, by);

    const double out[4] = {bj[0], bj[1], by[0], by[1]};
    const double desired[4] = {j0, j1, y0, y1};
    const double rtol = 2e-6;
    for (int i = 0; i < 4; ++i) {
        const double error = xsf::extended_relative_error(out[i], desired[i]);
        CAPTURE(x, i, out[i], desired[i], error, rtol);
        REQUIRE(error <= rtol);
    }
}

TEST_CASE("specfun::lamv x > 2**31 / 0.9 gh-351", "[lamv][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(test_case{2.4e9, -5.194543092337725e-28}, test_case{3e9, -5.483360479520362e-28});

    double vm = 0.0, vl[3] = {0.0, 0.0, 0.0}, dl[3] = {0.0, 0.0, 0.0};
    xsf::specfun::lamv(2.5, x, &vm, vl, dl);

    const double error = xsf::extended_relative_error(vl[2], desired);
    const double rtol = 2e-6;
    CAPTURE(x, vl[2], desired, error, rtol);
    REQUIRE(error <= rtol);
}

// Regression test for lamv with x < 0: Lambda_v is even, so the derivative must be odd.

/*
import mpmath

mpmath.mp.dps = 30

v = mpmath.mpf("1.5")
lam = lambda x: mpmath.gamma(v + 1) * (2 / x) ** v * mpmath.besselj(v, x)
for x in [-1, -20]:
    x = mpmath.mpf(x)
    # Lambda_v is even in x, so Lambda_v' is odd
    print(f"lamv: {float(x)!r}, {float(lam(-x))!r}, {float(-mpmath.diff(lam, -x))!r}")
*/

TEST_CASE("specfun::lamv x < 0", "[lamv][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, vl1, dl1] = GENERATE(
        test_case{-1.0, 0.9035060368192703, 0.1861051560341216},
        test_case{-20.0, -0.0027182609945775797, -0.007254828529643844}
    );

    double vm = 0.0, vl[2] = {0.0, 0.0}, dl[2] = {0.0, 0.0};
    xsf::specfun::lamv(1.5, x, &vm, vl, dl);

    const double rtol = 1e-13;
    const double vl_error = xsf::extended_relative_error(vl[1], vl1);
    const double dl_error = xsf::extended_relative_error(dl[1], dl1);
    CAPTURE(x, vl[1], vl1, vl_error, dl[1], dl1, dl_error, rtol);
    REQUIRE(vl_error <= rtol);
    REQUIRE(dl_error <= rtol);
}

TEST_CASE("specfun::lamv parity", "[lamv][xsf_tests]") {
    // Lambda_v(-x) = Lambda_v(x) and Lambda_v'(-x) = -Lambda_v'(x)
    SET_FP_FORMAT()
    const double v = GENERATE(0.5, 1.5, 4.2, 10.7);
    const double x = GENERATE(0.5, 7.0, 20.0, 50.0);

    // lamv writes up to index max((int)v, 1)
    double vm_pos = 0.0, vm_neg = 0.0;
    double vl_pos[16] = {}, dl_pos[16] = {}, vl_neg[16] = {}, dl_neg[16] = {};
    xsf::specfun::lamv(v, x, &vm_pos, vl_pos, dl_pos);
    xsf::specfun::lamv(v, -x, &vm_neg, vl_neg, dl_neg);

    CAPTURE(v, x, vm_pos, vm_neg);
    REQUIRE(vm_neg == vm_pos);
    for (int k = 0; k <= static_cast<int>(vm_pos); k++) {
        CAPTURE(k, vl_pos[k], vl_neg[k], dl_pos[k], dl_neg[k]);
        REQUIRE(vl_neg[k] == vl_pos[k]);
        REQUIRE(dl_neg[k] == -dl_pos[k]);
    }
}

// Regression test for lamv with large x: the backward recurrence overflowed from about x = 4476, giving nan.
// The second case checks the stored values that are rescaled too.

/*
import mpmath

mpmath.mp.dps = 40

lam = lambda v, x: mpmath.gamma(v + 1) * (2 / x) ** v * mpmath.besselj(v, x)
for v, x, k in [(4600.5, 5000, 2), (1480.5, 800, 856)]:
    x, o = mpmath.mpf(x), mpmath.mpf(k) + mpmath.mpf("0.5")
    print(f"lamv: {v!r}, {float(x)!r}, {k}, {float(lam(o, x))!r}, {float(-x / (2 * (o + 1)) * lam(o + 1, x))!r}")
*/

TEST_CASE("specfun::lamv large x", "[lamv][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, int, double, double>;
    auto [v, x, k, vl_k, dl_k] = GENERATE(
        test_case{4600.5, 5000.0, 2, 1.1854482230005149e-10, -1.870246475852011e-11},
        test_case{1480.5, 800.0, 856, 2.1794710709981976e-95, -1.4909006138395703e-95}
    );

    double vm = 0.0;
    std::vector<double> vl((int)v + 2, 0.0), dl((int)v + 2, 0.0);
    xsf::specfun::lamv(v, x, &vm, vl.data(), dl.data());

    const double rtol = 1e-12;
    const double vl_error = xsf::extended_relative_error(vl[k], vl_k);
    const double dl_error = xsf::extended_relative_error(dl[k], dl_k);
    CAPTURE(v, x, k, vl[k], vl_k, vl_error, dl[k], dl_k, dl_error, rtol);
    REQUIRE(vl_error <= rtol);
    REQUIRE(dl_error <= rtol);
}
