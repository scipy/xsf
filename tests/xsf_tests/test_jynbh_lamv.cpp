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
