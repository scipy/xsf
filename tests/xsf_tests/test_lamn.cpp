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
