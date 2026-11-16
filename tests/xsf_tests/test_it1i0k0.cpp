#include "../testing_utils.h"

#include <tuple>
#include <xsf/cpu/bessel.h>

TEST_CASE("it1i0k0 gh-379", "[it1i0k0][xsf_tests]") {
    // Regression test for gh-379: the K0 integral lost up to 2e-11 near x = 12, and the I0
    // integral up to 1.3e-7 at x = 20.

    /*
    import mpmath as mp

    mp.mp.dps = 40

    for x in [12.0, 18.91, 20.0, 45.0]:
        ti = x * mp.hyp1f2(0.5, 1, 1.5, mp.mpf(x) ** 2 / 4)
        tk = mp.pi / 2 - mp.quad(lambda t: mp.besselk(0, t), [x, x + 5, x + 20, mp.inf])
        print(f"{x!r}, {float(ti)!r}, {float(tk)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, ref_i0, ref_k0] = GENERATE(
        test_case{12.0, 19894.18273673456, 1.5707942066465768}, test_case{18.91, 15508847.272708228, 1.570796325082015},
        test_case{20.0, 44758598.91385575, 1.570796326234015},
        test_case{45.0, 2.1075225881890097e+18, 1.5707963267948966}
    );
    double i0int, k0int;
    xsf::it1i0k0(x, i0int, k0int);

    const double rtol = 5e-15;
    const auto rel_error_i0 = xsf::extended_relative_error(i0int, ref_i0);
    const auto rel_error_k0 = xsf::extended_relative_error(k0int, ref_k0);
    CAPTURE(x, i0int, ref_i0, rel_error_i0, k0int, ref_k0, rel_error_k0, rtol);
    REQUIRE(rel_error_i0 <= rtol);
    REQUIRE(rel_error_k0 <= rtol);
}
