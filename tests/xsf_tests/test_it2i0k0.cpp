#include "../testing_utils.h"
#include <tuple>
#include <xsf/cpu/bessel.h>

TEST_CASE("it2i0k0 large x", "[it2i0k0][xsf_tests]") {
    // import mpmath as mp
    // mp.mp.dps = 40

    // def it2i0k0(x):
    //     X = mp.mpf(x)
    //     with mp.workdps(int(x / 2.3) + 40):
    //         ti = +mp.nsum(lambda k: (X / 2) ** (2 * k) / (2 * k * mp.factorial(k) ** 2), [1, mp.inf])
    //     # K0 integrand scaled by exp(x) so the quadrature works with values near 1
    //     g = lambda u: mp.besselk(0, X + u) * mp.exp(X) / (X + u)
    //     tk = mp.quad(g, [0, 1, 5, 20, 60, 150, mp.inf]) * mp.exp(-X)
    //     return ti, tk

    // for x in [30, 40, 100, 300]:
    //     ti, tk = it2i0k0(x)
    //     print(x, float(ti), float(tk))
    using test_case = std::tuple<double, double, double, double>;
    auto [x, ref_i0, ref_k0, rtol] = GENERATE(
        test_case{30.0, 27486628817.36698, 6.780268362546424e-16, 1e-6},
        test_case{40.0, 387326525512364.4, 2.0242046145399133e-20, 1e-7},
        test_case{100.0, 1.090289679929662e+40, 4.5885214327746765e-47, 1e-7},
        test_case{300.0, 1.4994738835332027e+126, 1.2350782737918484e-134, 1e-7}
    );
    double i0int, k0int;
    xsf::it2i0k0(x, i0int, k0int);

    const auto rel_error_i0 = xsf::extended_relative_error(i0int, ref_i0);
    const auto rel_error_k0 = xsf::extended_relative_error(k0int, ref_k0);

    CAPTURE(x, i0int, ref_i0, rtol, rel_error_i0);
    CAPTURE(x, k0int, ref_k0, rtol, rel_error_k0);
    REQUIRE(rel_error_i0 <= rtol);
    REQUIRE(rel_error_k0 <= rtol);
}
