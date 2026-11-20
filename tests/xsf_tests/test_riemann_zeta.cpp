#include "../testing_utils.h"

#include <xsf/zeta.h>

TEST_CASE("riemann_zeta complex z near 0", "[riemann_zeta][xsf_tests]") {
    // Regression test: the reflection formula evaluated zeta(1 - z) next to its pole,
    // losing accuracy for non-real abs(z) < 1/2.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    from mpmath import mp, mpc, zeta

    mp.dps = 50
    for z in [1e-6 - 2e-6j, 1e-5 + 1e-5j, -1e-3 + 3e-3j]:
        w = zeta(mpc(z))
        print(f"{z!r}: ({float(w.real)!r}, {float(w.imag)!r})")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<std::complex<double>, std::complex<double>>;
    auto [z, desired] = GENERATE(
        test_case{{1e-6, -2e-6}, {-0.5000009189355237, 1.8378810791202557e-06}},
        test_case{{1e-5, 1e-5}, {-0.50000918938533, -9.189585969693889e-06}},
        test_case{{-1e-3, 3e-3}, {-0.4990730620890677, -0.002750778612088271}}
    );

    const auto out = xsf::riemann_zeta(z);
    // The imaginary part is checked separately: it is much smaller than the real part.
    const double re_error = xsf::extended_relative_error(out.real(), desired.real());
    const double im_error = xsf::extended_relative_error(out.imag(), desired.imag());
    const double rtol = 5e-15;
    CAPTURE(z, out, desired, re_error, im_error, rtol);
    REQUIRE(re_error <= rtol);
    REQUIRE(im_error <= rtol);
}
