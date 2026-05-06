#include "../testing_utils.h"
#include <xsf/cephes/j0.h>

TEST_CASE("j0 right tail gh-large-input", "[j0][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    // Reference values computed with mpmath with 1000 digits of precision:
    //
    // import math
    // import mpmath as mp
    // mp.mp.dps = 1000
    // xs = [
    //     math.nextafter(10.0, 0.0), 10.0, math.nextafter(10.0, math.inf),
    //     100.0, 1e4, 1e8, 1e12, 1e13, 1e15, 1e20,
    // ]
    // for x in xs:
    //     print(x, mp.nstr(mp.besselj(0, x), 18))
    auto [x, ref, rtol] = GENERATE(
        test_case{9.999999999999998, -0.24593576445134826, 5e-15}, test_case{10.0, -0.24593576445134834, 5e-15},
        test_case{10.000000000000002, -0.24593576445134841, 5e-15}, test_case{100.0, 0.019985850304223122, 5e-15},
        test_case{1e4, -0.0070961603533888015, 5e-15}, test_case{1e8, 3.206029534041208e-05, 1e-15},
        test_case{1e12, 1.0167125050040682e-07, 1e-15}, test_case{1e13, 1.1926484739665653e-07, 1e-15},
        test_case{1e15, 6.156638646885021e-09, 1e-15}, test_case{1e20, 6.698009040703424e-12, 2e-15}
    );
    const double w = xsf::cephes::j0(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
