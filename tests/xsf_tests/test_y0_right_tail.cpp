#include "../testing_utils.h"
#include <xsf/cephes/j0.h>

TEST_CASE("y0 right tail gh-large-input", "[y0][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
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
    //     print(x, mp.nstr(mp.bessely(0, x), 18))
    auto [x, ref, rtol] = GENERATE(
        test_case{9.999999999999998, 0.055671167283599834, 5e-15}, test_case{10.0, 0.05567116728359939, 5e-15},
        test_case{10.000000000000002, 0.05567116728359895, 5e-15}, test_case{100.0, -0.07724431336508315, 5e-15},
        test_case{1e4, 0.003647805558986606, 5e-15}, test_case{1e8, 0.00007306391165521707, 1e-15},
        test_case{1e12, -7.913802683850949e-7, 1e-15}, test_case{1e13, -2.2234629165383075e-7, 1e-15},
        test_case{1e15, 2.4468665123771323e-8, 1e-15}, test_case{1e20, -7.95068198242545e-11, 2e-15}
    );
    const double w = xsf::cephes::y0(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
