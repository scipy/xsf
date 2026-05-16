#include "../testing_utils.h"
#include <xsf/cephes/j0.h>
#include <xsf/cephes/j1.h>

TEST_CASE("j0 right tail gh-large-input", "[j0][xsf_tests]") {
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

TEST_CASE("j1 right tail gh-large-input", "[j1][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    /*
    Reference values computed with mpmath with 1000 digits of precision:
    import mpmath as mp
    import numpy as np

    mp.mp.dps = 1000

    # Two special cases around 10
    x_cases = [np.nextafter(10, 100), np.nextafter(10, 0.0)]
    # Logspaced values from 1e2 to 1e20
    x_cases += list(np.logspace(2, 20, num=20))

    for x in x_cases:
        y = mp.besselj(1, mp.mpf(str(x)))
        print(f"x = {x}, J1(x) = {float(y)}")
    */
    auto [x, ref, rtol] = GENERATE(
        test_case(std::nextafter(10, 100), 0.043472746168860994, 5e-15),
        test_case(std::nextafter(10, 0.0), 0.04347274616886188, 5e-15), test_case{100.0, -0.07714535201411216, 1e-15},
        test_case{885.8667904100823, -0.020092316011146107, 1e-15},
        test_case{7847.5997035146065, -0.0069653536476014105, 1e-15},
        test_case{69519.27961775605, 0.0029402478891245158, 1e-15},
        test_case{615848.2110660254, 0.0008651783371245079, 1e-15},
        test_case{5455594.781168515, -0.0003415758031327532, 1e-15},
        test_case{48329302.38571752, 0.00011473305310448148, 1e-15},
        test_case{428133239.8719396, 3.647690393737593e-05, 1e-15},
        test_case{3792690190.7322383, 1.0661301225390192e-05, 1e-15},
        test_case{33598182862.83774, 2.535659724115016e-06, 1e-15},
        test_case{297635144163.1313, 1.3721936572672073e-06, 1e-15},
        test_case{2636650898730.3555, 1.801055949369533e-09, 5e-14},
        test_case{23357214690901.215, 6.120507921031749e-08, 1e-15},
        test_case{206913808111479.0, -4.805366849955263e-08, 1e-15},
        test_case{1832980710832430.0, 1.375098982967035e-08, 1e-15},
        test_case{1.6237767391887176e+16, 4.295602361218747e-09, 1e-15},
        test_case{1.438449888287654e+17, -1.968878305887983e-09, 1e-15},
        test_case{1.2742749857031322e+18, 4.3550674582051375e-10, 1e-15},
        test_case{1.1288378916846838e+19, 1.904229128619541e-10, 1e-15}, test_case{1e+20, -7.95068198242545e-11, 1e-15}
    );
    const double w = xsf::cephes::j1(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
