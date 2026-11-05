#include "../testing_utils.h"

#include <xsf/cpu/par_cyl.h>

TEST_CASE("pbdv and pbvv derivative when |v| + 1 rounds up", "[pbdv][pbvv][xsf_tests]") {
    // Regression test for gh-401: the work buffer was one element too short.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import math

    import mpmath

    mpmath.mp.dps = 50


    def pbdv_d(v, x):
        return mpmath.diff(lambda t: mpmath.pcfd(v, t), x)


    def pbvv_d(v, x):
        return mpmath.diff(lambda t: mpmath.pcfv(-mpmath.mpf(v) - 0.5, t), x)


    for fn, v, x in [
        (pbdv_d, math.nextafter(2, 0), 2.0),
        (pbdv_d, -math.nextafter(2, 0), 2.0),
        (pbvv_d, math.nextafter(1, 0), 2.0),
        (pbvv_d, -math.nextafter(1, 0), 2.0),
    ]:
        print(f"{fn.__name__}({v!r}, {x!r}) = {float(fn(v, x))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<bool, double, double, double>;
    auto [is_pbdv, v, x, desired] = GENERATE(
        test_case{true, 1.9999999999999998, 2.0, 0.367879441171442},
        test_case{true, -1.9999999999999998, 2.0, -0.09715978862055018},
        test_case{false, 0.9999999999999999, 2.0, 0.780820901800019},
        test_case{false, -0.9999999999999999, 2.0, 2.168875102838455}
    );

    double f, d;
    if (is_pbdv) {
        xsf::pbdv(v, x, f, d);
    } else {
        xsf::pbvv(v, x, f, d);
    }
    const double error = xsf::extended_relative_error(d, desired);
    const double rtol = 1e-14;
    CAPTURE(is_pbdv, v, x, d, desired, error, rtol);
    REQUIRE(error <= rtol);
}
