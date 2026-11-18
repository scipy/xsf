#include "../testing_utils.h"

#include <xsf/cephes/hyp2f1.h>

TEST_CASE("cephes::hyp2f1 with b == -a", "[cephes::hyp2f1][xsf_tests]") {
    // Regression test for gh-460 (https://github.com/scipy/xsf/issues/460):
    // hyp2f1(n, -n, c, x) lost all digits because `a` and `b` were not swapped.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for a, c, x in [(195, 1, 0.9), (13, 0.5, 0.9), (195, -200.5, 0.01)]:
        print(f"hyp2f1({a}, {-a}, {c}, {x}) = {float(mpmath.hyp2f1(a, -a, c, x))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double, double>;
    auto [a, b, c, x, desired] = GENERATE(
        test_case{195.0, -195.0, 1.0, 0.9, -0.019154323187421253}, test_case{13.0, -13.0, 0.5, 0.9, 0.4895143985151992},
        test_case{195.0, -195.0, -200.5, 0.01, 6.722964015601578}
    );

    const double out = xsf::cephes::hyp2f1(a, b, c, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-12;
    CAPTURE(a, b, c, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
