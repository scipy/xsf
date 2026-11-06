#include "../testing_utils.h"

#include <xsf/cpu/exp.h>
#include <xsf/erf.h>

// Regression tests for scipy/scipy#26556: exp2, exp10 and erfc returned 0
// instead of a subnormal result, e.g. exp2(-1024.0000000000002), exp10(-310)
// and erfc(27) were all 0.0.

/*
The expected values are exp2(x), exp10(x) and erfc(x) computed with mpmath and
correctly rounded to double, including in the subnormal range. This is the
script:

import math
import mpmath

mpmath.mp.prec = 300


def to_double(v):
    # Round to the nearest double; mpmath.nint rounds half to even.
    tiny = mpmath.mpf(2) ** -1074
    if abs(v) < mpmath.mpf(2) ** -1022:
        return int(mpmath.nint(v / tiny)) * 5e-324
    return float(v)


cases = {
    "exp2": (lambda x: mpmath.power(2, x), [
        -1022.0, -1023.5, -1024.0, math.nextafter(-1024.0, -math.inf), -1074.0,
        -1074.5, math.nextafter(-1075.0, 0), -1075.0, -1076.0,
    ]),
    "exp10": (lambda x: mpmath.power(10, x), [
        -307.0, -308.0, -310.0, -323.3, -323.60724533877976, -323.6072453387798,
        -324.0,
    ]),
    "erfc": (mpmath.erfc, [
        26.5, 26.543258454250978, 26.54325845425098, 27.0, 27.2,
        27.226017111108362, 27.226017111108366, 28.0,
    ]),
}
for name, (f, xs) in cases.items():
    for x in xs:
        print(f"{name}({x!r}) = {to_double(f(mpmath.mpf(x))).hex()}")
*/

namespace {

// |out - desired| in units of the spacing of doubles just above |desired|,
// which is the smallest subnormal 2**-1074 when desired is subnormal or 0.
double ulp_error(double out, double desired) {
    const double abs_desired = std::abs(desired);
    const double ulp = std::nextafter(abs_desired, std::numeric_limits<double>::infinity()) - abs_desired;
    return std::abs(out - desired) / ulp;
}

} // namespace

TEST_CASE("exp2 subnormal results", "[exp2][xsf_tests]") {
    // 2**x is subnormal for -1075 < x < -1022; it rounds to 0 for x <= -1075.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, desired, max_ulps] = GENERATE(
        test_case{-1022.0, 0x1p-1022, 0}, test_case{-1023.5, 0x0.5a827999fcef3p-1022, 1},
        test_case{-1024.0, 0x0.4p-1022, 0}, test_case{-1024.0000000000002, 0x0.3ffffffffff4fp-1022, 1},
        test_case{-1074.0, 0x0.0000000000001p-1022, 0}, test_case{-1074.5, 0x0.0000000000001p-1022, 0},
        test_case{-1074.9999999999998, 0x0.0000000000001p-1022, 0}, test_case{-1075.0, 0.0, 0},
        test_case{-1076.0, 0.0, 0}
    );

    const double out = xsf::exp2(x);
    const double error = ulp_error(out, desired);
    CAPTURE(x, out, desired, error, max_ulps);
    REQUIRE(error <= max_ulps);
}

TEST_CASE("exp10 subnormal results", "[exp10][xsf_tests]") {
    // 10**x is subnormal for log10(2**-1075) = -323.607245... < x < -307.65;
    // -323.60724533877976 is the smallest double with a nonzero result.
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, desired, max_ulps] = GENERATE(
        test_case{-307.0, 0x1.1fa182c40c60dp-1020, 1}, test_case{-308.0, 0x0.730d67819e8d2p-1022, 1},
        test_case{-310.0, 0x0.012688b70e62bp-1022, 1}, test_case{-323.3, 0x0.0000000000001p-1022, 0},
        test_case{-323.60724533877976, 0x0.0000000000001p-1022, 0}, test_case{-323.6072453387798, 0.0, 0},
        test_case{-324.0, 0.0, 0}
    );

    const double out = xsf::exp10(x);
    const double error = ulp_error(out, desired);
    CAPTURE(x, out, desired, error, max_ulps);
    REQUIRE(error <= max_ulps);
}

TEST_CASE("erfc subnormal results", "[erfc][xsf_tests]") {
    // erfc(x) < DBL_MIN for x >= 26.54325845425098, and erfc(x) rounds to 0
    // for x > 27.226017111108362. erfc(26.5) and erfc(26.543258454250978)
    // are normal and use the unchanged Cephes branch, which evaluates
    // exp(-x * x) with a rounded argument; its error near x = 26.5 is up to
    // about 510 ulp (relative error 1.1e-13).
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [x, desired, max_ulps] = GENERATE(
        test_case{26.5, 0x1.3df6725a60cf5p-1019, 600}, test_case{26.543258454250978, 0x1.0000000000332p-1022, 600},
        test_case{26.54325845425098, 0x0.fffffffffffe0p-1022, 8}, test_case{27.0, 0x0.0000000019e0fp-1022, 1},
        test_case{27.2, 0x0.0000000000002p-1022, 0}, test_case{27.226017111108362, 0x0.0000000000001p-1022, 0},
        test_case{27.226017111108366, 0.0, 0}, test_case{28.0, 0.0, 0}
    );

    const double out = xsf::erfc(x);
    const double error = ulp_error(out, desired);
    CAPTURE(x, out, desired, error, max_ulps);
    REQUIRE(error <= max_ulps);
}
