#include "../testing_utils.h"

#include <xsf/kelvin.h>

TEST_CASE("kei is continuous down to the smallest subnormal", "[kelvin][kei][xsf_tests]") {
    // klvna computes log(x / 2) for the series branch. That halving is exact while
    // x is normal, but not once x is subnormal: at the smallest subnormal the
    // quotient underflows to zero and the logarithm becomes -inf. It multiplies a
    // series term that has itself underflowed to zero, so the product is
    // inf * 0 = NaN -- while both neighbours, x == 0.0 and the next representable
    // value up, are correct.
    const double rtol = 1e-15;
    // The same literal klvna uses, so this checks the value the implementation
    // targets rather than a constant from somewhere else.
    const double expected = -0.25 * 3.141592653589793;

    const auto x = GENERATE(std::numeric_limits<double>::denorm_min(), 1e-323, 1e-320, 1e-300, 1e-200, 1e-100);

    const auto out = xsf::kei(x);
    const auto rel_error = xsf::extended_relative_error(out, expected);
    CAPTURE(x, out, expected, rtol, rel_error);
    REQUIRE(std::isfinite(out));
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("kelvin functions agree either side of the subnormal boundary", "[kelvin][xsf_tests]") {
    // The four outputs the series produces for a vanishing argument. Pinned
    // together so a future change to the small-x path cannot fix one and leave
    // the rest to drift.
    const double denorm_min = std::numeric_limits<double>::denorm_min();

    REQUIRE(xsf::ber(denorm_min) == 1.0);
    REQUIRE(xsf::bei(denorm_min) == 0.0);
    REQUIRE(std::isfinite(xsf::kei(denorm_min)));
    // ker diverges as x -> 0, so only its sign and finiteness are meaningful.
    REQUIRE(xsf::ker(denorm_min) > 0.0);
}

TEST_CASE("ker is strictly decreasing across the subnormal range", "[kelvin][ker][xsf_tests]") {
    // The inexact halving reaches further than the underflow does. An odd multiple
    // of 2^-1074 rounds onto the neighbouring even one, so before the fix ker
    // returned the same value for three consecutive representable inputs -- worst
    // relative error 3.87e-04 over the first 25 multiples, against 2.17e-16 after.
    //
    // ker is strictly decreasing, so equal outputs for distinct inputs are a defect
    // on their own terms and this needs no reference table.
    //
    // The range matters and the bound is deliberate. Here the inputs differ by a
    // ratio, so ker moves by log(1 + 1/k) >= 0.039 between neighbours -- ker(3u) and
    // ker(4u) are 2.5e12 ulps apart, and equal outputs cannot be right. Higher up
    // that stops holding: adjacent doubles near DBL_MIN differ by a relative 2.2e-16,
    // so the true ker values are 2.2e-16 apart on a result of 708.5, which is 512x
    // below one ulp there. They round to the same double, and equal outputs are then
    // correct rather than a defect. So this invariant is asserted only over small
    // multiples of denorm_min, where the spacing is coarse enough to resolve.
    const double u = std::numeric_limits<double>::denorm_min();
    double previous = std::numeric_limits<double>::infinity();
    for (int k = 1; k <= 25; k++) {
        const double x = k * u;
        const double out = xsf::ker(x);
        CAPTURE(k, x, out, previous);
        REQUIRE(std::isfinite(out));
        REQUIRE(out < previous);
        previous = out;
    }
}

TEST_CASE("ker subnormal-input", "[kelvin][ker][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    // Reference values computed with mpmath with 60 digits of precision:
    //
    // import mpmath as mp
    // mp.mp.dps = 60
    // u = 4.9406564584124654e-324
    // for k in [1, 2, 3, 4, 5, 7, 9, 17, 25]:
    //     print(k * u, mp.nstr(mp.ker(0, mp.mpf(k * u)), 18))
    auto [x, ref] = GENERATE(
        test_case{4.9406564584124654e-324, 744.556003437039675},
        test_case{9.8813129168249309e-324, 743.862856256479729},
        test_case{1.4821969375237396e-323, 743.457391148371565},
        test_case{1.9762625833649862e-323, 743.169709075919784},
        test_case{2.4703282292062327e-323, 742.946565524605574},
        test_case{3.4584595208887258e-323, 742.610093287984361},
        test_case{4.4465908125712189e-323, 742.358778859703455},
        test_case{8.3991159793011913e-323, 741.722790092983459},
        test_case{1.2351641146031164e-322, 741.337127612171474},
        test_case{9.9998886718268301e-321, 736.943172406632319}, test_case{2.2250738585072014e-308, 708.512350047922519}
    );
    const double w = xsf::ker(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rel_error);
    REQUIRE(rel_error <= 5e-16);
}
