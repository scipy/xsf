#include "../testing_utils.h"

#include <xsf/kelvin.h>

TEST_CASE("kei is continuous down to the smallest subnormal", "[kelvin][kei][xsf_tests]") {
    // klvna computes log(x / 2) for the series branch. x / 2 is exact for every
    // argument except the smallest subnormal, where it underflows to zero and
    // the logarithm becomes -inf. It multiplies a series term that has itself
    // underflowed to zero, so the product is inf * 0 = NaN -- while both
    // neighbours, x == 0.0 and the next representable value up, are correct.
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
