#include "../testing_utils.h"
#include <xsf/cephes/j0.h>

TEST_CASE("j0 right tail gh-large-input", "[j0][xsf_tests]") {
    // Reference value computed with mpmath: mp.besselj(0, 1e15)
    const double x = 1e15;
    const double ref = 6.156638646885021e-09;
    const double rtol = 1e-13;

    const double w = xsf::cephes::j0(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
