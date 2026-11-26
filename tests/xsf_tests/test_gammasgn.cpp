#include <catch2/catch_test_macros.hpp>

#include <xsf/cephes/gamma.h>
#include <xsf/gamma.h>

// Regression: the sign of Gamma for negative non-integer x was taken from
// floor(x) cast to int, which overflows for |x| > INT_MAX. Gamma(x) for x in
// (-n, -n + 1) has the sign (-1)^n.
TEST_CASE("gammasgn large negative x", "[gammasgn][xsf_tests]") {
    REQUIRE(xsf::gammasgn(-2147483648.5) == -1.0); // n = 2^31 + 1
    REQUIRE(xsf::gammasgn(-3000000000.5) == -1.0); // n = 3000000001
    REQUIRE(xsf::gammasgn(-4294967297.5) == 1.0);  // n = 2^32 + 2
    REQUIRE(xsf::gammasgn(-1e12 - 0.25) == -1.0);  // n = 10^12 + 1
    REQUIRE(xsf::gammasgn(-1e15 - 0.5) == -1.0);   // n = 10^15 + 1
}

TEST_CASE("lgam_sgn sign for large negative x", "[lgam_sgn][xsf_tests]") {
    int sign = 0;
    xsf::cephes::detail::lgam_sgn(-3000000000.5, &sign);
    REQUIRE(sign == -1);
    xsf::cephes::detail::lgam_sgn(-4294967297.5, &sign);
    REQUIRE(sign == 1);
    xsf::cephes::detail::lgam_sgn(-1e12 - 0.25, &sign);
    REQUIRE(sign == -1);
}
