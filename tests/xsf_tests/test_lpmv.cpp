#include "../testing_utils.h"
#include <tuple>
#include <xsf/specfun/specfun.h>
#include <cmath>

TEST_CASE("lpmv endpoint x=-1 gives signed infinity for m=0", "[pmv][xsf_tests]") {
  const double m = 0.0;
  const double x = -1.0;
  constexpr double HUGE_NEG = -0.5 * std::numeric_limits<double>::max();
  constexpr double HUGE_POS = 0.5 * std::numeric_limits<double>::max();

  SECTION("nu = 0.5 -> -inf") {
    const double nu = 0.5;
    const double w = xsf::specfun::lpmv(x, m, nu);

    CAPTURE(x, m, nu, w);
    REQUIRE(w < -1e298);
    REQUIRE(std::signbit(w));  // negative infinity
  }

  SECTION("nu = 1.5 -> +inf") {
    const double nu = 1.5;
    const double w = xsf::specfun::lpmv(x, m, nu);

    CAPTURE(x, m, nu, w);
    REQUIRE(w > 1e298);
    REQUIRE(!std::signbit(w)); // positive infinity
  }
}
