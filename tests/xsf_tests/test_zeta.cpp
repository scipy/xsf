#include "../testing_utils.h"
#include <xsf/zeta.h>

TEST_CASE("zeta(x, q=1) == riemann_zeta(x) for float, double and complex", "[zeta][xsf_tests]") {
    SECTION("double") {
        for (double x : linspace(-20.0, 20.0, 81)) {
            REQUIRE(xsf::zeta(x, 1.0) == xsf::riemann_zeta(x));
        }
    }

    SECTION("float") {
        for (float x : linspace(-20.0f, 20.0f, 81)) {
            REQUIRE(xsf::zeta(x, 1.0f) == xsf::riemann_zeta(x));
        }
    }

    SECTION("complex") {
        using std::complex;
        for (double re : linspace(-5.0, 5.0, 21)) {
            for (double im : linspace(-10.0, 10.0, 41)) {
                complex<double> z(re, im);
                REQUIRE(xsf::zeta(z, 1.0) == xsf::riemann_zeta(z));
            }
        }
    }
}
