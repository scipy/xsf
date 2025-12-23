#include "../testing_utils.h"
#include <tuple>
#include <xsf/zeta.h>
#include <catch2/generators/catch_generators_range.hpp>

TEST_CASE("zeta(x, q=1) matches riemann_zeta for all types", "[zeta][xsf_tests]") {
    SECTION("double") {
        double x = GENERATE(range(-10.0, 10.0, 0.1));
        REQUIRE(xsf::zeta(x, 1.0) == xsf::riemann_zeta(x));
    }
    
    SECTION("float") {
        float x = GENERATE(range(-10.0f, 10.0f, 0.5f));
        REQUIRE(xsf::zeta(x, 1.0f) == xsf::riemann_zeta(x));
    }
    
    SECTION("complex") {
        using std::complex;
        double re = GENERATE(range(0.5, 5.0, 0.5));
        double im = GENERATE(range(-2.0, 2.0, 0.5));
        complex<double> z(re, im);
        REQUIRE(xsf::zeta(z, 1.0) == xsf::riemann_zeta(z));
    }
}
