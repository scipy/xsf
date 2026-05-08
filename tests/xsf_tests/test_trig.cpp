#include "../testing_utils.h"
#include <cmath>
#include <xsf/trig.h>

void check_zeros(double (*f)(double)) {
    std::vector<double> zeros{-0.0, 0.0};
    for (double x : zeros) {
        double result = f(x);
        CAPTURE(x, result);
        REQUIRE(std::signbit(result) == std::signbit(x));
        REQUIRE(result == 0.0);
    }
}

TEST_CASE("sindg IEEE scipy/20731", "[sindg][xsf_tests]") { check_zeros(&xsf::sindg); }

TEST_CASE("sinpi IEEE scipy/20731", "[sinpi][xsf_tests]") { check_zeros(&xsf::sinpi); }

TEST_CASE("tandg IEEE scipy/20731", "[tandg][xsf_tests]") { check_zeros(&xsf::tandg); }
