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

TEST_CASE("sindg IEEE pm zero scipy/20731", "[sindg][xsf_tests]") { check_zeros(&xsf::sindg); }

TEST_CASE("sinpi IEEE pm zero scipy/20731", "[sinpi][xsf_tests]") { check_zeros(&xsf::sinpi); }

TEST_CASE("tandg IEEE pm zero scipy/20731", "[tandg][xsf_tests]") { check_zeros(&xsf::tandg); }

TEST_CASE("tandg IEEE zero sign at multiples of 180 scipy/20731", "[tandg][xsf_tests]") {
    double x, result;
    for (int n : {1, 2}) {
        x = (2 * n + 1) * 180.0;
        result = xsf::tandg(x);
        CAPTURE(n, x, result);
        REQUIRE(std::signbit(result));
        REQUIRE(result == 0);
        x = (2 * n + 2) * 180.0;
        result = xsf::tandg(x);
        CAPTURE(n, x, result);
        REQUIRE(!std::signbit(result));
        REQUIRE(result == 0);
    }
    for (int n : {-1, -2}) {
        x = (2 * n - 1) * 180.0;
        result = xsf::tandg(x);
        CAPTURE(n, x, result);
        REQUIRE(!std::signbit(result));
        REQUIRE(result == 0);
        x = (2 * n - 2) * 180.0;
        result = xsf::tandg(x);
        CAPTURE(n, x, result);
        REQUIRE(std::signbit(result));
        REQUIRE(result == 0);
    }
}

TEST_CASE("cotdg IEEE zero sign at multiples of 90 scipy/20731", "[cotdg][xsf_tests]") {
    double x, result;
    for (int n : {1, 2}) {
        x = (4 * n - 3) * 90.0;
        result = xsf::cotdg(x);
        CAPTURE(n, x, result);
        REQUIRE(!std::signbit(result));
        REQUIRE(result == 0);
        x = (4 * n - 1) * 90.0;
        result = xsf::cotdg(x);
        CAPTURE(n, x, result);
        REQUIRE(std::signbit(result));
        REQUIRE(result == 0);
    }
    for (int n : {-1, -2}) {
        x = (4 * n + 3) * 90.0;
        result = xsf::cotdg(x);
        CAPTURE(n, x, result);
        REQUIRE(std::signbit(result));
        REQUIRE(result == 0);
        x = (4 * n + 1) * 90.0;
        result = xsf::cotdg(x);
        CAPTURE(n, x, result);
        REQUIRE(!std::signbit(result));
        REQUIRE(result == 0);
    }
}
