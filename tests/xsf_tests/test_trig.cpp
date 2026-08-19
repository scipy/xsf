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

TEST_CASE("tandg IEEE infinity sign scipy/20731", "[tandg][xsf_tests]") {
    double x, result;
    for (int n : {-2, -1, 0, 1, 2}) {
        x = (2 * n + 1) * 90.0;
        result = xsf::tandg(x);
        CAPTURE(n, x, result);
        REQUIRE(std::isinf(result));
        REQUIRE(std::signbit(result) == (n % 2 != 0));
    }
}

TEST_CASE("cotdg IEEE infinity sign scipy/20731", "[cotdg][xsf_tests]") {
    double x, result;
    for (double n : {-2.0, -1.0, -0.0, 0.0, 1.0, 2.0}) {
        x = n * 180.0;
        result = xsf::cotdg(x);
        CAPTURE(n, x, result);
        REQUIRE(std::isinf(result));
        if (std::signbit(n)) {
            REQUIRE(std::signbit(result) == (static_cast<int>(n) % 2 == 0));
        } else {
            REQUIRE(std::signbit(result) == (static_cast<int>(n) % 2 != 0));
        }
    }
}

/* These returned 0.0 for |x| > 1e14 because `x - 45 * floor(x / 45)` loses
 * accuracy there. An exact fmod by 360 first removes the limit. scipy/scipy#20723
 */
TEST_CASE("degree trig large arguments scipy/20723", "[sindg][cosdg][tandg][cotdg][xsf_tests]") {
    /* {x, x mod 360}. The reduction is exact, so the two must agree exactly. */
    std::vector<std::pair<double, double>> cases{
        {1e14 + 1.0, 281.0}, /* just past the old 1e14 cutoff */
        {2e14, 200.0},       /* the example in the issue */
        {1e15, 280.0},
        {9007199254740991.0, 31.0}, /* 2**53 - 1, the last odd integer */
        {1e16, 280.0},              /* beyond 2**53 */
        {1e17, 280.0},
        {1e18, 280.0},
        {1e300, 0.0},
        {360000000000030.0, 30.0},
        {360000000000090.0, 90.0}, /* pole of tandg */
        {360000000000123.0, 123.0},
        {360000000000180.0, 180.0}, /* pole of cotdg */
    };
    for (auto [x, reduced] : cases) {
        for (double sign : {1.0, -1.0}) {
            double y = sign * x;
            double y_reduced = sign * reduced;
            CAPTURE(x, reduced, sign);
            REQUIRE(xsf::sindg(y) == xsf::sindg(y_reduced));
            REQUIRE(xsf::cosdg(y) == xsf::cosdg(y_reduced));
            REQUIRE(xsf::tandg(y) == xsf::tandg(y_reduced));
            REQUIRE(xsf::cotdg(y) == xsf::cotdg(y_reduced));
        }
    }
    /* sindg(2e14) used to be 0.0; sindg(2e14 mod 360) = sindg(200) */
    REQUIRE(std::abs(xsf::sindg(2e14) + 0.3420201433256687) < 1e-15);
}
