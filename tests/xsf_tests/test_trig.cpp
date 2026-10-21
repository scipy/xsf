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

TEST_CASE("sinpi and cospi complex large imaginary part", "[sinpi][cospi][xsf_tests]") {
    // Reference values from mpmath with 1000 digits of precision:
    //
    // import mpmath as mp
    // mp.mp.dps = 1000
    // print(mp.nstr(mp.sinpi(0.25-225.0j),18))
    // print(mp.nstr(mp.cospi(0.25+225.0j),18))
    const double a = 3.41298689670248524e+306;
    const double rtol = 1e-12;
    std::complex<double> w;

    w = xsf::sinpi(std::complex<double>(0.25, -225.0));
    CAPTURE(w);
    REQUIRE(xsf::extended_relative_error(w, std::complex<double>(a, -a)) <= rtol);

    w = xsf::cospi(std::complex<double>(0.25, 225.0));
    CAPTURE(w);
    REQUIRE(xsf::extended_relative_error(w, std::complex<double>(a, -a)) <= rtol);
}

TEST_CASE("cospi complex imaginary part overflow sign", "[cospi][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    std::complex<double> w;

    w = xsf::cospi(std::complex<double>(0.25, 1000.0));
    CAPTURE(w);
    REQUIRE(w.real() == inf);
    REQUIRE(w.imag() == -inf);

    w = xsf::cospi(std::complex<double>(0.5, 1000.0));
    CAPTURE(w);
    REQUIRE(w.real() == 0.0);
    REQUIRE(w.imag() == -inf);

    w = xsf::cospi(std::complex<double>(0.0, 1000.0));
    CAPTURE(w);
    REQUIRE(w.real() == inf);
    REQUIRE(w.imag() == 0.0);
}
