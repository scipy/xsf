#include "../testing_utils.h"
#include <cmath>
#include <complex>
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

TEST_CASE("sinpi and cospi complex conjugate symmetry for large imaginary part", "[sinpi][cospi][xsf_tests]") {
    // For pi*|Im(z)| >= 700 the sign of Im(z) was dropped. sinpi and cospi
    // satisfy f(conj(z)) == conj(f(z)).
    for (double y : {225.0, 240.0, 500.0}) {
        for (double x : {0.25, -0.75, 0.5, 1.0}) {
            std::complex<double> z(x, y);
            auto sin_pos = xsf::sinpi(z);
            auto sin_neg = xsf::sinpi(std::conj(z));
            auto cos_pos = xsf::cospi(z);
            auto cos_neg = xsf::cospi(std::conj(z));
            CAPTURE(x, y, sin_pos, sin_neg, cos_pos, cos_neg);
            REQUIRE(sin_neg == std::conj(sin_pos));
            REQUIRE(std::signbit(sin_neg.imag()) != std::signbit(sin_pos.imag()));
            REQUIRE(cos_neg == std::conj(cos_pos));
            REQUIRE(std::signbit(cos_neg.imag()) != std::signbit(cos_pos.imag()));
        }
    }
}

TEST_CASE("sinpi and cospi complex large imaginary part", "[sinpi][cospi][xsf_tests]") {
    // Reference values computed with mpmath.
    const double a = 3.4129868967024855e+306;
    const double tol = 1e-12;
    auto result = xsf::sinpi(std::complex<double>(0.25, -225.0));
    CAPTURE(result);
    REQUIRE(xsf::extended_relative_error(result, std::complex<double>(a, -a)) <= tol);
    result = xsf::cospi(std::complex<double>(0.25, 225.0));
    CAPTURE(result);
    REQUIRE(xsf::extended_relative_error(result, std::complex<double>(a, -a)) <= tol);
}

TEST_CASE("sinpi and cospi complex zero factor stays zero on overflow", "[sinpi][cospi][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    // y = 240: exp(pi*|y|/2) is finite but the final product overflows.
    // y = 500: exp(pi*|y|/2) itself overflows.
    for (double y : {240.0, 500.0}) {
        CAPTURE(y);
        auto result = xsf::cospi(std::complex<double>(0.5, -y));
        CAPTURE(result);
        REQUIRE(result.real() == 0.0);
        REQUIRE(result.imag() == inf);
        result = xsf::cospi(std::complex<double>(0.5, y));
        CAPTURE(result);
        REQUIRE(result.real() == 0.0);
        REQUIRE(result.imag() == -inf);
        result = xsf::sinpi(std::complex<double>(0.5, y));
        CAPTURE(result);
        REQUIRE(result.real() == inf);
        REQUIRE(result.imag() == 0.0);
        REQUIRE(!std::signbit(result.imag()));
        result = xsf::sinpi(std::complex<double>(0.5, -y));
        CAPTURE(result);
        REQUIRE(result.real() == inf);
        REQUIRE(result.imag() == 0.0);
        REQUIRE(std::signbit(result.imag()));
    }
}
