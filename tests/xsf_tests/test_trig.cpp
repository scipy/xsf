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

TEST_CASE("sinpi and cospi complex non-finite real part and huge imaginary part", "[sinpi][cospi][xsf_tests]") {
    // Regression test for gh-469: this returned +-inf instead of NaN.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const std::complex<double> z = GENERATE_COPY(std::complex<double>(nan, 1000.0), std::complex<double>(inf, -1000.0));

    const std::complex<double> s = xsf::sinpi(z);
    const std::complex<double> c = xsf::cospi(z);
    CAPTURE(z, s, c);
    REQUIRE(std::isnan(s.real()));
    REQUIRE(std::isnan(s.imag()));
    REQUIRE(std::isnan(c.real()));
    REQUIRE(std::isnan(c.imag()));
}

TEST_CASE("sinpi and cospi complex subnormal real part and large imaginary part", "[sinpi][cospi][xsf_tests]") {
    // Regression test for gh-469: 0.5 * sinpi(x) lost precision for subnormal sinpi(x).
    // x = 113 * 2**-1074, so that sinpi(x) ~ 355 * 2**-1074 is accurate to 1e-7.
    /*
    import mpmath

    mpmath.mp.dps = 50

    x = 113 * 2.0**-1074
    print(repr(float(mpmath.sinpi(mpmath.mpc(x, 225)).real)))
    print(repr(float(mpmath.cospi(mpmath.mpc(x, 225)).imag)))
    */
    SET_FP_FORMAT()
    const std::complex<double> z(5.6e-322, 225.0);
    const double desired = 8.465694494121489e-15;
    const double rtol = 1e-6;

    const double s = xsf::sinpi(z).real();
    const double c = xsf::cospi(z).imag();
    const double error_s = xsf::extended_relative_error(s, desired);
    const double error_c = xsf::extended_relative_error(c, -desired);
    CAPTURE(z, s, c, desired, error_s, error_c, rtol);
    REQUIRE(error_s <= rtol);
    REQUIRE(error_c <= rtol);
}

TEST_CASE("tandg and cotdg near poles", "[tandg][cotdg][xsf_tests]") {
    // Regression test for gh-363: rounding of the argument was amplified near
    // the poles, e.g. cotdg(1e-15) was inf.
    /*
    import mpmath

    mpmath.mp.dps = 50

    for f, x in [("tandg", 89.9999999999), ("cotdg", 1e-10), ("cotdg", 1e-15), ("cotdg", 179.9999999999)]:
        t = mpmath.tan(mpmath.mpf(x) * mpmath.pi / 180)
        print(f, repr(x), repr(float(t if f == "tandg" else 1 / t)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<bool, double, double>;
    auto [cot, x, desired] = GENERATE(
        test_case{false, 89.9999999999, 572947570131.5804}, test_case{true, 1e-10, 572957795130.8232},
        test_case{true, 1e-15, 5.729577951308232e+16}, test_case{true, 179.9999999999, -573029000997.1478}
    );

    const double out = cot ? xsf::cotdg(x) : xsf::tandg(x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-15;
    CAPTURE(cot, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
