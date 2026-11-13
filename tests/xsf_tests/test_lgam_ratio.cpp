#include "../testing_utils.h"

#include <xsf/beta.h>
#include <xsf/cephes/poch.h>

// Regression tests for the loss of precision in lgam(a) + lgam(b) - lgam(a + b) for large
// arguments in beta and betaln, and in lgam(a + m) - lgam(a) in poch (scipy/scipy#26423).
// These go through cephes::detail::lgam_ratio, a port of jorenham/xsf-rust#415 and #420.
// The relative errors before the fix are given in the comments.

/*
If you want to check the expected values below with mpmath, this is the script
(the precision is chosen so that lgamma keeps >= 150 bits after the cancellation and
a + b is exact):

import math
import mpmath

def prec_for(*xs):
    m = max(abs(x) for x in xs)
    e = max(math.frexp(x)[1] for x in xs) - min(math.frexp(x)[1] for x in xs)
    return max(int(math.log2(m * abs(math.log(m)) + m + 10)) + 200, e + 120, 300)

def beta(a, b):
    mpmath.mp.prec = prec_for(a, b)
    return float(mpmath.beta(a, b))

def betaln(a, b):
    mpmath.mp.prec = prec_for(a, b)
    a, b = mpmath.mpf(a), mpmath.mpf(b)
    return float(mpmath.re(mpmath.loggamma(a) + mpmath.loggamma(b) - mpmath.loggamma(a + b)))

def poch(a, m):
    mpmath.mp.prec = prec_for(a, m)
    return float(mpmath.rf(a, m))

print(beta(999999.0, 0.5).hex())  # etc.

The values were also checked against MPFR (gmpy2.lgamma at the same precision).
*/

TEST_CASE("beta large arguments", "[beta][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, b, desired, rtol] = GENERATE(
        // a > 0, |b| < a / 8
        test_case{999999.0, 0.5, 0x1.d0a3705328066p-10, 5e-15},     // 0.0017724549586901838, was 9.4e-11
        test_case{0.5, 999999.0, 0x1.d0a3705328066p-10, 5e-15},     // was 9.4e-11
        test_case{300000.0, -2.5, -0x1.530c2d22cf6acp+45, 5e-15},   // -46598342334957.34, was 4.1e-10
        test_case{500.25, 20.75, 0x1.46c0ca4433e03p-127, 4e-14},    // 7.501882951310635e-39, was 1.3e-13
        test_case{5000.0, 10.5, 0x1.0dfaebd45c60ap-109, 5e-14},     // 1.6248839887769473e-33, was 2.4e-12
        test_case{100000.0, -30.5, -0x1.0a8c4465caa10p+398, 2e-13}, // -6.721614555437173e+119, was 1.9e-10
        // a < 0, |b| < |a| / 8 (reflection)
        test_case{-1000000.25, 0.5, -0x1.d0a355ae3c1ecp-10, 5e-15}, // -0.0017724534077921917, was 3.5e-10
        // a > 0, -a <= b <= -a / 2 (reflection)
        test_case{5000.25, -4999.5, 0x1.3831f9d3bee57p-2, 5e-15}, // 0.3048781428803741, was 3.7e-12
        test_case{180.5, -179.75, 0x1.fab5246732ba8p-1, 5e-15}    // 0.989663255307609, was 3.0e-14
    );

    const double out = xsf::beta(a, b);
    const double error = xsf::extended_relative_error(out, desired);
    CAPTURE(a, b, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("betaln large arguments", "[betaln][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, b, desired, rtol] = GENERATE(
        // a > 0, |b| < a / 8
        test_case{1e100, 1e94, -0x1.1c1aba8044549p+316, 4e-15},     // -1.4815511057964107e+95, was 6.8e-9
        test_case{1e94, 1e100, -0x1.1c1aba8044549p+316, 4e-15},     // was 6.8e-9
        test_case{10000.0, 0.25, -0x1.03b9c210ec4afp+0, 4e-15},     // -1.0145531932178462, was 3.3e-12
        test_case{1e6, 1e3, -0x1.ee6ca1a9a3b65p+12, 4e-15},         // -7910.78946842146, was 1.4e-13
        test_case{1e20, 1e17, -0x1.5f327399989f7p+59, 4e-15},       // -7.908255112398755e+17, was 4.1e-13
        test_case{1e300, 1e295, -0x1.87d813caf151ep+983, 2e-14},    // -1.2512930464953562e+296, was 4.5e-10
        test_case{2000.5, -100.25, 0x1.8c08b23156e0cp+8, 4e-15},    // 396.03396900530856, was 6.3e-15
        test_case{1e10, -123456.789, 0x1.72cc0bed99f88p+20, 4e-15}, // 1518784.745508166, was 1.6e-11
        // a < 0, |b| < |a| / 8 (reflection)
        test_case{-1000000.25, 0.5, -0x1.95770a13a02d5p+2, 4e-15}, // -6.3353905860573905, was 5.5e-11
        test_case{-2000.75, 3.5, -0x1.966bc37ef36b7p+4, 4e-15},    // -25.4013094862955, was 4.6e-14
        // a > 0, -a <= b <= -a / 2 (reflection)
        test_case{5000.25, -4999.5, -0x1.30167c7e54188p+0, 4e-15} // -1.1878431137484728, was 3.1e-12
    );

    const double out = xsf::betaln(a, b);
    const double error = xsf::extended_relative_error(out, desired);
    CAPTURE(a, b, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("poch moderately large arguments", "[poch][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [a, m, desired] = GENERATE(
        // 10 < a <= 1e4
        test_case{5878.8618471964755, -0.029080395984806877, 0x1.8dcb5573dd045p-1}, // 0.7769419387755155, was 7.9e-13
        test_case{150.5, 0.75, 0x1.5789290d6bed0p+5},                               // 42.94197283254255, was 1.6e-13
        test_case{50.25, -0.5, 0x1.231693873d72ep-3},                               // 0.14213290459127265, was 1.4e-14
        test_case{10.5, 0.3, 0x1.0090ad578484cp+1},                                 // 2.004415195210674
        test_case{9999.5, 0.999, 0x1.359e979ed8fbfp+13},                            // 9907.824033446492, was 5.9e-12
        test_case{1234.5, 7.25, 0x1.65823c35d9853p+74}, // 2.6379488276984808e+22, was 2.3e-14
        // a < -10 (reflection)
        test_case{-96483701341.52873, -0.0010688260996181021, 0x1.f27fa4bcfc1dep-1}, // 0.9736300926870702, was 3.6e-4
        test_case{-10000.25, 0.5, -0x1.90028f5b1c849p+6},                            // -100.00249998437461, was 5.2e-13
        test_case{-50.75, 2.5, 0x1.142735c4aac15p+14},                               // 17673.802508037614, was 4.3e-14
        test_case{-20.3, -0.6, 0x1.ae3835b6152abp-2}                                 // 0.4201362984070253, was 1.7e-14
    );

    const double out = xsf::cephes::poch(a, m);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 4e-15;
    CAPTURE(a, m, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
