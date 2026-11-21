#include "../testing_utils.h"

#include <xsf/cpu/specfun/specfun.h>

#include <vector>

TEST_CASE("specfun::cerzo many zeros", "[cerzo][xsf_tests]") {
    // cerf overflowed (intermediate cosh(n y) = inf) for large |Im z|, so zeros 77 and up were NaN. Also, the deflation
    // product prod(z - zo[i]) overflowed, so zeros 306 and up were NaN even with cerf fixed.

    /*
    import mpmath

    mpmath.mp.dps = 40

    for n in [77, 100, 400]:
        # Newton from the asymptotic start value of the n-th zero
        pu = mpmath.sqrt(mpmath.pi * (4 * n - mpmath.mpf(1) / 2))
        pv = mpmath.pi * mpmath.sqrt(2 * n - mpmath.mpf(1) / 4)
        z0 = mpmath.mpc(
            0.5 * pu - 0.5 * mpmath.log(pv) / pu,
            0.5 * pu + 0.5 * mpmath.log(pv) / pu,
        )
        z = mpmath.findroot(mpmath.erf, z0)
        print(f"{n}: {{{float(z.real)!r}, {float(z.imag)!r}}}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, std::complex<double>>;
    auto [n, desired] = GENERATE(
        test_case{77, {15.481740840635416, 15.599573700766374}},
        test_case{100, {17.659970331721766, 17.767043069336868}},
        test_case{400, {35.411898828807104, 35.475194975359294}}
    );

    std::vector<std::complex<double>> zo(n);
    xsf::specfun::cerzo(n, zo.data());
    const double error = xsf::extended_relative_error(zo[n - 1], desired);
    const double rtol = 5e-15;
    CAPTURE(n, zo[n - 1], desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("specfun::cerf large imaginary part", "[cerf][xsf_tests]") {
    // The series terms exp(-n^2 / 4) * cosh(n y) overflowed (cosh(n y) is inf while the Gaussian factor is not 0) for y
    // = 20.2, giving NaN, and the 100 terms were too few for y = 50. At 0.5108... + 2j the imaginary series stopped
    // too early, giving +0.545 instead of -1.508.

    /*
    import mpmath

    mpmath.mp.dps = 40

    for z in [mpmath.mpc(20, 20.2), mpmath.mpc(50, 50), mpmath.mpc(0.5108369254183123, 2)]:
        w = mpmath.erf(z)
        print(f"{{{float(z.real)!r}, {float(z.imag)!r}}}: {{{float(w.real)!r}, {float(w.imag)!r}}}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<std::complex<double>, std::complex<double>>;
    auto [z, desired] = GENERATE(
        test_case{{20.0, 20.2}, {11.439694678497254, -60.68845812593305}},
        test_case{{50.0, 50.0}, {0.9935538435582172, -0.004702025221298015}},
        test_case{{0.5108369254183123, 2.0}, {13.614680673133714, -1.5084131364950326}}
    );

    std::complex<double> cer, cder;
    xsf::specfun::cerf(z, &cer, &cder);
    const double error = xsf::extended_relative_error(cer, desired);
    const double rtol = 1e-12;
    CAPTURE(z, cer, desired, error, rtol);
    REQUIRE(error <= rtol);
}
