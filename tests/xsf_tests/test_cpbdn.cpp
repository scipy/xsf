#include "../testing_utils.h"

#include <xsf/cpu/specfun/specfun.h>

TEST_CASE("specfun::cpbdn n < 0 complex z gh-537", "[cpbdn][xsf_tests]") {
    // Regression test for gh-537 (https://github.com/scipy/xsf/issues/537):
    // D_0 used exp(-z * conj(z) / 4), and D_{-2}(z) and D_{-2}'(z) were wrong or missing for n = -2.
    // The three z cover the three n < 0 branches (Re(z) <= 0, |z| <= 3, |z| > 3).

    /*
    import mpmath as mp

    mp.mp.dps = 40
    for z in [-1+1j, 1+1j, 3+3j]:
        z = mp.mpc(z)
        d = mp.pcfd(-2, z)
        dp = mp.diff(lambda t: mp.pcfd(-2, t), z)
        print(z, repr(float(d.real)), repr(float(d.imag)), repr(float(dp.real)), repr(float(dp.imag)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<std::complex<double>, std::complex<double>, std::complex<double>>;
    auto [z, desired, desired_d] = GENERATE(
        test_case{{-1.0, 1.0}, {1.0538664503741575, -3.106422389505304}, {-0.8547879387310727, 3.7849264397525917}},
        test_case{{1.0, 1.0}, {0.055834798062130854, -0.29509248409770567}, {-0.143243713580954, 0.3834115661495822}},
        test_case{
            {3.0, 3.0}, {0.050366708890998826, 0.019519011699366543}, {-0.06986988715321024, -0.09836702149570346}
        }
    );

    std::complex<double> cpb[4] = {};
    std::complex<double> cpd[4] = {};
    xsf::specfun::cpbdn(-2, z, cpb, cpd);

    const double error = xsf::extended_relative_error(cpb[2], desired);
    const double error_d = xsf::extended_relative_error(cpd[2], desired_d);
    const double rtol = 1e-14;
    CAPTURE(z, cpb[2], desired, error, cpd[2], desired_d, error_d, rtol);
    REQUIRE(error <= rtol);
    REQUIRE(error_d <= rtol);
}
