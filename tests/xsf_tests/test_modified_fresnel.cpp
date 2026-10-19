#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/cpu/fresnel.h>

TEST_CASE("modified fresnel negative x", "[modified_fresnel][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double, double>;
    // Reference values of F+ and K+ were computed with the Python library mpmath.
    // F- and K- are their complex conjugates.
    auto [x, ref_f_re, ref_f_im, ref_k_re, ref_k_im] = GENERATE(
        test_case{-5.0, 1.2381238350542128, 1.1545743498230725, 0.9505625997663579, 0.09329796120960479},
        test_case{-2.0, 1.0881185310909665, 1.4314335580015063, -0.760667562401997, 0.6711795473720403},
        test_case{-1.0, 1.5311813065580222, 0.9369253703811312, 0.3325082579142154, -0.9566308785853518},
        test_case{-0.3, 0.9264141597657186, 0.6356518628567259, 0.6102263532447829, -0.17153789274616973}
    );
    constexpr double rtol = 1e-14;
    const std::complex<double> ref_fp(ref_f_re, ref_f_im), ref_kp(ref_k_re, ref_k_im);
    std::complex<double> fp, kp, fm, km;
    xsf::modified_fresnel_plus(x, fp, kp);
    xsf::modified_fresnel_minus(x, fm, km);

    const auto rel_error_fp = xsf::extended_relative_error(fp, ref_fp);
    const auto rel_error_kp = xsf::extended_relative_error(kp, ref_kp);
    const auto rel_error_fm = xsf::extended_relative_error(fm, std::conj(ref_fp));
    const auto rel_error_km = xsf::extended_relative_error(km, std::conj(ref_kp));

    CAPTURE(x, fp, ref_fp, rel_error_fp, kp, ref_kp, rel_error_kp);
    CAPTURE(fm, rel_error_fm, km, rel_error_km);
    REQUIRE(rel_error_fp <= rtol);
    REQUIRE(rel_error_kp <= rtol);
    REQUIRE(rel_error_fm <= rtol);
    REQUIRE(rel_error_km <= rtol);
}
