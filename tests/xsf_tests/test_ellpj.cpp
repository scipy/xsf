#include "../testing_utils.h"
#include <iostream>
#include <tuple>
#include <xsf/cephes/ellpj.h>

TEST_CASE("ellpj negative m", "[ellpj][xsf_tests]") {
    // u_values = [-10, -6, -2, 2, 6, 10]
    // Reference values computed with mpmath

    using test_case = std::tuple<double, double, double, double, double>;
    // u, m, ref_sn, ref_cn, ref_dn

    auto [u, m, ref_sn, ref_cn, ref_dn] = GENERATE(
        // m = -2
        test_case{-10.0, -2.0, -0.649677834014239, 0.760209650024628, 1.3579994757064},
        test_case{-6.0, -2.0, -0.969931206391936, -0.243379240830198, 1.69750790580363},
        test_case{-2.0, -2.0, -0.348621324159513, -0.937263662125694, 1.11493213036376},
        test_case{2.0, -2.0, 0.348621324159513, -0.937263662125694, 1.11493213036376},
        test_case{6.0, -2.0, 0.969931206391936, -0.243379240830198, 1.69750790580363},
        test_case{10.0, -2.0, 0.649677834014239, 0.760209650024628, 1.3579994757064},

        // m = -5
        test_case{-10.0, -5.0, 0.495527059808153, -0.868592501117692, 1.49255999377263},
        test_case{-6.0, -5.0, 0.279115494901804, -0.960257538635194, 1.17878212468267},
        test_case{-2.0, -5.0, 0.0894599315921725, -0.995990421961739, 1.01981145159406},
        test_case{2.0, -5.0, -0.0894599315921725, -0.995990421961739, 1.01981145159406},
        test_case{6.0, -5.0, -0.279115494901804, -0.960257538635194, 1.17878212468267},
        test_case{10.0, -5.0, -0.495527059808153, -0.868592501117692, 1.49255999377263},

        // m = -10
        test_case{-10.0, -10.0, -0.673927384474905, 0.738797590991479, 2.35409880749553},
        test_case{-6.0, -10.0, 0.377239038384813, 0.926115925744991, 1.55662870358187},
        test_case{-2.0, -10.0, 0.518272042082259, -0.855215815099256, 1.91991122087485},
        test_case{2.0, -10.0, -0.518272042082259, -0.855215815099256, 1.91991122087485},
        test_case{6.0, -10.0, -0.377239038384813, 0.926115925744991, 1.55662870358187},
        test_case{10.0, -10.0, 0.673927384474905, 0.738797590991479, 2.35409880749553},

        // m = -0.5
        test_case{-10.0, -0.5, 0.993960849894643, 0.109735267242208, 1.22228437180619},
        test_case{-6.0, -0.5, -0.333655198670901, 0.942695183184831, 1.02745457116121},
        test_case{-2.0, -0.5, -0.766580871421719, -0.642147777050048, 1.13746345708987},
        test_case{2.0, -0.5, 0.766580871421719, -0.642147777050048, 1.13746345708987},
        test_case{6.0, -0.5, 0.333655198670901, 0.942695183184831, 1.02745457116121},
        test_case{10.0, -0.5, -0.993960849894643, 0.109735267242208, 1.22228437180619}
    );

    const double rtol = 1e-8;

    double sn, cn, dn, ph;
    int result = xsf::cephes::ellpj(u, m, &sn, &cn, &dn, &ph);

    REQUIRE(result == 0); // Should succeed

    const auto rel_error_sn = xsf::extended_relative_error(sn, ref_sn);
    const auto rel_error_cn = xsf::extended_relative_error(cn, ref_cn);
    const auto rel_error_dn = xsf::extended_relative_error(dn, ref_dn);
    const auto rel_error_ph = xsf::extended_relative_error(ph, std::atan2(ref_sn, ref_cn));

    CAPTURE(u, m, sn, ref_sn, rtol, rel_error_sn);
    CAPTURE(u, m, cn, ref_cn, rtol, rel_error_cn);
    CAPTURE(u, m, dn, ref_dn, rtol, rel_error_dn);

    REQUIRE(rel_error_sn <= rtol);
    REQUIRE(rel_error_cn <= rtol);
    REQUIRE(rel_error_dn <= rtol);
    REQUIRE(rel_error_ph <= rtol);
}
