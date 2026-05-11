#include "../testing_utils.h"
#include <xsf/amos/amos.h>

TEST_CASE("amos besj", "[amos][xsf_tests]") {
    using std::complex;

    const complex<double> z{70.0, -0.7};
    const double fnu = 0.0;
    const int kode = 1;
    const int n = 15;

    std::array<complex<double>, n> cy{};
    int ierr = 0;

    const int nz = xsf::amos::besj(z, fnu, kode, n, cy.data(), &ierr);

    REQUIRE(ierr == 0);
    REQUIRE(nz == 0);

    // generated using:
    // from mpmath import mp
    // mp.dps = 50
    // n = 15
    // z = 70.0 - 0.7j
    // for n in range(n):
    //     result = mp.besselj(n, z)
    //     print(f"complex<double>{{{result.real:.20f}, {result.imag:.20f}}},")
    const std::array<complex<double>, n> ref = {
        complex<double>{0.11908800062250928549, 0.00765768917663659615},
        complex<double>{0.01289516192962112526, -0.07187585288729215033},
        complex<double>{-0.11869907035957860579, -0.00970739567078706507},
        complex<double>{-0.01967174120900008651, 0.07125337877045338558},
        complex<double>{0.11695202149331912143, 0.01579735764816333185},
        complex<double>{0.03301829774384231522, -0.06931450090245772537},
        complex<double>{-0.11213658264393736068, -0.02565127481716347501},
        complex<double>{-0.05219582054512202598, 0.06472536427721690399},
        complex<double>{0.10156902456216669875, 0.03849067209178914920},
        complex<double>{0.07532130636025452583, -0.05569624151237091633},
        complex<double>{-0.08205942069459905224, -0.05261746672772313695},
        complex<double>{-0.09861419022253167551, 0.04042975075154153429},
        complex<double>{0.05094243624931449501, 0.06501278997661191646},
        complex<double>{0.11585554323564695001, -0.01796723780085772449},
        complex<double>{-0.00784795132702038573, -0.07125539059627898735},
    };

    for (int i = 0; i < n; ++i) {
        const auto rel_error = xsf::extended_relative_error(cy[i], ref[i]);
        CAPTURE(i, cy[i], ref[i], rel_error);
        REQUIRE(rel_error <= 1e-14);
    }
}

TEST_CASE("amos besh vectorized", "[amos][xsf_tests]") {
    // tests the functionality of amos to return multiple consecutive orders for besh
    // by comparing to the versions returning only a single order
    using std::complex;

    using test_case = std::tuple<complex<double>, double, int, int, int, double>;
    auto [z, fnu, kode, m, n, rtol] = GENERATE(
        test_case{complex{14.0, -3.0}, 1.0, 1, 1, 260, 4e-13} // gh-92
    );

    std::vector<complex<double>> cy(n);
    int ierr = 0;
    int nz;

    nz = xsf::amos::besh(z, fnu, kode, m, n, cy.data(), &ierr);

    REQUIRE(ierr == 0);

    complex<double> ref;

    for (int i = 0; i < n; ++i) {
        nz = xsf::amos::besh(z, fnu + i, kode, m, 1, &ref, &ierr);
        REQUIRE(ierr == 0);

        const auto rel_error = xsf::extended_relative_error(cy[i], ref);
        CAPTURE(i, cy[i], ref, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
