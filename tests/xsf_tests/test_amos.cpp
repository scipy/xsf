#include "../testing_utils.h"
#include <xsf/amos/amos.h>

TEST_CASE("amos besj vectorized", "[amos][xsf_tests]") {
    // tests the functionality of amos to return multiple consecutive orders for besj
    // by comparing to the versions returning only a single order
    using std::complex;

    using test_case = std::tuple<complex<double>, double, int, int, double>;
    auto [z, fnu, kode, n, rtol] = GENERATE(
        test_case{complex{70.0, -0.7}, 0.0, 1, 15, 5e-13},    // gh-138
        test_case{complex{-50.0, -50.0}, 0.0, 1, 1000, 5e-13} // gh-145
    );

    std::vector<complex<double>> cy(n);
    int ierr = 0;
    int nz;

    nz = xsf::amos::besj(z, fnu, kode, n, cy.data(), &ierr);

    REQUIRE(ierr == 0);

    complex<double> ref;

    for (int i = 0; i < n; ++i) {
        nz = xsf::amos::besj(z, fnu + i, kode, 1, &ref, &ierr);
        REQUIRE(ierr == 0);

        const auto rel_error = xsf::extended_relative_error(cy[i], ref);
        CAPTURE(i, cy[i], ref, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("amos besi vectorized", "[amos][xsf_tests]") {
    // tests the functionality of amos to return multiple consecutive orders for besi
    // by comparing to the versions returning only a single order
    using std::complex;

    using test_case = std::tuple<complex<double>, double, int, int, double>;
    auto [z, fnu, kode, n, rtol] = GENERATE(
        test_case{complex{0.0, -15.0}, 0.0, 1, 10, 1e-15} // gh-158
    );

    std::vector<complex<double>> cy(n);
    int ierr = 0;
    int nz;

    nz = xsf::amos::besi(z, fnu, kode, n, cy.data(), &ierr);

    REQUIRE(ierr == 0);

    complex<double> ref;

    for (int i = 0; i < n; ++i) {
        nz = xsf::amos::besi(z, fnu + i, kode, 1, &ref, &ierr);
        REQUIRE(ierr == 0);

        const auto rel_error = xsf::extended_relative_error(cy[i], ref);
        CAPTURE(i, cy[i], ref, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("amos seri buffer overflow gh-92", "[amos][xsf_tests]") {
    using std::complex;

    // parameters for besh which trigger overflow in seri
    const complex<double> z{14.0, -3.0};
    const double fnu = 1.0;
    const int m = 1;
    const int n = 260;
    const int kode = 1;

    const complex<double> sentinel{12345.67, 98765.43};
    std::vector<complex<double>> cy(n + 2, sentinel);

    int ierr = 0;
    int nz = xsf::amos::besh(z, fnu, kode, m, n, cy.data() + 1, &ierr);

    // check if guards elements were touched
    CAPTURE(cy[0], cy[n + 1]);
    CHECK(cy[0] == sentinel);
    CHECK(cy[n + 1] == sentinel);
}

TEST_CASE("amos asyi buffer overflow gh-158", "[amos][xsf_tests]") {
    using std::complex;

    // parameters for asyi (computed as in besi)
    const double tol = 2.220446049250313e-16;
    const double rl = 24.6;
    const double elim = 700.0;
    const double alim = 664.0;

    const complex<double> z{700.0, 10.0};
    const double fnu = 0.0;
    const int n = 5;
    const int kode = 1;

    // allocate n+1 elements, initialize extra to sentinel to detect overflow
    const complex<double> sentinel{12345.67, 98765.43};
    std::vector<complex<double>> y(n + 1, sentinel);

    int nz = xsf::amos::asyi(z, fnu, kode, n, y.data(), rl, tol, elim, alim);

    REQUIRE(nz == 0);

    // check if the extra element (index n) was touched
    CAPTURE(y[n]);
    CHECK(y[n] == sentinel);
}
