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
