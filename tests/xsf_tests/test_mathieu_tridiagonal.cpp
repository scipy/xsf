#include "../testing_utils.h"

#include <xsf/mathieu_very_new.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

#include <cstdio>
#include <stdlib.h>

/*
 *
 * The goal of these tests are to verify that my C/C++
 * impl of the Mathieu fcns has been carried over from
 * my Matlab impl correctly.  Therefore, main just calls
 * a bunch of golden value tests.  The GVs were generated
 * using the Matlab impl.  I did fairly extensive
 * validation of the Matlab impl.  Therefore, if the tests
 * here show the C impl matches the Matlab impl, then that
 * should serve as verification of the C impl's correctness.
 *
 * A secondary goal is to show how to call the various fcns
 * in my API.
 *
 */

TEST_CASE("make_matrix_ee", "[mathieu_tridiagonal][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.8284271247461903, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {0.0, 4.0, 16.0, 36.0, 64.0, 100.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix_ee(q, D_span, E_span);

    for (int i = 0; i < N - 1; ++i) {
        const double rel_error = xsf::extended_relative_error(E[i], E_expected[i]);
        CAPTURE(i, E[i], E_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    for (int i = 0; i < N; ++i) {
        const double rel_error = xsf::extended_relative_error(D[i], D_expected[i]);
        CAPTURE(i, D[i], D_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("make_matrix_eo", "[mathieu_tridiagonal][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {3.0, 9.0, 25.0, 49.0, 81.0, 121.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix_eo(q, D_span, E_span);

    for (int i = 0; i < N - 1; ++i) {
        const double rel_error = xsf::extended_relative_error(E[i], E_expected[i]);
        CAPTURE(i, E[i], E_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    for (int i = 0; i < N; ++i) {
        const double rel_error = xsf::extended_relative_error(D[i], D_expected[i]);
        CAPTURE(i, D[i], D_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("make_matrix_oe", "[mathieu_tridiagonal][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {4.0, 16.0, 36.0, 64.0, 100.0, 144.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix_oe(q, D_span, E_span);

    for (int i = 0; i < N - 1; ++i) {
        const double rel_error = xsf::extended_relative_error(E[i], E_expected[i]);
        CAPTURE(i, E[i], E_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    for (int i = 0; i < N; ++i) {
        const double rel_error = xsf::extended_relative_error(D[i], D_expected[i]);
        CAPTURE(i, D[i], D_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("make_matrix_oo", "[mathieu_tridiagonal][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N - 1, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {-1.0, 9.0, 25.0, 49.0, 81.0, 121.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix_oo(q, D_span, E_span);

    for (int i = 0; i < N - 1; ++i) {
        const double rel_error = xsf::extended_relative_error(E[i], E_expected[i]);
        CAPTURE(i, E[i], E_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    for (int i = 0; i < N; ++i) {
        const double rel_error = xsf::extended_relative_error(D[i], D_expected[i]);
        CAPTURE(i, D[i], D_expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
