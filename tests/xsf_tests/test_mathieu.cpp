#include "../testing_utils.h"

#include <xsf/mathieu.h>
#include <xsf/third_party/kokkos/mdspan.hpp>


auto constexpr Even = xsf::mathieu::Parity::Even;
auto constexpr Odd = xsf::mathieu::Parity::Odd;


TEST_CASE("make_matrix_ee", "[mathieu][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.8284271247461903, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {0.0, 4.0, 16.0, 36.0, 64.0, 100.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix<Even, Even>(q, D_span, E_span);

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

TEST_CASE("make_matrix_eo", "[mathieu][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {3.0, 9.0, 25.0, 49.0, 81.0, 121.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix<Even, Odd>(q, D_span, E_span);

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

TEST_CASE("make_matrix_oe", "[mathieu][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {4.0, 16.0, 36.0, 64.0, 100.0, 144.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix<Odd, Even>(q, D_span, E_span);

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

TEST_CASE("make_matrix_oo", "[mathieu][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N - 1, 0.0);
    double q = 2.0;

    std::vector<double> E_expected = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> D_expected = {-1.0, 9.0, 25.0, 49.0, 81.0, 121.0};
    double rtol = 1e-15;

    std::mdspan E_span(E.data(), E.size());
    std::mdspan D_span(D.data(), D.size());

    xsf::mathieu::make_matrix<Odd, Odd>(q, D_span, E_span);

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
