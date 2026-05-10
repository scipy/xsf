#include "../testing_utils.h"

#include <cmath>
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

TEST_CASE("sum_fourier_series_even_even", "[mathieu][xsf_texts]") {
    double q = 1.0;
    int m = 0;
    // Fourier coefficients
    std::vector<double> AA = {
        0.67298967231649987,    -0.30630358003683739,    0.018645559365419482,   -0.00051168367225324469,
        7.9398280312329386e-06, -7.9043992651665076e-08, 5.4720641212326225e-10, -2.7854566570351717e-12,
        1.0861510315639324e-14, -3.3476412051017353e-17, 8.3596341717010846e-20, -1.7255806451365471e-22,
        2.9934420820894599e-25, -4.4251976995614426e-28, 5.6411179512232202e-31, -6.264747486668504e-34,
        6.1152045985104173e-37, -5.2878909636082434e-40, 4.0787325607116287e-43, -2.8237182422819167e-46,
        1.7643226442884998e-49, -9.9992520186900771e-53, 5.1636902188609567e-56, -2.4397830047988059e-59,
        1.0587244524433085e-62,
    };

    std::mdspan AA_span(AA.data(), AA.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {0.38482782930129905, 0.45675426482776166, 0.6543520522319161,
                                    0.8892092001491974,  0.9984585148130367,  0.8892092001491975,
                                    0.6543520522319162,  0.4567542648277617,  0.38482782930129905};
    std::vector<double> expected_diff = {
        0.0,
        0.36076677721015715,
        0.609537848441086,
        0.5099312390005591,
        9.444966065193975e-17,
        -0.509931239000559,
        -0.6095378484410862,
        -0.36076677721015726,
        -1.1570532269097101e-16,
    };

    for (decltype(v.size()) i = 0; i < v.size(); i++) {
        double out, out_diff;
        xsf::mathieu::sum_fourier_series<Even, Even>(AA_span, v[i], out, out_diff);
        double rel_error = xsf::extended_relative_error(out, expected[i]);
        double rel_error_diff = xsf::extended_relative_error(out_diff, expected_diff[i]);
        CAPTURE(m, q, out, out_diff, expected[i], expected_diff[i], rel_error, rel_error_diff);
        REQUIRE(rel_error <= 1e-15);
        REQUIRE(rel_error_diff <= 1e-15);
    }
}

TEST_CASE("sum_fourier_series_even_odd", "[mathieu][xsf_texts]") {
    double q = 1.0;
    int m = 1;
    // Fourier coefficients
    std::vector<double> AA = {
        0.99020205940794281,    -0.13951147675023179,    0.0060343187093877355,  -0.00012804035971445581,
        1.6180502678122694e-06, -1.358166371415527e-08,  8.1260951356738996e-11, -3.6417448610665401e-13,
        1.2682903335036881e-15, -3.5314786030658603e-18, 8.0418242153716169e-21, -1.5255596482175257e-23,
        2.4521128059772998e-26,
    };

    std::mdspan AA_span(AA.data(), AA.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {0.8565984655568863, 0.8592462553771193,   0.794471810010097,
                                    0.5134490713977,    8.81634913001536e-17, -0.5134490713976999,
                                    -0.794471810010097, -0.8592462553771193,  -0.8565984655568863};

    std::vector<double> expected_diff = {
        0.0,
        -0.01978501858311021,
        -0.38353947703721303,
        -1.062605871301982,
        -1.439819078636166,
        -1.0626058713019821,
        -0.3835394770372132,
        -0.019785018583110256,
        1.4779993885759145e-17
    };

    for (decltype(v.size()) i = 0; i < v.size(); i++) {
        // i = 4 corresponds to v near pi / 2, a tricky case because pi/2 is a root.
        double rtol = (i != 4) ? 1e-14 : 1e-6;
        // similarly, i = 8 corresponds to v near pi: a root of the derivative.
        double rtol_diff = (i != 8) ? 1e-14 : 5e-4;
        double out, out_diff;
        xsf::mathieu::sum_fourier_series<Even, Odd>(AA_span, v[i], out, out_diff);
        double rel_error = xsf::extended_relative_error(out, expected[i]);
        double rel_error_diff = xsf::extended_relative_error(out_diff, expected_diff[i]);
        CAPTURE(m, q, i, out, out_diff, expected[i], expected_diff[i], rel_error, rel_error_diff);
        REQUIRE(rel_error <= rtol);
        REQUIRE(rel_error_diff <= rtol_diff);
    }
}

TEST_CASE("sum_fourier_series_odd_even", "[mathieu][xsf_texts]") {
    double q = 1.0;
    int m = 2;
    // Fourier coefficients
    std::vector<double> BB = {
        0.99657191561800729,    -0.082690780921752541,   0.0025787417609396878,  -4.2927111309754299e-05,
        4.4680445101583735e-07, -3.1896885143918772e-09, 1.6606127571064949e-11, -6.5876457417632196e-14,
        2.0581218664476538e-16, -5.1962159920821663e-19, 1.0823618374177297e-21, -1.8919714509406692e-24,
        2.8192809544910441e-27,
    };

    std::mdspan BB_span(BB.data(), BB.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {
        0.0,
        0.6238151115773868,
        0.9939936206449127,
        0.7891966670415149,
        1.4326697921279722e-16,
        -0.7891966670415148,
        -0.9939936206449127,
        -0.6238151115773869,
        -2.0543623603136102e-16
    };
    std::vector<double> expected_diff = {1.6775141712238477,  1.3987651018998324, 0.33041974507174054,
                                         -1.398078268120984,  -2.339727328933457, -1.3980782681209847,
                                         0.33041974507174005, 1.3987651018998322, 1.6775141712238477};

    for (decltype(v.size()) i = 0; i < v.size(); i++) {
        double out, out_diff;
        xsf::mathieu::sum_fourier_series<Odd, Even>(BB_span, v[i], out, out_diff);
        double rel_error = xsf::extended_relative_error(out, expected[i]);
        double rel_error_diff = xsf::extended_relative_error(out_diff, expected_diff[i]);
        CAPTURE(m, q, out, out_diff, expected[i], expected_diff[i], rel_error, rel_error_diff);
        REQUIRE(rel_error <= 1e-15);
        REQUIRE(rel_error_diff <= 1e-15);
    }
}

TEST_CASE("sum_fourier_series_odd_odd", "[mathieu][xsf_texts]") {
    double q = 1.0;
    int m = 1;
    // Fourier coefficients
    std::vector<double> BB = {0.99396796139893406,     -0.109583791872277,      0.0043676488669854262,
                              -8.8957922980728422e-05, 1.0968648409691742e-06,  -9.0571888918196522e-09,
                              5.3559301601145843e-11,  -2.3792841444697411e-13, 8.2297564493057137e-16,
                              -2.2790292278436223e-18, 5.1665972365745621e-21,  -9.764718967266412e-24,
                              1.5645823771962644e-26};

    std::mdspan BB_span(BB.data(), BB.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {
        0.0,
        0.28313357595976163,
        0.6223293212236893,
        0.9584879270333326,
        1.1080094660370052,
        0.9584879270333327,
        0.6223293212236894,
        0.28313357595976174,
        8.4064887257814e-17
    };

    std::vector<double> expected_diff = {0.686441897503372,   0.784708039274585,     0.9194285732565753,
                                         0.7045120739186153,  1.282120820069729e-16, -0.7045120739186151,
                                         -0.9194285732565753, -0.784708039274585,    -0.686441897503372};

    for (decltype(v.size()) i = 0; i < v.size(); i++) {
        // i = 8 corresponds to v near pi, a tricky case because pi is a root.
        double rtol = (i != 8) ? 1e-14 : 1e-6;
        // similarly, i = 4 corresponds to v near pi/2: a root of the derivative.
        double rtol_diff = (i != 4) ? 1e-14 : 1e-5;
        double out, out_diff;
        xsf::mathieu::sum_fourier_series<Odd, Odd>(BB_span, v[i], out, out_diff);
        double rel_error = xsf::extended_relative_error(out, expected[i]);
        double rel_error_diff = xsf::extended_relative_error(out_diff, expected_diff[i]);
        CAPTURE(m, q, i, out, out_diff, expected[i], expected_diff[i], rel_error, rel_error_diff);
        REQUIRE(rel_error <= rtol);
        REQUIRE(rel_error_diff <= rtol_diff);
    }
}
