/* These tests are adapted from the test code written by Stuart Brorson for
 * https://github.com/scipy/xsf/pull/99. It includes only tests for the
 * kernels featured in xsf for generating the recurrence matrices and
 * calculating Fourier sums. xsf does not contain full implementations of Mathieu
 * functions, only kernels which can potentially be used in SciPy or CuPy as
 * part of a full CPU or GPU implementation. Generating the Fourier coefficients
 * requires finding eigenvectors of a symmetric tridiagonal matrix, something which
 * is not currently available in xsf, and would likely require a different approach
 * on CPU vs GPU.
 *
 * Note that the SciPy implementations rely on
 * in-ufunc caching to allow for reuse of the same computed Fourier coefficients
 * across different values of the angle `x`. */

#include "../testing_utils.h"

#include <cmath>
#include <xsf/mathieu.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

auto constexpr Even = xsf::mathieu::Parity::Even;
auto constexpr Odd = xsf::mathieu::Parity::Odd;

TEST_CASE("make_matrix_ee", "[mathieu][xsf_tests]") {
    int N = 6;

    std::vector<double> D(N, 0.0);
    std::vector<double> E(N - 1, 0.0);
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
    std::vector<double> E(N - 1, 0.0);
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
    std::vector<double> E(N - 1, 0.0);
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

TEST_CASE("sum_fourier_series_even_even_q1_m0", "[mathieu][xsf_tests]") {
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

    /* Reference values computed with Wolfram Engine. */
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

TEST_CASE("sum_fourier_series_even_even_q1_m2", "[mathieu][xsf_tests]") {
    double q = 1.0;
    int m = 2;
    // Fourier coefficients
    std::vector<double> AA = {0.21692794675226462,     0.94825734682087837,     -0.081767008744172018,
                              0.0025865874071663215,   -4.3385838939211658e-05, 4.5372455146196267e-07,
                              -3.2496292879186506e-09, 1.6958298519778925e-11,  -6.7394972996424881e-14,
                              2.1085561518676298e-16,  -5.3296620523337814e-19, 1.1112098872910467e-21,
                              -1.9473431483638467e-24};

    std::mdspan AA_span(AA.data(), AA.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with Wolfram Engine. */
    std::vector<double> expected = {1.0859619368890114,   0.8856612185226889,  0.29865157290705935,
                                    -0.45171855334041594, -0.8157268390500982, -0.4517185533404161,
                                    0.2986515729070591,   0.8856612185226888,  1.0859619368890114};
    std::vector<double> expected_diff = {
        0.0,
        -1.0249411570023645,
        -1.8809997062068615,
        -1.6790771489646379,
        -3.182392406897388e-16,
        1.6790771489646377,
        1.8809997062068615,
        1.0249411570023648,
        3.1536401524284287e-16
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

TEST_CASE("sum_fourier_series_even_even_q1_m4", "[mathieu][xsf_tests]") {
    double q = 1.0;
    int m = 4;
    // Fourier coefficients
    std::vector<double> AA = {0.0052120044021059962,   0.083568404740583246,    0.99522416279527204,
                              -0.04989764213088127,    0.0010405258106442553,   -1.2393356664384822e-05,
                              9.6852903203680751e-08,  -5.3818530031774657e-10, 2.242785256765688e-12,
                              -7.2826318578220097e-15, 1.8966963764012226e-17,  -4.0530785034823585e-20,
                              7.2380970879747685e-23,  -1.0984019276257108e-25};

    std::mdspan AA_span(AA.data(), AA.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {1.0351351585780133,   0.09855498844065817, -0.9889717294331822,
                                    -0.09021203125324925, 0.9678184211483233,  -0.09021203125324875,
                                    -0.9889717294331822,  0.0985549884406577,  1.0351351585780133};
    std::vector<double> expected_diff = {
        0.0,
        -3.8874689332135333,
        -0.4663987362342739,
        4.074322044678967,
        1.0687171233875382e-15,
        -4.074322044678967,
        0.46639873623427197,
        3.8874689332135333,
        1.779033783159018e-15
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

TEST_CASE("sum_fourier_series_even_odd_q1_m1", "[mathieu][xsf_tests]") {
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
        /* i = 4 corresponds to v = M_PI / 2, a tricky case because it's near a root.
         * Even Wolfram warns of precision error here and the reference probably isn't
         * too accurate either. */
        double rtol = (i != 4) ? 1e-14 : 1e-6;
        // similarly, i = 8 corresponds to v = M_PI: near a root of the derivative.
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

TEST_CASE("sum_fourier_series_even_odd_q1_m3", "[mathieu][xsf_tests]") {
    double q = 1.0;
    int m = 3;
    // Fourier coefficients
    std::vector<double> AA = {0.1396156546048497,     0.98825110013683626,     -0.06216755513992854,
                              0.0015577824722684032,  -2.1662134262712558e-05, 1.935581567739718e-07,
                              -1.210366359154677e-09, 5.6056745161052482e-12,  -2.0026074434904442e-14,
                              5.6905300532462578e-17, -1.3174971993046381e-19, 2.5340382753735713e-22,
                              -4.114231210571251e-25, 5.7227812216111021e-28};

    std::mdspan AA_span(AA.data(), AA.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {1.067235512293139,   0.5295465974599358,     -0.5550307133804574,
                                    -0.9176193661208333, -1.927028889183907e-16, 0.9176193661208332,
                                    0.5550307133804577,  -0.5295465974599356,    -1.067235512293139};

    std::vector<double> expected_diff = {
        0.0,
        -2.509572450109097,
        -2.407068868211543,
        0.8763661055907831,
        3.1470770029784583,
        0.8763661055907842,
        -2.4070688682115424,
        -2.509572450109098,
        -9.25133290801537e-16
    };

    for (decltype(v.size()) i = 0; i < v.size(); i++) {
        /* i = 4 corresponds to v = M_PI / 2, a tricky case because it's near a root.
         * Even Wolfram warns of precision error here and the reference probably isn't
         * too accurate either. */
        double rtol = (i != 4) ? 1e-14 : 1e-5;
        // similarly, i = 8 corresponds to v = M_PI: near a root of the derivative.
        double rtol_diff = (i != 8) ? 1e-14 : 1e-5;
        double out, out_diff;
        xsf::mathieu::sum_fourier_series<Even, Odd>(AA_span, v[i], out, out_diff);
        double rel_error = xsf::extended_relative_error(out, expected[i]);
        double rel_error_diff = xsf::extended_relative_error(out_diff, expected_diff[i]);
        CAPTURE(m, q, i, out, out_diff, expected[i], expected_diff[i], rel_error, rel_error_diff);
        REQUIRE(rel_error <= rtol);
        REQUIRE(rel_error_diff <= rtol_diff);
    }
}

TEST_CASE("sum_fourier_series_odd_even_q1_m2", "[mathieu][xsf_tests]") {
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

TEST_CASE("sum_fourier_series_odd_even_q1_m4", "[mathieu][xsf_tests]") {
    double q = 1.0;
    int m = 4;
    // Fourier coefficients
    std::vector<double> BB = {
        0.082716278298929441,    0.9953225020162535,     -0.049900414382434623,   0.0010405649080668831,
        -1.2393695048768228e-05, 9.6854894971322646e-08, -5.3819378931634619e-10, 2.2428125738735586e-12,
        -7.2827001694483099e-15, 1.8967099081376478e-17, -4.0530999514897637e-20, 7.2381242151049455e-23,
        -1.0967424411688601e-25, 1.4299733865812567e-28,
    };

    std::mdspan BB_span(BB.data(), BB.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {
        0.0,
        1.018535489113575,
        0.1326042995245019,
        -0.9721093212091415,
        -2.525043858964328e-16,
        0.9721093212091416,
        -0.1326042995245014,
        -1.018535489113575,
        -4.721652450323951e-16
    };
    std::vector<double> expected_diff = {3.8555218154420925,   0.3204511230639429, -3.972966651023333,
                                         -0.33710016152124356, 4.123709563806253,  -0.3371001615212415,
                                         -3.972966651023333,   0.3204511230639411, 3.8555218154420925};

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

TEST_CASE("sum_fourier_series_odd_odd_q1_m1", "[mathieu][xsf_tests]") {
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
        // i = 8 corresponds to v = M_PI, a tricky case because it's near a root.
        double rtol = (i != 8) ? 1e-14 : 1e-6;
        // similarly, i = 4 corresponds to v = M_PI / 2: near a root of the derivative.
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

TEST_CASE("sum_fourier_series_odd_odd_q1_m3", "[mathieu][xsf_tests]") {
    double q = 1.0;
    int m = 3;
    // Fourier coefficients
    std::vector<double> BB = {0.10964247330123811,     0.9920165102302132,     -0.062284339384168448,
                              0.0015595116565660678,   -2.167694632566829e-05, 1.9363750953333752e-07,
                              -1.210630689061176e-09,  5.6061034523423122e-12, -2.0025415533333008e-14,
                              5.6898475991588914e-17,  -1.31724578162025e-19,  2.5334054567707489e-22,
                              -4.1129992356066497e-25, 5.7208238435360791e-28};

    std::mdspan BB_span(BB.data(), BB.size());

    const std::vector<double> v = linspace(0.0, M_PI, 9);

    /* Reference values computed with wolfram engine. */
    std::vector<double> expected = {
        0.0,
        0.9015237982782264,
        0.8219142851414616,
        -0.2530357599220488,
        -0.9462397597698017,
        -0.2530357599220492,
        0.8219142851414614,
        0.9015237982782265,
        3.410633878148768e-16
    };

    std::vector<double> expected_diff = {2.784993910508225,   1.3494510247766995,     -1.7990677315460222,
                                         -2.9993699766876163, -6.401112565382882e-16, 2.999369976687616,
                                         1.799067731546023,   -1.3494510247766986,    -2.784993910508225};

    for (decltype(v.size()) i = 0; i < v.size(); i++) {
        // i = 8 corresponds to v = M_PI, close to a root.
        double rtol = (i != 8) ? 1e-14 : 1e-5;
        // i = 4 corresponds to v = M_PI / 2, close to a root of the derivative.
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
