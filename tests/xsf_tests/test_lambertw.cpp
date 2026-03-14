#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/lambertw.h>

TEST_CASE("lambertw", "[lambertw][xsf_tests]") {
    using test_case = std::tuple<std::complex<double>, int, std::complex<double>>;
    using std::complex;
    // Reference values were computed with the Python library mpmath.
    auto [z, k, ref] = GENERATE(
        test_case{complex{0.0, 0.0}, 0, complex{0.0, 0.0}},
        test_case{complex{1.0, 0.0}, 0, complex{0.56714329040978387, 0.0}},
        test_case{complex{M_E, 0.0}, 0, complex{1.0, 0.0}},
        test_case{complex{-M_PI/2, 0.0}, 0, complex{0.0, M_PI/2}},
        test_case{complex{-M_LN2/2, 0.0}, 0, complex{-M_LN2, 0.0}},
        test_case{complex{-0.0001, 0.0}, 0, complex{-0.00010001000150026672, 0.0}},
        test_case{complex{-1.0/M_E, 0.0}, 0, complex{-1.0, 8.220079714836618e-9}},
        test_case{complex{-0.36787944, 0.0}, 0, complex{-0.9999201984841515, 0.0}},
        test_case{complex{-0.36787945, 0.0}, 0, complex{-0.9999999840009949, 0.00021908221081537528}},
        test_case{complex{0.25, 0.0}, -1, complex{-3.008998009970046, -4.076529788991597}},
        test_case{complex{-0.25, 0.0}, -1, complex{-2.15329236411035, 0.0}},
        test_case{complex{-1.0/M_E, 0.0}, -1, complex{-1.0, -8.220079714836618e-9}},
        test_case{complex{-0.36787945, 0.0}, -1, complex{-0.9999999840009949, -0.00021908221081537528}},
        test_case{complex{-0.36787944, 0.0}, -1, complex{-1.0000798057615958, 0.0}},
        test_case{complex{2.0, 3.0}, -2, complex{-1.0163597810788894, -9.910584876118309}},
        test_case{complex{2.0, 3.0}, -1, complex{-0.0315828083898750, -3.721107986637061}},
        test_case{complex{2.0, 3.0},  0, complex{ 1.0900765344857908,  0.530139720774839}},
        test_case{complex{2.0, 3.0},  1, complex{-0.4462717121285748,  5.615883357759838}},
        test_case{complex{2.0, 3.0},  2, complex{-1.1972601078812846, 11.877910111617442}}
    );

    const double rtol = 1e-9;
    const complex w = xsf::lambertw(z, k, rtol);
    const auto rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(z, k, w, ref, rel_error);
    REQUIRE(rel_error <= rtol);
}
