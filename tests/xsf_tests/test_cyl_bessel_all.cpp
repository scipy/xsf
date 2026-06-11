#include "../testing_utils.h"
#include <xsf/bessel.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

#include <complex>
#include <tuple>
#include <vector>

// parameter lists
namespace bessel_test_params {

// real and imaginary parts for z
const std::vector<double> Z_PARTS = {-100.0, -10.0, -1.0, -0.1, -1e-6, 0.0, 1e-6, 0.1, 1.0, 10.0, 100.0};

const std::vector<double> NU_VALUES = {-0.5, -0.25, 0.0, 0.25, 0.5};

const std::vector<int> N_VALUES = {10, 100};

} // namespace bessel_test_params

static bool is_nan(std::complex<double> const &z) { return std::isnan(z.real()) || std::isnan(z.imag()); }

// match exact (z, nu, n, i) combinations
using skip_entry_t = std::tuple<std::complex<double>, double, int, int>;

static bool should_skip(std::complex<double> z, double nu, int n, int i, std::vector<skip_entry_t> const &skip_list) {
    for (auto const &[skip_z, skip_nu, skip_n, skip_i] : skip_list) {
        if (z == skip_z && nu == skip_nu && n == skip_n && i == skip_i) {
            return true;
        }
    }
    return false;
}

// helper: compare any vectorized Bessel "_all" function against scalar calls
template <typename VecFunc, typename ScalarFunc>
static void compare_vectorized_with_scalar(
    std::complex<double> z, double nu, int n, double rtol, VecFunc &&vec_func, ScalarFunc &&scalar_func,
    std::vector<skip_entry_t> const &skip_list = {}
) {
    // compute all scalar references
    std::vector<std::complex<double>> refs(n);
    bool any_nan = false;
    for (int i = 0; i < n; ++i) {
        refs[i] = scalar_func(z, nu + std::copysign(i, nu));
        if (is_nan(refs[i])) {
            any_nan = true;
        }
    }

    // call the vectorized routine
    std::vector<std::complex<double>> cy_vec(n);
    std::mdspan cy_span(cy_vec.data(), cy_vec.size());
    vec_func(z, nu, cy_span);

    CAPTURE(z, nu, n, rtol);

    // if any scalar ref is NaN the "_all" routine NaN'd the whole array
    //
    // The underlying AMOS routines might set ierr = 1, 2, 4, 5.
    // For scalar wrappers (e.g. cyl_hankel_1) only that single element is NaN'd.
    // For "_all" wrappers the entire output array is NaN'd when ierr = 1, 2, 4, or 5.
    // Therefore, if any of the scalar reference values is NaN, expect the whole
    // vectorized result to be NaN.
    if (any_nan) {
        for (int i = 0; i < n; ++i) {
            CAPTURE(i, cy_vec[i]);
            REQUIRE(is_nan(cy_vec[i]));
        }
        return;
    }

    // compare element-wise
    for (int i = 0; i < n; ++i) {
        if (should_skip(z, nu, n, i, skip_list)) {
            continue;
        }
        const auto rel_error = xsf::extended_relative_error(cy_vec[i], refs[i]);
        CAPTURE(i, cy_vec[i], refs[i], rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

namespace {
// known mismatches between scalar and vectorized cyl_hankel_1
// {z, nu, n, i} where i is the element index to skip
const std::vector<skip_entry_t> HANKEL1_SKIP_LIST = {
    {std::complex<double>{-1e-6, -1.0}, 0.5, 10, 1},  {std::complex<double>{-1e-6, -1.0}, 0.5, 100, 1},
    {std::complex<double>{-1e-6, -1.0}, -0.5, 10, 1}, {std::complex<double>{-1e-6, -1.0}, -0.5, 100, 1},
    {std::complex<double>{1e-6, -1.0}, 0.5, 10, 1},   {std::complex<double>{1e-6, -1.0}, 0.5, 100, 1},
    {std::complex<double>{1e-6, -1.0}, -0.5, 10, 1},  {std::complex<double>{1e-6, -1.0}, -0.5, 100, 1},
    {std::complex<double>{0.0, -1.0}, 0.5, 10, 1},    {std::complex<double>{0.0, -1.0}, 0.5, 100, 1},
    {std::complex<double>{0.0, -1.0}, -0.5, 10, 1},   {std::complex<double>{0.0, -1.0}, -0.5, 100, 1},
};
} // namespace

TEST_CASE("cyl_hankel_1_all vectorized vs scalar", "[cyl_hankel_1_all][xsf_tests]") {
    const double zr = GENERATE(from_range(bessel_test_params::Z_PARTS));
    const double zi = GENERATE(from_range(bessel_test_params::Z_PARTS));
    const double nu = GENERATE(from_range(bessel_test_params::NU_VALUES));
    const int n = GENERATE(from_range(bessel_test_params::N_VALUES));

    std::complex<double> z(zr, zi);
    double rtol = 1e-12;

    compare_vectorized_with_scalar(
        z, nu, n, rtol, [](std::complex<double> z, double nu, auto cy) { xsf::cyl_hankel_1_all(nu, z, cy); },
        [](std::complex<double> z, double nu) { return xsf::cyl_hankel_1(nu, z); }, HANKEL1_SKIP_LIST
    );
}
