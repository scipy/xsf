#include "../testing_utils.h"
#define MDSPAN_USE_PAREN_OPERATOR 1
#include <xsf/stats.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

/*
// Reference values computed with scipy.stats._hypotests._pval_cvm_2samp_exact

import numpy as np
from scipy import stats

rng = np.random.default_rng(seed=42)

list_m = rng.integers(3, 30, size=5)
list_n = rng.integers(3, 30, size=5)
rtol = 1e-10

for m, n in zip(list_m, list_n):
    x = rng.standard_normal(m)
    y = rng.standard_normal(n)
    res = stats.cramervonmises_2samp(x, y, method="exact")
    T = res.statistic
    # Convert normalized statistic T to the unnormalized U
    U = m * n * (m + n) * T + m * n * (4 * m * n - 1) / 6
    p_value = stats._hypotests._pval_cvm_2samp_exact(U, m, n)
    assert np.isclose(res.pvalue, p_value, rtol=rtol), "The p-values do not match!"
    print(f"U={U}, m={m}, n={n}, p-value={p_value}")
*/
TEST_CASE("pval_cvm_2samp_exact test", "[pval_cvm_2samp_exact][xsf_tests]") {
    using test_case = std::tuple<double, int, int, double, double>;
    auto [s, m, n, pval_expected, rtol] = GENERATE(
        test_case{8862.0, 14, 8, 0.2679738562091503, 1e-10},
        test_case{3491.0000000000005, 14, 5, 0.34657722738218094, 1e-10},
        test_case{12559.0, 5, 26, 0.11812654860485784, 1e-10},
        test_case{8901.0, 23, 5, 0.9907610907610908, 1e-10} //, test_case{119376.0, 20, 21, 0.5716351061359124, 1e-10}
    );

    const int64_t lcm = xsf::detail::lcm(m, n);
    const int64_t K = (m + n) * lcm * lcm + 1;

    std::vector<int64_t> buf1((m + 1) * K, 0);
    std::vector<int64_t> buf2((m + 1) * K, 0);

    using mdspan_2d = std::mdspan<int64_t, std::dextents<size_t, 2>>;
    mdspan_2d gs(buf1.data(), static_cast<size_t>(m + 1), static_cast<size_t>(K));
    mdspan_2d next_gs(buf2.data(), static_cast<size_t>(m + 1), static_cast<size_t>(K));

    auto pval = xsf::pval_cvm_2samp_exact(s, m, n, gs, next_gs);
    const double rel_error = xsf::extended_relative_error(pval, pval_expected);
    CAPTURE(s, m, n, K, pval, pval_expected, rel_error);
    REQUIRE(rel_error <= rtol);
}
