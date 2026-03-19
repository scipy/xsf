#include "../testing_utils.h"
#include <xsf/stats.h>

/*
# Reference values computed with scipy.stats._hypotests._pval_cvm_2samp_exact

import numpy as np
from scipy.stats._stats import _kendall_dis

assert np.isclose(
    0,
    _kendall_dis(np.array([1, 2, 3, 4]), np.array([1, 2, 3, 4])),
)
assert np.isclose(
    6,
    _kendall_dis(np.array([1, 2, 3, 4]), np.array([4, 3, 2, 1])),
)
assert np.isclose(
    1,
    _kendall_dis(np.array([1, 2, 3, 4]), np.array([1, 3, 2, 4])),
)
*/
TEST_CASE("kendall_dis test", "[kendall_dis][xsf_tests]") {
    using test_case = std::tuple<std::vector<intptr_t>, std::vector<intptr_t>, int64_t>;
    auto [x, y, expected] = GENERATE(
        test_case{{1, 2, 3, 4}, {1, 2, 3, 4}, 0}, test_case{{1, 2, 3, 4}, {4, 3, 2, 1}, 6},
        test_case{{1, 2, 3, 4}, {1, 3, 2, 4}, 1}
    );
    const auto result = xsf::kendall_dis(x, y);
    CAPTURE(x, y, result, expected);
    REQUIRE(result == expected);
}
