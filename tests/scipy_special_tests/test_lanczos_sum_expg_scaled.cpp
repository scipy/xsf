#include "../testing_utils.h"

#include <xsf/cephes/lanczos.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "lanczos_sum_expg_scaled"};

TEST_CASE(
    "cephes::lanczos_sum_expg_scaled d->d scipy_special_tests",
    "[cephes::lanczos_sum_expg_scaled][d->d][scipy_special_tests]"
) {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(xsf_test_cases<double, std::tuple<double, bool>, double>(
        tables_path / "In_d-d.parquet", tables_path / "Out_d-d.parquet",
        tables_path / ("Err_d-d_" + get_platform_str() + ".parquet")
    ));

    auto z = input;
    auto [desired, fallback] = output;
    auto out = xsf::cephes::lanczos_sum_expg_scaled(z);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(z, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
