#include "../testing_utils.h"

#include <xsf/beta.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "beta"};

TEST_CASE("beta dd->d scipy_special_tests", "[beta][dd->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<std::tuple<double, double>, std::tuple<double, bool>, double>(
            tables_path / "In_d_d-d.parquet", tables_path / "Out_d_d-d.parquet",
            tables_path / ("Err_d_d-d_" + get_platform_str() + ".parquet")
        )
    );

    auto [a, b] = input;
    auto [desired, fallback] = output;
    auto out = xsf::beta(a, b);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(a, b, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
