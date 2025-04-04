#include "../testing_utils.h"

#include <xsf/stats.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "gdtr"};

TEST_CASE("gdtr ddd->d scipy_special_tests", "[gdtr][ddd->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] =
        GENERATE(xsf_test_cases<std::tuple<double, double, double>, std::tuple<double, bool>, double>(
            tables_path / "In_d_d_d-d.parquet", tables_path / "Out_d_d_d-d.parquet",
            tables_path / ("Err_d_d_d-d_" + get_platform_str() + ".parquet")
        ));

    auto [a, b, x] = input;
    auto [desired, fallback] = output;
    auto out = xsf::gdtr(a, b, x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(a, b, x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
