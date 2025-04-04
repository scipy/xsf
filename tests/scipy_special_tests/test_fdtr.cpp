#include "../testing_utils.h"

#include <xsf/stats.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "fdtr"};

TEST_CASE("fdtr ddd->d scipy_special_tests", "[fdtr][ddd->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] =
        GENERATE(xsf_test_cases<std::tuple<double, double, double>, std::tuple<double, bool>, double>(
            tables_path / "In_d_d_d-d.parquet", tables_path / "Out_d_d_d-d.parquet",
            tables_path / ("Err_d_d_d-d_" + get_platform_str() + ".parquet")
        ));

    auto [dfn, dfd, x] = input;
    auto [desired, fallback] = output;
    auto out = xsf::fdtr(dfn, dfd, x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(dfn, dfd, x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
