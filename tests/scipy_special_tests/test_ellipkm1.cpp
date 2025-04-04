#include "../testing_utils.h"

#include <xsf/ellip.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "ellipkm1"};

TEST_CASE("ellipkm1 d->d scipy_special_tests", "[ellipkm1][d->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(xsf_test_cases<double, std::tuple<double, bool>, double>(
        tables_path / "In_d-d.parquet", tables_path / "Out_d-d.parquet",
        tables_path / ("Err_d-d_" + get_platform_str() + ".parquet")
    ));

    auto p = input;
    auto [desired, fallback] = output;
    auto out = xsf::ellipkm1(p);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(p, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
