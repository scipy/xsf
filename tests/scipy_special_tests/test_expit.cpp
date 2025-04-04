#include "../testing_utils.h"

#include <xsf/log_exp.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "expit"};

TEST_CASE("expit f->f scipy_special_tests", "[expit][f->f][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(xsf_test_cases<float, std::tuple<float, bool>, float>(
        tables_path / "In_f-f.parquet", tables_path / "Out_f-f.parquet",
        tables_path / ("Err_f-f_" + get_platform_str() + ".parquet")
    ));

    auto x = input;
    auto [desired, fallback] = output;
    auto out = xsf::expit(x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
