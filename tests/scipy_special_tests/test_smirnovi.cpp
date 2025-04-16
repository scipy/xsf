#include "../testing_utils.h"

#include <xsf/stats.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "smirnovi"};


TEST_CASE("smirnovi pd->d scipy_special_tests", "[smirnovi][pd->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] =
        GENERATE(xsf_test_cases<std::tuple<std::ptrdiff_t, double>, std::tuple<double, bool>, double>(
            tables_path / "In_p_d-d.parquet", tables_path / "Out_p_d-d.parquet",
            tables_path / ("Err_p_d-d_" + get_platform_str() + ".parquet")
        ));

    auto [n, p] = input;
    auto [desired, fallback] = output;
    auto out = xsf::smirnovi(n, p);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(n, p, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
