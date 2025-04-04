#include "../testing_utils.h"

#include <xsf/specfun.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "hyp1f1"};

TEST_CASE("hyp1f1 ddD->D scipy_special_tests", "[hyp1f1][ddD->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] =
        GENERATE(xsf_test_cases<
                 std::tuple<double, double, std::complex<double>>, std::tuple<std::complex<double>, bool>, double>(
            tables_path / "In_d_d_cd-cd.parquet", tables_path / "Out_d_d_cd-cd.parquet",
            tables_path / ("Err_d_d_cd-cd_" + get_platform_str() + ".parquet")
        ));

    auto [a, b, z] = input;
    auto [desired, fallback] = output;
    auto out = xsf::hyp1f1(a, b, z);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(a, b, z, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
