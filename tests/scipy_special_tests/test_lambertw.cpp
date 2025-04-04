#include "../testing_utils.h"

#include <xsf/lambertw.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "lambertw"};

TEST_CASE("lambertw Dpd->D scipy_special_tests", "[lambertw][Dpd->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<
            std::tuple<std::complex<double>, std::ptrdiff_t, double>, std::tuple<std::complex<double>, bool>, double>(
            tables_path / "In_cd_p_d-cd.parquet", tables_path / "Out_cd_p_d-cd.parquet",
            tables_path / ("Err_cd_p_d-cd_" + get_platform_str() + ".parquet")
        )
    );

    auto [z, k, eps] = input;
    auto [desired, fallback] = output;
    auto out = xsf::lambertw(z, k, eps);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(z, k, eps, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
