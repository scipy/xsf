#include "../testing_utils.h"

#include <xsf/gamma.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "gamma"};

TEST_CASE("gamma D->D scipy_special_tests", "[gamma][D->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<std::complex<double>, std::tuple<std::complex<double>, bool>, double>(
            tables_path / "In_cd-cd.parquet", tables_path / "Out_cd-cd.parquet",
            tables_path / ("Err_cd-cd_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired, fallback] = output;
    auto out = xsf::gamma(x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}

TEST_CASE("gamma D->D huge negative inputs underflow to zero", "[gamma][D->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto z = GENERATE(
        std::complex<double>{-8.98847e307, -10.0}, std::complex<double>{-8.98847e307, 1e-30},
        std::complex<double>{-8.98847e307, 2.99808e154}
    );

    auto out = xsf::gamma(z);
    CAPTURE(z, out);
    REQUIRE(out.real() == 0.0);
    REQUIRE(out.imag() == 0.0);
}

TEST_CASE("gamma D->D positive real axis overflows", "[gamma][D->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto z = GENERATE(std::complex<double>{2.99808e154, 0.0}, std::complex<double>{8.98847e307, 0.0});

    auto out = xsf::gamma(z);
    CAPTURE(z);
    CAPTURE(out);
    REQUIRE(out.real() == std::numeric_limits<double>::infinity());
    REQUIRE(out.imag() == 0.0);
}

TEST_CASE("gamma D->D balanced huge inputs overflow", "[gamma][D->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto z = GENERATE(std::complex<double>{2.99808e154, 2.99808e154}, std::complex<double>{2.99808e154, -2.99808e154});

    auto out = xsf::gamma(z);
    CAPTURE(z, out);
    REQUIRE(out.real() == -std::numeric_limits<double>::infinity());
    REQUIRE(out.imag() == std::copysign(std::numeric_limits<double>::infinity(), z.imag()));
}

TEST_CASE("gamma d->d scipy_special_tests", "[gamma][d->d][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<double, std::tuple<double, bool>, double>(
            tables_path / "In_d-d.parquet", tables_path / "Out_d-d.parquet",
            tables_path / ("Err_d-d_" + get_platform_str() + ".parquet")
        )
    );

    auto x = input;
    auto [desired, fallback] = output;
    auto out = xsf::gamma(x);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(x, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
