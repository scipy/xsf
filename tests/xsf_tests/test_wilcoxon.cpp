#include "../testing_utils.h"
#include <xsf/stats.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

namespace {

std::vector<double> wilcoxon_pmf_table(int n) {
    std::vector<double> table(n * (n + 1) / 2 + 1, 0.0);
    std::mdspan table_span(table.data(), table.size());
    xsf::wilcoxon_pmf_all(n, table_span);
    return table;
}

std::vector<double> wilcoxon_cdf_table(int n) {
    std::vector<double> table(n * (n + 1) / 2 + 1, 0.0);
    std::mdspan table_span(table.data(), table.size());
    xsf::wilcoxon_cdf_all(n, table_span);
    return table;
}

} // namespace

// Values computed from scipy.stats._wilcoxon
//
// # Example for CDF and n = 5
// from scipy.stats._wilcoxon import WilcoxonDistribution
// import numpy as np
// np.set_printoptions(precision=17)
// n = 5
// wd = WilcoxonDistribution(n)
// wd.cdf([5, 6, 7, 8, 9])

TEST_CASE("wilcoxon_cdf exact small cases", "[wilcoxon_cdf]") {
    {
        // n = 0
        std::vector<double> cdf_table = wilcoxon_cdf_table(0);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_cdf(0.0, cdf_span) == 1.0);
    }

    { // n = 1
        std::vector<double> cdf_table = wilcoxon_cdf_table(1);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_cdf(0.0, cdf_span) == 0.5);
        REQUIRE(xsf::wilcoxon_cdf(1.0, cdf_span) == 1.0);
    }

    { // n = 2
        std::vector<double> cdf_table = wilcoxon_cdf_table(2);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_cdf(1.0, cdf_span) == 0.5);
        REQUIRE(xsf::wilcoxon_cdf(2.0, cdf_span) == 0.75);
        REQUIRE(xsf::wilcoxon_cdf(3.0, cdf_span) == 1.0);
        REQUIRE(xsf::wilcoxon_cdf(4.0, cdf_span) == 1.0);
    }

    { // n = 3
        std::vector<double> cdf_table = wilcoxon_cdf_table(3);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_cdf(3.0, cdf_span) == 0.625);
        REQUIRE(xsf::wilcoxon_cdf(4.0, cdf_span) == 0.75);
        REQUIRE(xsf::wilcoxon_cdf(5.0, cdf_span) == 0.875);
        REQUIRE(xsf::wilcoxon_cdf(6.0, cdf_span) == 1.0);
        REQUIRE(xsf::wilcoxon_cdf(7.0, cdf_span) == 1.0);
    }

    { // n = 5
        std::vector<double> cdf_table = wilcoxon_cdf_table(5);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_cdf(5.0, cdf_span) == 0.3125);
        REQUIRE(xsf::wilcoxon_cdf(6.0, cdf_span) == 0.40625);
        REQUIRE(xsf::wilcoxon_cdf(7.0, cdf_span) == 0.5);
        REQUIRE(xsf::wilcoxon_cdf(8.0, cdf_span) == 0.59375);
        REQUIRE(xsf::wilcoxon_cdf(9.0, cdf_span) == 0.6875);
    }
}

TEST_CASE("wilcoxon_cdf edge cases", "[wilcoxon_cdf]") {
    std::vector<double> cdf_table = wilcoxon_cdf_table(2);
    std::mdspan cdf_span(cdf_table.data(), cdf_table.size());

    REQUIRE(xsf::wilcoxon_cdf(4.0, cdf_span) == 1.0);

    REQUIRE(xsf::wilcoxon_cdf(1.9, cdf_span) == xsf::wilcoxon_cdf(1.0, cdf_span));
    REQUIRE(xsf::wilcoxon_cdf(-0.5, cdf_span) == xsf::wilcoxon_cdf(0.0, cdf_span));

    REQUIRE(std::isnan(xsf::wilcoxon_cdf(std::numeric_limits<double>::quiet_NaN(), cdf_span)));
}

TEST_CASE("wilcoxon_sf exact small cases", "[wilcoxon_sf]") {
    {
        // n = 0
        std::vector<double> cdf_table = wilcoxon_cdf_table(0);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_sf(0.0, cdf_span) == 1.0);
        REQUIRE(xsf::wilcoxon_sf(1.0, cdf_span) == 0.0);
    }

    { // n = 1
        std::vector<double> cdf_table = wilcoxon_cdf_table(1);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_sf(0.0, cdf_span) == 1.0);
        REQUIRE(xsf::wilcoxon_sf(1.0, cdf_span) == 0.5);
        REQUIRE(xsf::wilcoxon_sf(2.0, cdf_span) == 0.0);
    }

    { // n = 2
        std::vector<double> cdf_table = wilcoxon_cdf_table(2);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_sf(0.0, cdf_span) == 1.0);
        REQUIRE(xsf::wilcoxon_sf(1.0, cdf_span) == 0.75);
        REQUIRE(xsf::wilcoxon_sf(2.0, cdf_span) == 0.5);
        REQUIRE(xsf::wilcoxon_sf(3.0, cdf_span) == 0.25);
        REQUIRE(xsf::wilcoxon_sf(4.0, cdf_span) == 0.0);
    }

    { // n = 3
        std::vector<double> cdf_table = wilcoxon_cdf_table(3);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_sf(3.0, cdf_span) == 0.625);
        REQUIRE(xsf::wilcoxon_sf(4.0, cdf_span) == 0.375);
        REQUIRE(xsf::wilcoxon_sf(5.0, cdf_span) == 0.25);
        REQUIRE(xsf::wilcoxon_sf(6.0, cdf_span) == 0.125);
        REQUIRE(xsf::wilcoxon_sf(7.0, cdf_span) == 0.0);
    }

    { // n = 5
        std::vector<double> cdf_table = wilcoxon_cdf_table(5);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::wilcoxon_sf(5.0, cdf_span) == 0.78125);
        REQUIRE(xsf::wilcoxon_sf(6.0, cdf_span) == 0.6875);
        REQUIRE(xsf::wilcoxon_sf(7.0, cdf_span) == 0.59375);
        REQUIRE(xsf::wilcoxon_sf(8.0, cdf_span) == 0.5);
        REQUIRE(xsf::wilcoxon_sf(9.0, cdf_span) == 0.40625);
    }
}

TEST_CASE("wilcoxon_sf edge cases", "[wilcoxon_sf]") {
    std::vector<double> cdf_table = wilcoxon_cdf_table(2);
    std::mdspan cdf_span(cdf_table.data(), cdf_table.size());

    REQUIRE(xsf::wilcoxon_sf(-1.0, cdf_span) == 1.0);
    REQUIRE(xsf::wilcoxon_sf(5.0, cdf_span) == 0.0);

    REQUIRE(xsf::wilcoxon_sf(2.9, cdf_span) == xsf::wilcoxon_sf(2.0, cdf_span));
    REQUIRE(xsf::wilcoxon_sf(-0.5, cdf_span) == xsf::wilcoxon_sf(0.0, cdf_span));

    REQUIRE(std::isnan(xsf::wilcoxon_sf(std::numeric_limits<double>::quiet_NaN(), cdf_span)));
}

TEST_CASE("wilcoxon_pmf_all exact small cases", "[wilcoxon_pmf]") {
    REQUIRE(wilcoxon_pmf_table(0) == std::vector<double>{1.0});
    REQUIRE(wilcoxon_pmf_table(1) == std::vector<double>{0.5, 0.5});
    REQUIRE(wilcoxon_pmf_table(2) == std::vector<double>{0.25, 0.25, 0.25, 0.25});
    REQUIRE(wilcoxon_pmf_table(3) == std::vector<double>{0.125, 0.125, 0.125, 0.25, 0.125, 0.125, 0.125});
    const std::vector<double> expected_n5 = {0.03125, 0.03125, 0.03125, 0.0625, 0.0625, 0.09375, 0.09375, 0.09375,
                                             0.09375, 0.09375, 0.09375, 0.0625, 0.0625, 0.03125, 0.03125, 0.03125};
    REQUIRE(wilcoxon_pmf_table(5) == expected_n5);

    std::vector<double> pmf_table = wilcoxon_pmf_table(3);
    std::mdspan pmf_span(pmf_table.data(), pmf_table.size());
    REQUIRE(xsf::wilcoxon_pmf(3.0, pmf_span) == 0.25);
    REQUIRE(xsf::wilcoxon_pmf(-1.0, pmf_span) == 0.0);
    REQUIRE(xsf::wilcoxon_pmf(7.0, pmf_span) == 0.0);
    REQUIRE(std::isnan(xsf::wilcoxon_pmf(std::numeric_limits<double>::quiet_NaN(), pmf_span)));
}

TEST_CASE("wilcoxon_cdf_all exact small cases", "[wilcoxon_cdf]") {
    REQUIRE(wilcoxon_cdf_table(3) == std::vector<double>{0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 1.0});
}
