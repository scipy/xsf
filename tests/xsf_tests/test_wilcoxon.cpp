#include "../testing_utils.h"
#include <xsf/stats.h>

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
        const std::vector<double> pmf = xsf::wilcoxon_pmf(0);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_cdf(0.0, cdf_table, 0) == 1.0);
    }

    {   // n = 1
        const std::vector<double> pmf = xsf::wilcoxon_pmf(1);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_cdf(0.0, cdf_table, 1) == 0.5);
        REQUIRE(xsf::wilcoxon_cdf(1.0, cdf_table, 1) == 1.0);
    }

    {   // n = 2
        const std::vector<double> pmf = xsf::wilcoxon_pmf(2);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_cdf(1.0, cdf_table, 2) == 0.5);
        REQUIRE(xsf::wilcoxon_cdf(2.0, cdf_table, 2) == 0.75);
        REQUIRE(xsf::wilcoxon_cdf(3.0, cdf_table, 2) == 1.0);
        REQUIRE(xsf::wilcoxon_cdf(4.0, cdf_table, 2) == 1.0);
    }

    {   // n = 3
        const std::vector<double> pmf = xsf::wilcoxon_pmf(3);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_cdf(3.0, cdf_table, 3) == 0.625);
        REQUIRE(xsf::wilcoxon_cdf(4.0, cdf_table, 3) == 0.75);
        REQUIRE(xsf::wilcoxon_cdf(5.0, cdf_table, 3) == 0.875);
        REQUIRE(xsf::wilcoxon_cdf(6.0, cdf_table, 3) == 1.0);
        REQUIRE(xsf::wilcoxon_cdf(7.0, cdf_table, 3) == 1.0);
    }

    {   // n = 5
        const std::vector<double> pmf = xsf::wilcoxon_pmf(5);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_cdf(5.0, cdf_table, 5) == 0.3125);
        REQUIRE(xsf::wilcoxon_cdf(6.0, cdf_table, 5) == 0.40625);
        REQUIRE(xsf::wilcoxon_cdf(7.0, cdf_table, 5) == 0.5);
        REQUIRE(xsf::wilcoxon_cdf(8.0, cdf_table, 5) == 0.59375);
        REQUIRE(xsf::wilcoxon_cdf(9.0, cdf_table, 5) == 0.6875);
    }
}

TEST_CASE("wilcoxon_cdf edge cases", "[wilcoxon_cdf]") {
    const std::vector<double> pmf = xsf::wilcoxon_pmf(2);
    const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);

    REQUIRE(xsf::wilcoxon_cdf(4.0, cdf_table, 2) == 1.0);

    REQUIRE(xsf::wilcoxon_cdf(1.9, cdf_table, 2) == xsf::wilcoxon_cdf(1.0, cdf_table, 2));
    REQUIRE(xsf::wilcoxon_cdf(-0.5, cdf_table, 2) == xsf::wilcoxon_cdf(0.0, cdf_table, 2));

    REQUIRE(std::isnan(xsf::wilcoxon_cdf(0.0, cdf_table, -1)));
    REQUIRE(std::isnan(xsf::wilcoxon_cdf(std::numeric_limits<double>::quiet_NaN(), cdf_table, 2)));
}

TEST_CASE("wilcoxon_sf exact small cases", "[wilcoxon_sf]") {
    {
        // n = 0
        const std::vector<double> pmf = xsf::wilcoxon_pmf(0);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_sf(0.0, cdf_table, 0) == 1.0);
        REQUIRE(xsf::wilcoxon_sf(1.0, cdf_table, 0) == 0.0);
    }

    {   // n = 1
        const std::vector<double> pmf = xsf::wilcoxon_pmf(1);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_sf(0.0, cdf_table, 1) == 1.0);
        REQUIRE(xsf::wilcoxon_sf(1.0, cdf_table, 1) == 0.5);
        REQUIRE(xsf::wilcoxon_sf(2.0, cdf_table, 1) == 0.0);
    }

    {   // n = 2
        const std::vector<double> pmf = xsf::wilcoxon_pmf(2);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_sf(0.0, cdf_table, 2) == 1.0);
        REQUIRE(xsf::wilcoxon_sf(1.0, cdf_table, 2) == 0.75);
        REQUIRE(xsf::wilcoxon_sf(2.0, cdf_table, 2) == 0.5);
        REQUIRE(xsf::wilcoxon_sf(3.0, cdf_table, 2) == 0.25);
        REQUIRE(xsf::wilcoxon_sf(4.0, cdf_table, 2) == 0.0);
    }

    {   // n = 3
        const std::vector<double> pmf = xsf::wilcoxon_pmf(3);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_sf(3.0, cdf_table, 3) == 0.625);
        REQUIRE(xsf::wilcoxon_sf(4.0, cdf_table, 3) == 0.375);
        REQUIRE(xsf::wilcoxon_sf(5.0, cdf_table, 3) == 0.25);
        REQUIRE(xsf::wilcoxon_sf(6.0, cdf_table, 3) == 0.125);
        REQUIRE(xsf::wilcoxon_sf(7.0, cdf_table, 3) == 0.0);
    }

    {   // n = 5
        const std::vector<double> pmf = xsf::wilcoxon_pmf(5);
        const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);
        REQUIRE(xsf::wilcoxon_sf(5.0, cdf_table, 5) == 0.78125);
        REQUIRE(xsf::wilcoxon_sf(6.0, cdf_table, 5) == 0.6875);
        REQUIRE(xsf::wilcoxon_sf(7.0, cdf_table, 5) == 0.59375);
        REQUIRE(xsf::wilcoxon_sf(8.0, cdf_table, 5) == 0.5);
        REQUIRE(xsf::wilcoxon_sf(9.0, cdf_table, 5) == 0.40625);
    }
}

TEST_CASE("wilcoxon_sf edge cases", "[wilcoxon_sf]") {
    const std::vector<double> pmf = xsf::wilcoxon_pmf(2);
    const std::vector<double> cdf_table = xsf::wilcoxon_cdf_table(pmf);

    REQUIRE(xsf::wilcoxon_sf(-1.0, cdf_table, 2) == 1.0);
    REQUIRE(xsf::wilcoxon_sf(5.0, cdf_table, 2) == 0.0);

    REQUIRE(xsf::wilcoxon_sf(2.9, cdf_table, 2) == xsf::wilcoxon_sf(2.0, cdf_table, 2));
    REQUIRE(xsf::wilcoxon_sf(-0.5, cdf_table, 2) == xsf::wilcoxon_sf(0.0, cdf_table, 2));

    REQUIRE(std::isnan(xsf::wilcoxon_sf(0.0, cdf_table, -1)));
    REQUIRE(std::isnan(xsf::wilcoxon_sf(std::numeric_limits<double>::quiet_NaN(), cdf_table, 2)));
}
