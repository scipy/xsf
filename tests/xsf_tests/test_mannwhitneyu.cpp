#include "../testing_utils.h"
#include <xsf/stats.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

namespace {

std::vector<double> mannwhitneyu_pmf_table(int n1, int n2) {
    std::vector<double> table(n1 * n2 + 1, 0.0);
    std::mdspan table_span(table.data(), table.size());
    xsf::mannwhitneyu_pmf_all(n1, n2, table_span);
    return table;
}

std::vector<double> mannwhitneyu_cdf_table(int n1, int n2) {
    std::vector<double> table(n1 * n2 + 1, 0.0);
    std::mdspan table_span(table.data(), table.size());
    xsf::mannwhitneyu_cdf_all(n1, n2, table_span);
    return table;
}

std::vector<double> divide(const std::vector<double> &counts, double total) {
    std::vector<double> out;
    for (double c : counts) {
        out.push_back(c / total);
    }
    return out;
}

} // namespace

// Reference counts from brute-force enumeration of all binom(n1 + n2, n1) rank
// assignments, e.g. in Python:
//
// from itertools import combinations
// n1, n2 = 3, 4
// counts = [0] * (n1 * n2 + 1)
// for S in combinations(range(n1 + n2), n1):
//     counts[sum(S) - n1 * (n1 - 1) // 2] += 1

TEST_CASE("mannwhitneyu_pmf_all exact small cases", "[mannwhitneyu_pmf]") {
    REQUIRE(mannwhitneyu_pmf_table(0, 0) == std::vector<double>{1.0});
    REQUIRE(mannwhitneyu_pmf_table(0, 3) == std::vector<double>{1.0});
    REQUIRE(mannwhitneyu_pmf_table(1, 1) == std::vector<double>{0.5, 0.5});
    REQUIRE(mannwhitneyu_pmf_table(1, 3) == std::vector<double>{0.25, 0.25, 0.25, 0.25});
    REQUIRE(mannwhitneyu_pmf_table(2, 2) == divide({1, 1, 2, 1, 1}, 6.0));
    REQUIRE(mannwhitneyu_pmf_table(3, 4) == divide({1, 1, 2, 3, 4, 4, 5, 4, 4, 3, 2, 1, 1}, 35.0));
    REQUIRE(mannwhitneyu_pmf_table(4, 4) == divide({1, 1, 2, 3, 5, 5, 7, 7, 8, 7, 7, 5, 5, 3, 2, 1, 1}, 70.0));
    REQUIRE(
        mannwhitneyu_pmf_table(5, 7) == divide(
                                            {1,  1,  2,  3,  5,  7,  10, 13, 17, 21, 26, 30, 35, 39, 43, 46, 48, 49,
                                             49, 48, 46, 43, 39, 35, 30, 26, 21, 17, 13, 10, 7,  5,  3,  2,  1,  1},
                                            792.0
                                        )
    );

    // The distribution does not depend on the order of n1 and n2.
    REQUIRE(mannwhitneyu_pmf_table(4, 3) == mannwhitneyu_pmf_table(3, 4));
    REQUIRE(mannwhitneyu_pmf_table(7, 5) == mannwhitneyu_pmf_table(5, 7));

    std::vector<double> pmf_table = mannwhitneyu_pmf_table(2, 2);
    std::mdspan pmf_span(pmf_table.data(), pmf_table.size());
    REQUIRE(xsf::mannwhitneyu_pmf(2.0, pmf_span) == 2.0 / 6.0);
    REQUIRE(xsf::mannwhitneyu_pmf(-1.0, pmf_span) == 0.0);
    REQUIRE(xsf::mannwhitneyu_pmf(5.0, pmf_span) == 0.0);
    REQUIRE(std::isnan(xsf::mannwhitneyu_pmf(std::numeric_limits<double>::quiet_NaN(), pmf_span)));
}

TEST_CASE("mannwhitneyu_cdf exact small cases", "[mannwhitneyu_cdf]") {
    { // n1 = n2 = 2
        std::vector<double> cdf_table = mannwhitneyu_cdf_table(2, 2);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::mannwhitneyu_cdf(0.0, cdf_span) == 1.0 / 6.0);
        REQUIRE(xsf::mannwhitneyu_cdf(1.0, cdf_span) == 2.0 / 6.0);
        REQUIRE(xsf::mannwhitneyu_cdf(2.0, cdf_span) == 4.0 / 6.0);
        REQUIRE(xsf::mannwhitneyu_cdf(4.0, cdf_span) == 1.0);
        REQUIRE(xsf::mannwhitneyu_cdf(5.0, cdf_span) == 1.0);
        REQUIRE(xsf::mannwhitneyu_cdf(-1.0, cdf_span) == 0.0);
        REQUIRE(std::isnan(xsf::mannwhitneyu_cdf(std::numeric_limits<double>::quiet_NaN(), cdf_span)));
    }

    { // n1 = n2 = 4: P(U <= 0) = 1 / binom(8, 4) = 1 / 70
        std::vector<double> cdf_table = mannwhitneyu_cdf_table(4, 4);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::mannwhitneyu_cdf(0.0, cdf_span) == 1.0 / 70.0);
        // P(U <= 3) = (1 + 1 + 2 + 3) / 70
        REQUIRE(xsf::extended_relative_error(xsf::mannwhitneyu_cdf(3.0, cdf_span), 7.0 / 70.0) <= 1e-15);
    }

    { // n1 = 5, n2 = 7: P(U <= 8) = (1 + 1 + 2 + 3 + 5 + 7 + 10 + 13 + 17) / 792
        std::vector<double> cdf_table = mannwhitneyu_cdf_table(5, 7);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::extended_relative_error(xsf::mannwhitneyu_cdf(8.0, cdf_span), 59.0 / 792.0) <= 1e-15);
    }
}

TEST_CASE("mannwhitneyu_sf exact small cases", "[mannwhitneyu_sf]") {
    { // n1 = n2 = 2; sf includes the mass at k: P(U >= k)
        std::vector<double> cdf_table = mannwhitneyu_cdf_table(2, 2);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::mannwhitneyu_sf(0.0, cdf_span) == 1.0);
        REQUIRE(xsf::extended_relative_error(xsf::mannwhitneyu_sf(1.0, cdf_span), 5.0 / 6.0) <= 1e-15);
        REQUIRE(xsf::mannwhitneyu_sf(2.0, cdf_span) == 4.0 / 6.0);
        REQUIRE(xsf::mannwhitneyu_sf(3.0, cdf_span) == 2.0 / 6.0);
        REQUIRE(xsf::mannwhitneyu_sf(4.0, cdf_span) == 1.0 / 6.0);
        REQUIRE(xsf::mannwhitneyu_sf(5.0, cdf_span) == 0.0);
        REQUIRE(xsf::mannwhitneyu_sf(-1.0, cdf_span) == 1.0);
        REQUIRE(std::isnan(xsf::mannwhitneyu_sf(std::numeric_limits<double>::quiet_NaN(), cdf_span)));
    }

    { // n1 = n2 = 4: P(U >= 16) = 1 / 70
        std::vector<double> cdf_table = mannwhitneyu_cdf_table(4, 4);
        std::mdspan cdf_span(cdf_table.data(), cdf_table.size());
        REQUIRE(xsf::mannwhitneyu_sf(16.0, cdf_span) == 1.0 / 70.0);
        REQUIRE(xsf::mannwhitneyu_sf(16.0, cdf_span) == xsf::mannwhitneyu_cdf(0.0, cdf_span));
    }
}

TEST_CASE("mannwhitneyu large n1, n2", "[mannwhitneyu_sf]") {
    // Here binom(n1 + n2, n1) is far above 2^64, so the counts leave the exact
    // integer range and the floating-point path is used.
    // Reference values from an exact computation with Python integers and
    // fractions.Fraction (same recurrence), rounded to double.
    std::vector<double> cdf_table = mannwhitneyu_cdf_table(50, 50);
    std::mdspan cdf_span(cdf_table.data(), cdf_table.size());

    auto [k, expected] = GENERATE(
        table<long long int, double>(
            {{0, 9.911653021418339e-30},
             {1, 1.9823306042836678e-29},
             {100, 1.6154986072803412e-20},
             {500, 3.5967866378289605e-08},
             {900, 0.007761972615356092},
             {1250, 0.5013688955007638}}
        )
    );
    CAPTURE(k, expected);
    // P(U <= k) = P(U >= n1 * n2 - k) by symmetry
    REQUIRE(xsf::extended_relative_error(xsf::mannwhitneyu_cdf(k, cdf_span), expected) <= 1e-14);
    REQUIRE(xsf::extended_relative_error(xsf::mannwhitneyu_sf(2500 - k, cdf_span), expected) <= 1e-14);
}

TEST_CASE("mannwhitneyu large unequal n1, n2", "[mannwhitneyu_sf]") {
    // Very unequal sample sizes. binom(110, 10) ~ 4.7e13 < 2^53, so the counts stay
    // in the exact integer path.
    // Reference values from an independent exact computation with Python integers
    // (counting partitions that fit in a 10 x 100 box), rounded to double.
    std::vector<double> cdf_table = mannwhitneyu_cdf_table(10, 100);
    std::mdspan cdf_span(cdf_table.data(), cdf_table.size());

    auto [k, expected] = GENERATE(
        table<long long int, double>(
            {{0, 2.1323036126913318e-14},
             {1, 4.2646072253826637e-14},
             {10, 2.963902021640951e-12},
             {100, 1.7477197978477222e-06},
             {250, 0.004017232966163219},
             {400, 0.15287598344036554},
             {500, 0.5020423598360503}}
        )
    );
    CAPTURE(k, expected);
    // P(U <= k) = P(U >= n1 * n2 - k) by symmetry
    REQUIRE(xsf::extended_relative_error(xsf::mannwhitneyu_cdf(k, cdf_span), expected) <= 1e-14);
    REQUIRE(xsf::extended_relative_error(xsf::mannwhitneyu_sf(1000 - k, cdf_span), expected) <= 1e-14);
    // Swapping the sample sizes gives the same distribution.
    std::vector<double> swapped = mannwhitneyu_cdf_table(100, 10);
    REQUIRE(swapped == cdf_table);
}

TEST_CASE("mannwhitneyu invalid arguments", "[mannwhitneyu_pmf]") {
    // A wrongly sized output is left untouched.
    std::vector<double> table(4, -1.0);
    std::mdspan table_span(table.data(), table.size());
    xsf::mannwhitneyu_pmf_all(2, 2, table_span);
    REQUIRE(table == std::vector<double>{-1.0, -1.0, -1.0, -1.0});
}
