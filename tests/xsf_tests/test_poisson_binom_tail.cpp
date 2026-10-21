#include <catch2/catch_test_macros.hpp>
#include <xsf/stats.h>
#include <xsf/third_party/kokkos/mdspan.hpp>
#include <cfenv>
#include <cmath>
#include <limits>
#include <vector>

TEST_CASE("poisson_binom direct tails: scipy issue 26409", "[poisson_binom][xsf_tests]") {
    std::vector<std::vector<double>> cases(3);
    for (int i = 0; i < 32; ++i) cases[0].push_back((1 + i % 9) / 128.0);
    cases[1] = std::vector<double>(256, 1.0 / 128);
    cases[2] = std::vector<double>(16, 1e-20);
    const int k[] = {24, 250, 1};
    // Exact binary64-input rational convolution, followed by 120-digit log.
    const double expected_sf[] = {1.196353171051496e-29, 0, 1.6e-19};
    const double expected_logsf[] = {-66.59568979136812, -1186.4216412494768, -43.27911313764113};
    for (int i = 0; i < 3; ++i) {
        auto &p = cases[i];
        std::vector<double> sf(p.size() + 1), logsf(p.size() + 1);
        xsf::poisson_binom_sf_all(std::mdspan(p.data(), p.size()), std::mdspan(sf.data(), sf.size()));
        xsf::poisson_binom_logsf_all(std::mdspan(p.data(), p.size()), std::mdspan(logsf.data(), logsf.size()));
        CAPTURE(i);
        REQUIRE(std::abs(sf[k[i] - 1] - expected_sf[i]) <= 2e-11 * expected_sf[i]);
        REQUIRE(std::abs(logsf[k[i] - 1] - expected_logsf[i]) <= 2e-11 * std::abs(expected_logsf[i]));
    }
}

template <typename T>
void check_poisson_binom_tail_contract() {
    const T inf = std::numeric_limits<T>::infinity();
    std::vector<T> p{0, 1, 0, 1}, sf(5, 99), logsf(5, 99);
    auto ps = std::mdspan(p.data(), p.size());
    auto ss = std::mdspan(sf.data(), sf.size());
    auto ls = std::mdspan(logsf.data(), logsf.size());
    std::feclearexcept(FE_ALL_EXCEPT);
    xsf::poisson_binom_sf_all(ps, ss);
    xsf::poisson_binom_logsf_all(ps, ls);
    REQUIRE(std::fetestexcept(FE_DIVBYZERO | FE_INVALID) == 0);
    for (int k = 0; k < 5; ++k) {
        REQUIRE(sf[k] == (k < 2 ? T(1) : T(0)));
        REQUIRE(logsf[k] == (k < 2 ? T(0) : -inf));
    }
    for (T bad : {T(-.1), T(1.1), inf, -inf, std::numeric_limits<T>::quiet_NaN()}) {
        p[2] = bad;
        xsf::poisson_binom_sf_all(ps, ss);
        xsf::poisson_binom_logsf_all(ps, ls);
        for (int k = 0; k < 5; ++k) {
            REQUIRE(std::isnan(sf[k]));
            REQUIRE(std::isnan(logsf[k]));
        }
    }
    p[2] = 0;
    std::vector<T> short_output(4, 17);
    xsf::poisson_binom_sf_all(ps, std::mdspan(short_output.data(), short_output.size()));
    xsf::poisson_binom_logsf_all(ps, std::mdspan(short_output.data(), short_output.size()));
    for (auto v : short_output) REQUIRE(v == T(17));
    p.clear();
    sf.resize(1); logsf.resize(1);
    xsf::poisson_binom_sf_all(std::mdspan(p.data(), p.size()), std::mdspan(sf.data(), sf.size()));
    xsf::poisson_binom_logsf_all(std::mdspan(p.data(), p.size()), std::mdspan(logsf.data(), logsf.size()));
    REQUIRE(sf[0] == 0);
    REQUIRE(logsf[0] == -inf);
}

TEST_CASE("poisson_binom tail endpoints and domain", "[poisson_binom][xsf_tests]") {
    check_poisson_binom_tail_contract<double>();
    check_poisson_binom_tail_contract<float>();
}

TEST_CASE("poisson_binom log tail near one and subnormal SF", "[poisson_binom][xsf_tests]") {
    std::vector<double> p(64, .5), sf(65), logsf(65);
    xsf::poisson_binom_sf_all(std::mdspan(p.data(), p.size()), std::mdspan(sf.data(), sf.size()));
    xsf::poisson_binom_logsf_all(std::mdspan(p.data(), p.size()), std::mdspan(logsf.data(), logsf.size()));
    REQUIRE(sf[0] == 1);
    const double expected = std::log1p(-std::ldexp(1., -64));
    REQUIRE(std::abs(logsf[0] - expected) <= 2e-11 * std::abs(expected));
    p = {std::numeric_limits<double>::denorm_min(), .5};
    sf.resize(3); logsf.resize(3);
    xsf::poisson_binom_sf_all(std::mdspan(p.data(), p.size()), std::mdspan(sf.data(), sf.size()));
    xsf::poisson_binom_logsf_all(std::mdspan(p.data(), p.size()), std::mdspan(logsf.data(), logsf.size()));
    REQUIRE(sf[1] == 0); // exact half-denorm rounds to even (zero)
    REQUIRE(std::isfinite(logsf[1]));
}
