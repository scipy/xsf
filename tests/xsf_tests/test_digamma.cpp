#include "../testing_utils.h"
#include <xsf/digamma.h>

#include <array>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <complex>
#include <limits>

using Catch::Matchers::IsNaN;
using Catch::Matchers::WithinAbs;

// Regression: digamma must short-circuit non-finite complex inputs
// before std::abs(z), which can return NaN for inf operands under
// libstdc++ and false-route every subsequent dispatch comparison.

TEMPLATE_TEST_CASE("digamma non-finite complex inputs",
                   "[digamma][xsf_tests][nonfinite]", float, double) {
    using T = TestType;
    using cd = std::complex<T>;
    const T inf = std::numeric_limits<T>::infinity();
    const T qnan = std::numeric_limits<T>::quiet_NaN();

    SECTION("+inf real -> +inf + 0j") {
        for (T im : std::array<T, 4>{T{0}, T{1}, T{-1}, T{1e30}}) {
            cd r = xsf::digamma(cd{inf, im});
            CAPTURE(im, r);
            REQUIRE(r.real() == inf);
            REQUIRE(r.imag() == T{0});
        }
    }

    SECTION("-inf real -> NaN") {
        for (T im : std::array<T, 3>{T{0}, T{1}, T{-1e30}}) {
            cd r = xsf::digamma(cd{-inf, im});
            CAPTURE(im, r);
            REQUIRE_THAT(r.real(), IsNaN());
            REQUIRE_THAT(r.imag(), IsNaN());
        }
    }

    SECTION("imag inf -> NaN") {
        for (cd z : std::array<cd, 7>{{
                {T{0}, inf}, {T{0}, -inf}, {T{1}, inf}, {T{-1}, -inf},
                {inf, inf}, {-inf, inf}, {inf, -inf}}}) {
            cd r = xsf::digamma(z);
            CAPTURE(z, r);
            REQUIRE_THAT(r.real(), IsNaN());
            REQUIRE_THAT(r.imag(), IsNaN());
        }
    }

    SECTION("any NaN component -> NaN") {
        for (cd z : std::array<cd, 8>{{
                {qnan, T{0}}, {T{0}, qnan}, {qnan, qnan},
                {qnan, T{1}}, {T{1}, qnan},
                {qnan, inf}, {inf, qnan}, {-inf, qnan}}}) {
            cd r = xsf::digamma(z);
            CAPTURE(z, r);
            REQUIRE_THAT(r.real(), IsNaN());
            REQUIRE_THAT(r.imag(), IsNaN());
        }
    }

    SECTION("finite inputs still work") {
        cd v1 = xsf::digamma(cd{T{1}, T{0}});
        REQUIRE_THAT(v1.real(), WithinAbs(-0.5772156649015329, 1e-6));
        cd v2 = xsf::digamma(cd{T{2}, T{0}});
        REQUIRE_THAT(v2.real(), WithinAbs(0.4227843350984671, 1e-6));
        cd v3 = xsf::digamma(cd{T{1}, T{1}});
        REQUIRE(std::isfinite(v3.real()));
        REQUIRE(std::isfinite(v3.imag()));
    }
}

TEMPLATE_TEST_CASE("digamma non-finite enumeration",
                   "[digamma][xsf_tests][nonfinite]", float, double) {
    using T = TestType;
    using cd = std::complex<T>;
    const T inf = std::numeric_limits<T>::infinity();
    const T qnan = std::numeric_limits<T>::quiet_NaN();

    const std::array<T, 4> values{T{1.5}, inf, -inf, qnan};
    for (T re : values) {
        for (T im : values) {
            cd r = xsf::digamma(cd{re, im});
            CAPTURE(re, im, r);
            if (re == inf && std::isfinite(im)) {
                REQUIRE(r.real() == inf);
                REQUIRE(r.imag() == T{0});
            } else if (!std::isfinite(re) || !std::isfinite(im)) {
                REQUIRE_THAT(r.real(), IsNaN());
                REQUIRE_THAT(r.imag(), IsNaN());
            } else {
                REQUIRE(std::isfinite(r.real()));
                REQUIRE(std::isfinite(r.imag()));
            }
        }
    }
}
