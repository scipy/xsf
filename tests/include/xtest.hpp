#pragma once
#include <type_traits>
#include <cmath>
#include <limits>


namespace xsf::test {

// Constant
constexpr float Inf32 = std::numeric_limits<float>::infinity();
constexpr float NaN32 = std::numeric_limits<float>::quiet_NaN();
constexpr double Inf64 = std::numeric_limits<double>::infinity();
constexpr double NaN64 = std::numeric_limits<double>::quiet_NaN();

/*

copy form: julia's isapprox:
    https://github.com/JuliaLang/julia/blob/af9e6e3167f7e444140c81326a2c3cf058ddba1a/base/floatfuncs.jl#L220

*/

// Floating-point types
template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
constexpr T rel_tol_default() {
    return std::sqrt(std::numeric_limits<T>::epsilon());
};
// Real (non-floating-point) types
template <typename T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
constexpr T rel_tol_default() { return 0; };

template <typename T>
bool isapprox(T x, T y, T rel_tol=rel_tol_default<T>()) {
    assert(rel_tol >= 0);
    return (x == y)
        // rel_error <= rel_tol
        || std::abs(x - y) <= (rel_tol * std::max(std::abs(x), std::abs(y)))
        // NaN == x == y
        || std::isnan(x) && std::isnan(y);
}

} // xsf::test