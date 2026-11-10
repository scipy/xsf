/*

Implementation of spherical Bessel functions and modified spherical Bessel
functions of the first and second kinds, as well as their derivatives.

Author: Tadeusz Pudlik

Distributed under the same license as SciPy.

I attempt to correctly handle the edge cases (0 and infinity), but this is
tricky: the values of the functions often depend on the direction in which
the limit is taken. At zero, I follow the convention of numpy (1.9.2),
which treats zero differently depending on its type:

    >>> np.cos(0)/0
    inf
    >>> np.cos(0+0j)/(0+0j)
    inf + nan*j

So, real zero is assumed to be "positive zero", while complex zero has an
unspecified sign and produces nans.  Similarly, complex infinity is taken to
represent the "point at infinity", an ambiguity which for some functions
makes `nan` the correct return value.

Translated to C++ by SciPy developers in 2024.

*/

#pragma once

#include "../error.h"
#include "amos.h"
#include "bessel.h"

namespace xsf {

namespace detail {

    // Leading term x^n / (2n + 1)!! of j_n(x) and i_n(x) as x -> 0, https://dlmf.nist.gov/10.52.E1
    template <typename T>
    T sph_bessel_tiny(long n, T x) {
        T res = 1;
        for (long k = 1; k <= n && res != static_cast<T>(0); ++k) {
            res *= x / (2 * static_cast<remove_complex_t<T>>(k) + 1);
        }
        return res;
    }

} // namespace detail

template <typename T>
T sph_bessel_j(long n, T x) {
    if (std::isnan(x)) {
        return x;
    }

    if (n < 0) {
        set_error("spherical_jn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    // For real negative x, the principal-branch product
    // sqrt(M_PI_2/x) * J_{n+1/2}(x) introduces a spurious sign relative
    // to the (entire) function j_n. Use the reflection identity
    // j_n(-x) = (-1)^n * j_n(x) explicitly (DLMF 10.47.14).
    if (x < 0) {
        T r = sph_bessel_j(n, -x);
        return (n % 2 == 0) ? r : -r;
    }

    if ((x == std::numeric_limits<T>::infinity()) || (x == -std::numeric_limits<T>::infinity())) {
        return 0;
    }

    if (x == 0) {
        if (n == 0) {
            return 1;
        }

        return 0;
    }

    // Avoid underflow of J_{n + 1/2}(x) as in sph_bessel_i; the first correction term,
    // -x^2 / (4n + 6) (https://dlmf.nist.gov/10.53.E1), is below eps / 4 here.
    // x < 1 keeps the product from overflowing for huge n.
    if (x < 1 && x * x < (static_cast<T>(n) + static_cast<T>(1.5)) * std::numeric_limits<T>::epsilon()) {
        return detail::sph_bessel_tiny(n, x);
    }

    if ((n > 0) && (n >= x)) {
        return std::sqrt(M_PI_2 / x) * cyl_bessel_j(n + 1 / static_cast<T>(2), x);
    }

    T s0 = std::sin(x) / x;
    if (n == 0) {
        return s0;
    }

    T s1 = (s0 - std::cos(x)) / x;
    if (n == 1) {
        return s1;
    }

    T sn;
    for (int i = 0; i < n - 1; ++i) {
        sn = (2 * i + 3) * s1 / x - s0;
        s0 = s1;
        s1 = sn;
        if (std::isinf(sn)) {
            // Overflow occurred already : terminate recurrence.
            return sn;
        }
    }

    return sn;
}

template <typename T>
std::complex<T> sph_bessel_j(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_jn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::real(z) == std::numeric_limits<T>::infinity() || std::real(z) == -std::numeric_limits<T>::infinity()) {
        // https://dlmf.nist.gov/10.52.E3
        if (std::imag(z) == 0) {
            return 0;
        }

        return std::complex<T>(1, 1) * std::numeric_limits<T>::infinity();
    }

    if ((std::real(z) == 0) && (std::imag(z) == 0)) {
        if (n == 0) {
            return 1;
        }

        return 0;
    }

    // For purely real negative z, the principal-branch product
    // sqrt(M_PI_2/z) * J_{n+1/2}(z) introduces a spurious sign relative
    // to the (entire) function j_n. Use the reflection identity
    // j_n(-z) = (-1)^n * j_n(z) explicitly (DLMF 10.47.14).
    if (std::imag(z) == 0 && std::real(z) < 0) {
        std::complex<T> r = sph_bessel_j(n, std::complex<T>{-std::real(z), static_cast<T>(0)});
        return (n % 2 == 0) ? r : -r;
    }

    std::complex<T> out = std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_j(n + 1 / static_cast<T>(2), z);
    if (std::imag(z) == 0) {
        return std::real(out); // Small imaginary part is spurious
    }

    return out;
}

template <typename T>
T sph_bessel_j_jac(long n, T z) {
    if (n == 0) {
        return -sph_bessel_j(1, z);
    }

    // Leading term n z^(n - 1) / (2n + 1)!! of the derivative of https://dlmf.nist.gov/10.52.E1,
    // as DLMF 10.51.2 fails when j_n(z) underflows; the first correction term is below eps / 2 here.
    // |z| < 1 keeps the product from overflowing for huge n, as in sph_bessel_j.
    using R = remove_complex_t<T>;
    if (std::abs(z) < 1 &&
        std::norm(z) < (static_cast<R>(n) + static_cast<R>(0.5)) * std::numeric_limits<R>::epsilon()) {
        return static_cast<R>(n) / (2 * static_cast<R>(n) + 1) * detail::sph_bessel_tiny(n - 1, z);
    }

    // DLMF 10.51.2
    return sph_bessel_j(n - 1, z) - static_cast<T>(n + 1) * sph_bessel_j(n, z) / z;
}

template <typename T>
T sph_bessel_y(long n, T x) {
    T s0, s1, sn;
    int idx;

    if (std::isnan(x)) {
        return x;
    }

    if (n < 0) {
        set_error("spherical_yn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (x < 0) {
        T r = sph_bessel_y(n, -x);
        return (n % 2 == 1) ? r : -r;
    }

    if (x == std::numeric_limits<T>::infinity() || x == -std::numeric_limits<T>::infinity()) {
        return 0;
    }

    if (x == 0) {
        return -std::numeric_limits<T>::infinity();
    }

    s0 = -cos(x) / x;
    if (n == 0) {
        return s0;
    }

    s1 = (s0 - sin(x)) / x;
    // If s1 overflows (tiny x), so does y_n(x); if s0 does too, the recurrence would give inf - inf.
    if (n == 1 || std::isinf(s1)) {
        return s1;
    }

    for (idx = 0; idx < n - 1; ++idx) {
        sn = (2 * idx + 3) * s1 / x - s0;
        s0 = s1;
        s1 = sn;
        if (std::isinf(sn)) {
            // Overflow occurred already: terminate recurrence.
            return sn;
        }
    }

    return sn;
}

inline float sph_bessel_y(long n, float x) { return sph_bessel_y(n, static_cast<double>(x)); }

template <typename T>
std::complex<T> sph_bessel_y(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_yn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::real(z) == 0 && std::imag(z) == 0) {
        // https://dlmf.nist.gov/10.52.E2
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::real(z) == std::numeric_limits<T>::infinity() || std::real(z) == -std::numeric_limits<T>::infinity()) {
        // https://dlmf.nist.gov/10.52.E3
        if (std::imag(z) == 0) {
            return 0;
        }

        return std::complex<T>(1, 1) * std::numeric_limits<T>::infinity();
    }

    // y_n(-z) = (-1)^(n+1) * y_n(z); see sph_bessel_j note above.
    if (std::imag(z) == 0 && std::real(z) < 0) {
        std::complex<T> r = sph_bessel_y(n, std::complex<T>{-std::real(z), static_cast<T>(0)});
        return (n % 2 == 0) ? -r : r;
    }

    return std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_y(n + 1 / static_cast<T>(2), z);
}

template <typename T>
T sph_bessel_y_jac(long n, T x) {
    if (n == 0) {
        return -sph_bessel_y(1, x);
    }

    // DLMF 10.51.2; for real x, y_{n-1}(x) is negligible if the second term overflows, so avoid inf - inf.
    T res = -static_cast<T>(n + 1) * sph_bessel_y(n, x) / x;
    if constexpr (cxx::is_floating_point_v<T>) {
        if (std::isinf(res)) {
            return res;
        }
    }

    return sph_bessel_y(n - 1, x) + res;
}

template <typename T>
T sph_bessel_i(long n, T x) {
    if (std::isnan(x)) {
        return x;
    }

    if (n < 0) {
        set_error("spherical_in", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (x < 0) {
        T r = sph_bessel_i(n, -x);
        return (n % 2 == 0) ? r : -r;
    }

    if (x == 0) {
        // https://dlmf.nist.gov/10.52.E1
        if (n == 0) {
            return 1;
        }
        return 0;
    }

    if (std::isinf(x)) {
        // https://dlmf.nist.gov/10.49.E8
        if (x == -std::numeric_limits<T>::infinity()) {
            return std::pow(-1, n) * std::numeric_limits<T>::infinity();
        }

        return std::numeric_limits<T>::infinity();
    }

    // Computing i_n through I_{n + 1/2} can underflow before multiplication
    // by the spherical normalization.  In the tiny-x region, use the leading
    // term of the power series i_n(x) = x^n / (2n + 1)!! * (1 + O(x^2)).
    // The cutoff bounds the omitted first correction by roughly eps / 4.
    if (x * x < (static_cast<T>(n) + static_cast<T>(1.5)) * std::numeric_limits<T>::epsilon()) {
        return detail::sph_bessel_tiny(n, x);
    }

    return sqrt(static_cast<T>(M_PI_2) / x) * cyl_bessel_i(n + 1 / static_cast<T>(2), x);
}

template <typename T>
std::complex<T> sph_bessel_i(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_in", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::abs(z) == 0) {
        // https://dlmf.nist.gov/10.52.E1
        if (n == 0) {
            return 1;
        }

        return 0;
    }

    if (std::isinf(std::real(z)) || std::isinf(std::imag(z))) {
        // https://dlmf.nist.gov/10.52.E5
        if (std::imag(z) == 0) {
            if (std::real(z) == -std::numeric_limits<T>::infinity()) {
                return std::pow(-1, n) * std::numeric_limits<T>::infinity();
            }

            return std::numeric_limits<T>::infinity();
        }

        return std::numeric_limits<T>::quiet_NaN();
    }

    // i_n(-z) = (-1)^n * i_n(z); see sph_bessel_j note above.
    if (std::imag(z) == 0 && std::real(z) < 0) {
        std::complex<T> r = sph_bessel_i(n, std::complex<T>{-std::real(z), static_cast<T>(0)});
        return (n % 2 == 0) ? r : -r;
    }

    return std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_i(n + 1 / static_cast<T>(2), z);
}

template <typename T>
T sph_bessel_i_jac(long n, T z) {
    if (n == 0) {
        return sph_bessel_i(1, z);
    }

    // Leading term n z^(n - 1) / (2n + 1)!! of the derivative of https://dlmf.nist.gov/10.52.E1,
    // as DLMF 10.51.5 fails when i_n(z) underflows; the first correction term is below eps / 2 here.
    // |z| < 1 keeps the product from overflowing for huge n, as in sph_bessel_j.
    using R = remove_complex_t<T>;
    if (std::abs(z) < 1 &&
        std::norm(z) < (static_cast<R>(n) + static_cast<R>(0.5)) * std::numeric_limits<R>::epsilon()) {
        return static_cast<R>(n) / (2 * static_cast<R>(n) + 1) * detail::sph_bessel_tiny(n - 1, z);
    }

    return sph_bessel_i(n - 1, z) - static_cast<T>(n + 1) * sph_bessel_i(n, z) / z;
}

template <typename T>
T sph_bessel_k(long n, T z) {
    if (std::isnan(z)) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_kn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (z == 0) {
        return std::numeric_limits<T>::infinity();
    }

    if (z < 0) {
        // https://dlmf.nist.gov/10.47#E17
        // and
        // https://dlmf.nist.gov/10.47.E11
        int sign = (std::abs(n) % 2 == 0) ? 1 : -1;
        return -(M_PI * sph_bessel_i(n, -z) + sign * sph_bessel_k(n, -z));
    }

    if (std::isinf(z)) {
        // https://dlmf.nist.gov/10.52.E6
        if (z == std::numeric_limits<T>::infinity()) {
            return 0;
        }

        return -std::numeric_limits<T>::infinity();
    }

    return std::sqrt(M_PI_2 / z) * cyl_bessel_k(n + 1 / static_cast<T>(2), z);
}

template <typename T>
std::complex<T> sph_bessel_k(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_kn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::abs(z) == 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::isinf(std::real(z)) || std::isinf(std::imag(z))) {
        // https://dlmf.nist.gov/10.52.E6
        if (std::imag(z) == 0) {
            if (std::real(z) == std::numeric_limits<T>::infinity()) {
                return 0;
            }

            return -std::numeric_limits<T>::infinity();
        }

        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::imag(z) == 0) {
        return {sph_bessel_k(n, std::real(z)), 0};
    }

    return std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_k(n + 1 / static_cast<T>(2), z);
}

template <typename T>
T sph_bessel_k_jac(long n, T x) {
    if (n == 0) {
        return -sph_bessel_k(1, x);
    }

    return -sph_bessel_k(n - 1, x) - static_cast<T>(n + 1) * sph_bessel_k(n, x) / x;
}

} // namespace xsf
