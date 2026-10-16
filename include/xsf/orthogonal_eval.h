#pragma once

#include "beta.h"
#include "binom.h"
#include "config.h"
#include "error.h"
#include "hyp2f1.h"

namespace xsf {

namespace detail {

    template <typename T>
    XSF_HOST_DEVICE inline T eval_jacobi(double n, double alpha, double beta, T x) {
        double a, b, c, d;
        T g;

        if (alpha == -1 && cxx::abs(beta) == 1) {
            if (n == 0) {
                return 1.0;
            } else if (n == 1) {
                return 0.5 * (1.0 + beta) * (x - 1.0);
            } else if (n > 1) {
                return ((n + beta) / (2.0 * n)) * (x - 1.0) * eval_jacobi(n - 1.0, 1.0, beta, x);
            }
        }

        d = binom(n + alpha, n);
        a = -n;
        b = n + alpha + beta + 1.0;
        c = alpha + 1.0;
        g = 0.5 * (1.0 - x);
        return d * hyp2f1(a, b, c, g);
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_jacobi_l(Int n, double alpha, double beta, double x) {
        Int kk;
        double p, d;
        double k, t;

        if (n < 0) {
            return eval_jacobi(n, alpha, beta, x);
        } else if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return 0.5 * (2.0 * (alpha + 1.0) + (alpha + beta + 2.0) * (x - 1.0));
        } else if (alpha == -1 && cxx::abs(beta) == 1) {
            return ((n + beta) / (2.0 * n)) * (x - 1.0) * eval_jacobi(n - 1.0, 1.0, beta, x);
        } else {
            d = (alpha + beta + 2.0) * (x - 1.0) / (2.0 * (alpha + 1.0));
            p = d + 1.0;
            for (kk = 0; kk < n - 1; kk++) {
                k = kk + 1.0;
                t = 2.0 * k + alpha + beta;
                d = ((t * (t + 1.0) * (t + 2.0)) * (x - 1.0) * p + 2.0 * k * (k + beta) * (t + 2.0) * d) /
                    (2.0 * (k + alpha + 1.0) * (k + alpha + beta + 1.0) * t);
                p = d + p;
            }
            return binom(n + alpha, n) * p;
        }
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_hermitenorm(Int n, double x) {
        if (cxx::isnan(x)) {
            return x;
        }

        if (n < 0) {
            set_error("eval_hermitenorm", SF_ERROR_DOMAIN, "polynomial only defined for nonnegative n");
            return cxx::numeric_limits<double>::quiet_NaN();
        } else if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return x;
        }

        double y3 = 0.0;
        double y2 = 1.0;
        for (Int k = n; k > 1; --k) {
            const double y1 = x * y2 - k * y3;
            y3 = y2;
            y2 = y1;
        }
        return x * y2 - y3;
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_hermite(Int n, double x) {
        if (n < 0) {
            set_error("eval_hermite", SF_ERROR_DOMAIN, "polynomial only defined for nonnegative n");
            return cxx::numeric_limits<double>::quiet_NaN();
        }
        return eval_hermitenorm(n, cxx::sqrt(2.0) * x) * cxx::pow(2.0, n / 2.0);
    }

    // Integer-degree L_n^(alpha)(x) by forward recurrence on the normalized values
    // p_k = L_k^(alpha)(x) / binom(k + alpha, k) and their differences d_k = p_k - p_{k-1}:
    //   d_{k+1} = (k d_k - x p_k) / (k + alpha + 1),  p_{k+1} = p_k + d_{k+1},
    // starting from p_0 = 1, d_0 = 0.
    template <typename Int>
    XSF_HOST_DEVICE inline double eval_genlaguerre_l(Int n, double alpha, double x) {
        if (alpha <= -1) {
            set_error("eval_genlaguerre", SF_ERROR_DOMAIN, "polynomial defined only for alpha > -1");
            return cxx::numeric_limits<double>::quiet_NaN();
        }

        if (cxx::isnan(alpha) || cxx::isnan(x)) {
            return cxx::numeric_limits<double>::quiet_NaN();
        }

        if (n < 0) {
            return 0.0;
        } else if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return -x + alpha + 1.0;
        }

        double d = -x / (alpha + 1.0);
        double p = d + 1.0;
        for (Int kk = 0; kk < n - 1; ++kk) {
            const double k = kk + 1.0;
            d = -x / (k + alpha + 1.0) * p + (k / (k + alpha + 1.0)) * d;
            p = d + p;
        }
        return binom(n + alpha, n) * p;
    }

    // Integer-degree C_n^(alpha)(x). For |x| >= 1e-5, forward recurrence on the normalized values
    // p_k = C_k^(alpha)(x) / binom(k + 2 alpha - 1, k) and their differences d_k = p_k - p_{k-1}:
    //   d_{k+1} = 2 (k + alpha) / (k + 2 alpha) (x - 1) p_k + k / (k + 2 alpha) d_k,  p_{k+1} = p_k + d_{k+1},
    // starting from p_1 = x, d_1 = x - 1. For |x| < 1e-5, the explicit power series in x is used instead.
    template <typename Int>
    XSF_HOST_DEVICE inline double eval_gegenbauer_l(Int n, double alpha, double x) {
        if (cxx::isnan(alpha) || cxx::isnan(x)) {
            return cxx::numeric_limits<double>::quiet_NaN();
        }

        if (alpha == 0.0) {
            return 0.0;
        }
        if (n < 0) {
            return 0.0;
        } else if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return 2.0 * alpha * x;
        }

        // Symmetry property: C_n^alpha(-x) = (-1)^n C_n^alpha(x) for integer n
        const int sign = (x < 0 && n % 2) ? -1 : 1;
        x = cxx::abs(x);

        if (x < 1e-5) {
            // Power series rather than recurrence due to loss of precision
            // http://functions.wolfram.com/Polynomials/GegenbauerC3/02/
            const Int a = n / 2;

            double d = (a % 2 == 0) ? 1.0 : -1.0;
            d /= cephes::beta(alpha, 1.0 + a);
            if (n == 2 * a) {
                d /= (a + alpha);
            } else {
                d *= 2.0 * x;
            }

            double p = 0.0;
            for (Int kk = 0; kk <= a; ++kk) {
                p += d;
                d *= -4.0 * x * x * (a - kk) * (alpha - a + kk + n) /
                     ((n + 1.0 - 2.0 * a + 2.0 * kk) * (n + 2.0 - 2.0 * a + 2.0 * kk));
                if (cxx::abs(d) < 1e-20 * cxx::abs(p)) {
                    // converged
                    break;
                }
            }
            return sign * p;
        }

        double d = x - 1.0;
        double p = x;
        for (Int kk = 0; kk < n - 1; ++kk) {
            const double k = kk + 1.0;
            d = (2.0 * (k + alpha) / (k + 2.0 * alpha)) * (x - 1.0) * p + (k / (k + 2.0 * alpha)) * d;
            p = d + p;
        }

        if (cxx::abs(alpha / n) < 1e-8) {
            // avoid loss of precision
            return sign * 2.0 * alpha / n * p;
        }
        return sign * binom(n + 2.0 * alpha - 1.0, n) * p;
    }

    template <typename T>
    XSF_HOST_DEVICE inline T eval_gegenbauer(double n, double alpha, T x) {
        // If n is an integer, use more stable `eval_gegenbauer_l`
        if constexpr (cxx::is_same_v<T, double>) {
            if (n >= 0 && n < 1e10 && n == static_cast<long long>(n)) {
                return eval_gegenbauer_l(static_cast<long long>(n), alpha, x);
            }
        }

        const double d = cephes::Gamma(n + 2.0 * alpha) / cephes::Gamma(1.0 + n) / cephes::Gamma(2.0 * alpha);
        const double a = -n;
        const double b = n + 2.0 * alpha;
        const double c = alpha + 0.5;
        const T g = (1.0 - x) / 2.0;
        return d * hyp2f1(a, b, c, g);
    }

    template <typename T>
    XSF_HOST_DEVICE inline T eval_legendre(double n, T x) {
        if constexpr (cxx::is_same_v<T, double>) {
            if (cxx::isnan(x)) {
                return x;
            }
        }
        return hyp2f1(-n, n + 1.0, 1.0, 0.5 * (1.0 - x));
    }

    // Integer-degree P_n(x) by forward recurrence on the values p_k = P_k(x) and their
    // differences d_k = p_k - p_{k-1}:
    //   d_{k+1} = ((2k + 1) (x - 1) p_k + k d_k) / (k + 1),  p_{k+1} = p_k + d_{k+1},
    // starting from p_1 = x, d_1 = x - 1. Near x = 0 an explicit power series is used instead.
    template <typename Int>
    XSF_HOST_DEVICE inline double eval_legendre_l(Int n, double x) {
        if (cxx::isnan(x)) {
            return x;
        }

        if (n < 0) {
            // P_{-n-1}(x) = P_n(x); add before negating to avoid signed overflow.
            n = -(n + 1);
        }

        if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return x;
        } else if (cxx::abs(x) < 1e-5) {
            // Power series rather than recurrence due to loss of precision
            // http://functions.wolfram.com/Polynomials/LegendreP/02/
            const Int a = n / 2;
            double d = a % 2 == 0 ? 1.0 : -1.0;
            if (n % 2 == 0) {
                d *= -2.0 / beta(a + 1.0, -0.5);
            } else {
                d *= 2.0 * x / beta(a + 1.0, 0.5);
            }

            double p = 0.0;
            for (Int kk = 0; kk <= a; ++kk) {
                p += d;
                d *= -2.0 * x * x * (a - kk) * (2.0 * n + 1.0 - 2.0 * a + 2.0 * kk) /
                     ((n + 1.0 - 2.0 * a + 2.0 * kk) * (n + 2.0 - 2.0 * a + 2.0 * kk));
                if (cxx::abs(d) < 1e-20 * cxx::abs(p)) {
                    // converged
                    break;
                }
            }
            return p;
        }

        double d = x - 1.0;
        double p = x;
        for (Int kk = 0; kk < n - 1; ++kk) {
            const double k = kk + 1.0;
            d = ((2.0 * k + 1.0) / (k + 1.0)) * (x - 1.0) * p + (k / (k + 1.0)) * d;
            p = d + p;
        }
        return p;
    }

    // Chebyshev 1st kind (T)

    template <typename T>
    XSF_HOST_DEVICE inline T eval_chebyt(double n, T x) {
        return hyp2f1(-n, n, 0.5, 0.5 * (1.0 - x));
    }

    // Use the Chebyshev T recurrence directly, see Mason & Handscomb, Chebyshev Polynomials, CRC Press (2003).
    template <typename Int>
    XSF_HOST_DEVICE inline double eval_chebyt_l(Int k, double x) {
        if (k < 0) {
            // symmetry: T_{-k}(x) = T_k(x)
            k = -k;
        }

        double b2 = 0.0;
        double b1 = -1.0;
        double b0 = 0.0;
        x = 2.0 * x;
        for (Int m = 0; m < k + 1; ++m) {
            b2 = b1;
            b1 = b0;
            b0 = x * b1 - b2;
        }
        return (b0 - b2) / 2.0;
    }

    // Chebyshev 2nd kind (U)

    template <typename T>
    XSF_HOST_DEVICE inline T eval_chebyu(double n, T x) {
        return (n + 1.0) * hyp2f1(-n, n + 2.0, 1.5, 0.5 * (1.0 - x));
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_chebyu_l(Int k, double x) {
        int sign = 1;
        if (k < 0) {
            if (k == -1) {
                return 0.0;
            }
            // symmetry: U_{-k-2}(x) = -U_k(x)
            k = -k - 2;
            sign = -1;
        }

        double b2 = 0.0;
        double b1 = -1.0;
        double b0 = 0.0;
        x = 2.0 * x;
        for (Int m = 0; m < k + 1; ++m) {
            b2 = b1;
            b1 = b0;
            b0 = x * b1 - b2;
        }
        return b0 * sign;
    }

    // Chebyshev S

    template <typename T>
    XSF_HOST_DEVICE inline T eval_chebys(double n, T x) {
        return eval_chebyu(n, 0.5 * x);
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_chebys_l(Int n, double x) {
        return eval_chebyu_l(n, 0.5 * x);
    }

    // Chebyshev C

    template <typename T>
    XSF_HOST_DEVICE inline T eval_chebyc(double n, T x) {
        return 2.0 * eval_chebyt(n, 0.5 * x);
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_chebyc_l(Int n, double x) {
        return 2.0 * eval_chebyt_l(n, 0.5 * x);
    }

    // Chebyshev 1st kind shifted

    template <typename T>
    XSF_HOST_DEVICE inline T eval_sh_chebyt(double n, T x) {
        return eval_chebyt(n, 2.0 * x - 1.0);
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_sh_chebyt_l(Int n, double x) {
        return eval_chebyt_l(n, 2.0 * x - 1.0);
    }

    // Chebyshev 2nd kind shifted

    template <typename T>
    XSF_HOST_DEVICE inline T eval_sh_chebyu(double n, T x) {
        return eval_chebyu(n, 2.0 * x - 1.0);
    }

    template <typename Int>
    XSF_HOST_DEVICE inline double eval_sh_chebyu_l(Int n, double x) {
        return eval_chebyu_l(n, 2.0 * x - 1.0);
    }

} // namespace detail

// Jacobi

XSF_HOST_DEVICE inline double eval_jacobi(double n, double alpha, double beta, double x) {
    return detail::eval_jacobi(n, alpha, beta, x);
}

XSF_HOST_DEVICE inline float eval_jacobi(float n, float alpha, float beta, float x) {
    return detail::eval_jacobi(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x)
    );
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_jacobi(double n, double alpha, double beta, cxx::complex<double> x) {
    return detail::eval_jacobi(n, alpha, beta, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_jacobi(float n, float alpha, float beta, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(detail::eval_jacobi(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(beta),
        static_cast<cxx::complex<double>>(x)
    ));
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_jacobi(Int n, double alpha, double beta, double x) {
    return detail::eval_jacobi_l(n, alpha, beta, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_jacobi(Int n, float alpha, float beta, float x) {
    return detail::eval_jacobi_l(n, static_cast<double>(alpha), static_cast<double>(beta), static_cast<double>(x));
}

// Shifted Jacobi

XSF_HOST_DEVICE inline double eval_sh_jacobi(double n, double p, double q, double x) {
    return detail::eval_jacobi(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_sh_jacobi(Int n, double p, double q, double x) {
    return detail::eval_jacobi_l(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

XSF_HOST_DEVICE inline float eval_sh_jacobi(float n, float p, float q, float x) {
    return detail::eval_jacobi(
               static_cast<double>(n), static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
               2.0 * static_cast<double>(x) - 1.0
           ) /
           binom(2.0 * static_cast<double>(n) + static_cast<double>(p) - 1.0, static_cast<double>(n));
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_sh_jacobi(Int n, float p, float q, float x) {
    return detail::eval_jacobi_l(
               n, static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
               2.0 * static_cast<double>(x) - 1.0
           ) /
           binom(2.0 * n + static_cast<double>(p) - 1.0, n);
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_sh_jacobi(double n, double p, double q, cxx::complex<double> x) {
    return detail::eval_jacobi(n, p - q, q - 1.0, 2.0 * x - 1.0) / binom(2.0 * n + p - 1.0, n);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_sh_jacobi(float n, float p, float q, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_jacobi(
            static_cast<double>(n), static_cast<double>(p) - static_cast<double>(q), static_cast<double>(q) - 1.0,
            2.0 * static_cast<cxx::complex<double>>(x) - 1.0
        ) /
        binom(2.0 * static_cast<double>(n) + static_cast<double>(p) - 1.0, static_cast<double>(n))
    );
}

// Hermite (probabilist's)

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_hermitenorm(Int n, double x) {
    return detail::eval_hermitenorm(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_hermitenorm(Int n, float x) {
    return detail::eval_hermitenorm(n, static_cast<double>(x));
}

// Hermite (physicist's)

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_hermite(Int n, double x) {
    return detail::eval_hermite(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_hermite(Int n, float x) {
    return detail::eval_hermite(n, static_cast<double>(x));
}

// Generalized Laguerre. Floating-degree and complex overloads are in cpu/orthogonal_eval.h.

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_genlaguerre(Int n, double alpha, double x) {
    return detail::eval_genlaguerre_l(n, alpha, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_genlaguerre(Int n, float alpha, float x) {
    return detail::eval_genlaguerre_l(n, static_cast<double>(alpha), static_cast<double>(x));
}

// Laguerre

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_laguerre(Int n, double x) {
    return detail::eval_genlaguerre_l(n, 0.0, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_laguerre(Int n, float x) {
    return detail::eval_genlaguerre_l(n, 0.0, static_cast<double>(x));
}

// Gegenbauer (ultraspherical)

XSF_HOST_DEVICE inline double eval_gegenbauer(double n, double alpha, double x) {
    return detail::eval_gegenbauer(n, alpha, x);
}

XSF_HOST_DEVICE inline float eval_gegenbauer(float n, float alpha, float x) {
    return detail::eval_gegenbauer(static_cast<double>(n), static_cast<double>(alpha), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_gegenbauer(double n, double alpha, cxx::complex<double> x) {
    return detail::eval_gegenbauer(n, alpha, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_gegenbauer(float n, float alpha, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(detail::eval_gegenbauer(
        static_cast<double>(n), static_cast<double>(alpha), static_cast<cxx::complex<double>>(x)
    ));
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_gegenbauer(Int n, double alpha, double x) {
    return detail::eval_gegenbauer_l(n, alpha, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_gegenbauer(Int n, float alpha, float x) {
    return detail::eval_gegenbauer_l(n, static_cast<double>(alpha), static_cast<double>(x));
}

// Legendre

XSF_HOST_DEVICE inline double eval_legendre(double n, double x) { return detail::eval_legendre(n, x); }

XSF_HOST_DEVICE inline float eval_legendre(float n, float x) {
    return detail::eval_legendre(static_cast<double>(n), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_legendre(double n, cxx::complex<double> x) {
    return detail::eval_legendre(n, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_legendre(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_legendre(static_cast<double>(n), static_cast<cxx::complex<double>>(x))
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_legendre(Int n, double x) {
    return detail::eval_legendre_l(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_legendre(Int n, float x) {
    return detail::eval_legendre_l(n, static_cast<double>(x));
}

// Shifted Legendre

XSF_HOST_DEVICE inline double eval_sh_legendre(double n, double x) { return detail::eval_legendre(n, 2.0 * x - 1.0); }

XSF_HOST_DEVICE inline float eval_sh_legendre(float n, float x) {
    return detail::eval_legendre(static_cast<double>(n), 2.0 * static_cast<double>(x) - 1.0);
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_sh_legendre(double n, cxx::complex<double> x) {
    return detail::eval_legendre(n, 2.0 * x - 1.0);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_sh_legendre(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_legendre(static_cast<double>(n), 2.0 * static_cast<cxx::complex<double>>(x) - 1.0)
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_sh_legendre(Int n, double x) {
    return detail::eval_legendre_l(n, 2.0 * x - 1.0);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_sh_legendre(Int n, float x) {
    return detail::eval_legendre_l(n, 2.0 * static_cast<double>(x) - 1.0);
}

// Chebyshev 1st kind (T)

XSF_HOST_DEVICE inline double eval_chebyt(double n, double x) { return detail::eval_chebyt(n, x); }

XSF_HOST_DEVICE inline float eval_chebyt(float n, float x) {
    return detail::eval_chebyt(static_cast<double>(n), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_chebyt(double n, cxx::complex<double> x) {
    return detail::eval_chebyt(n, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_chebyt(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_chebyt(static_cast<double>(n), static_cast<cxx::complex<double>>(x))
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_chebyt(Int n, double x) {
    return detail::eval_chebyt_l(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_chebyt(Int n, float x) {
    return detail::eval_chebyt_l(n, static_cast<double>(x));
}

// Chebyshev 2nd kind (U)

XSF_HOST_DEVICE inline double eval_chebyu(double n, double x) { return detail::eval_chebyu(n, x); }

XSF_HOST_DEVICE inline float eval_chebyu(float n, float x) {
    return detail::eval_chebyu(static_cast<double>(n), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_chebyu(double n, cxx::complex<double> x) {
    return detail::eval_chebyu(n, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_chebyu(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_chebyu(static_cast<double>(n), static_cast<cxx::complex<double>>(x))
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_chebyu(Int n, double x) {
    return detail::eval_chebyu_l(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_chebyu(Int n, float x) {
    return detail::eval_chebyu_l(n, static_cast<double>(x));
}

// Chebyshev S

XSF_HOST_DEVICE inline double eval_chebys(double n, double x) { return detail::eval_chebys(n, x); }

XSF_HOST_DEVICE inline float eval_chebys(float n, float x) {
    return detail::eval_chebys(static_cast<double>(n), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_chebys(double n, cxx::complex<double> x) {
    return detail::eval_chebys(n, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_chebys(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_chebys(static_cast<double>(n), static_cast<cxx::complex<double>>(x))
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_chebys(Int n, double x) {
    return detail::eval_chebys_l(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_chebys(Int n, float x) {
    return detail::eval_chebys_l(n, static_cast<double>(x));
}

// Chebyshev C

XSF_HOST_DEVICE inline double eval_chebyc(double n, double x) { return detail::eval_chebyc(n, x); }

XSF_HOST_DEVICE inline float eval_chebyc(float n, float x) {
    return detail::eval_chebyc(static_cast<double>(n), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_chebyc(double n, cxx::complex<double> x) {
    return detail::eval_chebyc(n, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_chebyc(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_chebyc(static_cast<double>(n), static_cast<cxx::complex<double>>(x))
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_chebyc(Int n, double x) {
    return detail::eval_chebyc_l(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_chebyc(Int n, float x) {
    return detail::eval_chebyc_l(n, static_cast<double>(x));
}

// Chebyshev 1st kind shifted

XSF_HOST_DEVICE inline double eval_sh_chebyt(double n, double x) { return detail::eval_sh_chebyt(n, x); }

XSF_HOST_DEVICE inline float eval_sh_chebyt(float n, float x) {
    return detail::eval_sh_chebyt(static_cast<double>(n), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_sh_chebyt(double n, cxx::complex<double> x) {
    return detail::eval_sh_chebyt(n, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_sh_chebyt(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_sh_chebyt(static_cast<double>(n), static_cast<cxx::complex<double>>(x))
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_sh_chebyt(Int n, double x) {
    return detail::eval_sh_chebyt_l(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_sh_chebyt(Int n, float x) {
    return detail::eval_sh_chebyt_l(n, static_cast<double>(x));
}

// Chebyshev 2nd kind shifted

XSF_HOST_DEVICE inline double eval_sh_chebyu(double n, double x) { return detail::eval_sh_chebyu(n, x); }

XSF_HOST_DEVICE inline float eval_sh_chebyu(float n, float x) {
    return detail::eval_sh_chebyu(static_cast<double>(n), static_cast<double>(x));
}

XSF_HOST_DEVICE inline cxx::complex<double> eval_sh_chebyu(double n, cxx::complex<double> x) {
    return detail::eval_sh_chebyu(n, x);
}

XSF_HOST_DEVICE inline cxx::complex<float> eval_sh_chebyu(float n, cxx::complex<float> x) {
    return static_cast<cxx::complex<float>>(
        detail::eval_sh_chebyu(static_cast<double>(n), static_cast<cxx::complex<double>>(x))
    );
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline double eval_sh_chebyu(Int n, double x) {
    return detail::eval_sh_chebyu_l(n, x);
}

template <typename Int, cxx::enable_if_t<cxx::is_integral_v<Int>, int> = 0>
XSF_HOST_DEVICE inline float eval_sh_chebyu(Int n, float x) {
    return detail::eval_sh_chebyu_l(n, static_cast<double>(x));
}

} // namespace xsf
