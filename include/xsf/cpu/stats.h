#include "../bessel.h"
#include "../cephes/kolmogorov.h"
#include "../config.h"
#include "../gamma.h"

namespace xsf {
namespace cpu {

    inline double kolmogorov(double x) { return cephes::kolmogorov(x); }

    inline float kolmogorov(float x) { return static_cast<float>(kolmogorov(static_cast<double>(x))); }

    inline double kolmogc(double x) { return cephes::kolmogc(x); }

    inline float kolmogc(float x) { return static_cast<float>(kolmogc(static_cast<double>(x))); }

    inline double kolmogi(double x) { return cephes::kolmogi(x); }

    inline float kolmogi(float x) { return static_cast<float>(kolmogi(static_cast<double>(x))); }

    inline double kolmogci(double x) { return cephes::kolmogci(x); }

    inline float kolmogci(float x) { return static_cast<float>(kolmogci(static_cast<double>(x))); }

    inline double kolmogp(double x) { return cephes::kolmogp(x); }

    inline float kolmogp(float x) { return static_cast<float>(kolmogp(static_cast<double>(x))); }

    inline double smirnov(int n, double x) { return cephes::smirnov(n, x); }

    inline float smirnov(int n, float x) { return static_cast<float>(smirnov(n, static_cast<double>(x))); }

    inline double smirnovc(int n, double x) { return cephes::smirnovc(n, x); }

    inline float smirnovc(int n, float x) { return static_cast<float>(smirnovc(n, static_cast<double>(x))); }

    inline double smirnovi(int n, double x) { return cephes::smirnovi(n, x); }

    inline float smirnovi(int n, float x) { return static_cast<float>(smirnovi(n, static_cast<double>(x))); }

    inline double smirnovci(int n, double x) { return cephes::smirnovci(n, x); }

    inline float smirnovci(int n, float x) { return static_cast<float>(smirnovci(n, static_cast<double>(x))); }

    inline double smirnovp(int n, double x) { return cephes::smirnovp(n, x); }

    inline float smirnovp(int n, float x) { return static_cast<float>(smirnovp(n, static_cast<double>(x))); }

    inline double cdf_cvm_inf(double x) {
        // CDF of the Cramer-von Mises test statistic (infinite sample limit).
        // Accurate for practical hypothesis testing but not expected to be accurate
        // for large values of x, e.g. x > 4, when the cdf is very close to 1.
        if (std::isnan(x)) {
            return x;
        }
        if (x <= 0) {
            return 0.0;
        }

        double total = 0.0;
        int k = 0;

        while (true) {
            double u = std::exp(gammaln(k + 0.5) - gammaln(k + 1.0)) / (std::pow(M_PI, 1.5) * std::sqrt(x));
            double y = 4.0 * k + 1.0;
            double q = (y * y) / (16.0 * x);
            double b = cyl_bessel_k(0.25, q);
            double z = u * std::sqrt(y) * std::exp(-q) * b;

            total += z;
            if (std::abs(z) < 1e-7) {
                break;
            }
            ++k;
        }

        return total;
    }

    // Helpers for psi1 computations
    inline double psi1_G(double y) {
        // Helper implementing the G-type Bessel combination appearing in equation 1.10
        // in Csörgő, S. and Faraway, J. (1996).
        double z = y * y / 4.0;
        double b = cyl_bessel_k(0.25, z) + cyl_bessel_k(0.75, z);
        return std::exp(-z) * std::pow(y / 2.0, 1.5) * b / std::sqrt(M_PI);
    }

    inline double psi1_H(double y) {
        // Helper implementing the H-type Bessel combination appearing in equation 1.10
        // in Csörgő, S. and Faraway, J. (1996).
        double z = y * y / 4.0;
        double b = 2.0 * cyl_bessel_k(0.25, z) + 3.0 * cyl_bessel_k(0.75, z) - cyl_bessel_k(1.25, z);
        return std::exp(-z) * std::pow(y / 2.0, 2.5) * b / std::sqrt(M_PI);
    }

    inline double psi1_term(int k, double x) {
        // Compute the k-th term of the psi1 series expansion in equation 1.10
        // of Csörgő, S. and Faraway, J. (1996).
        double m = 2.0 * k + 1.0;
        double sx = 2.0 * std::sqrt(x);
        double y1 = std::pow(x, 0.75);
        double y2 = std::pow(x, 1.25);
        double gamma_kp1_2 = gamma(k + 0.5);
        double gamma_kp3_2 = gamma(k + 1.5);
        double e1 = m * gamma_kp1_2 * psi1_G((4.0 * k + 3.0) / sx) / (9.0 * y1);
        double e2 = gamma_kp1_2 * psi1_H((4.0 * k + 1.0) / sx) / (72.0 * y2);
        double e3 = 2.0 * (m + 2.0) * gamma_kp3_2 * psi1_H((4.0 * k + 5.0) / sx) / (12.0 * y2);
        double e4 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k + 1.0) / sx) / (144.0 * y1);
        double e5 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k + 5.0) / sx) / (144.0 * y1);
        return e1 + e2 + e3 + e4 + e5;
    }

    inline double psi1_mod(double x) {
        // psi1 is defined in equation 1.10 in Csörgő, S. and Faraway, J. (1996).
        // This implements a modified version by excluding the term V(x) / 12
        // (here: cdf_cvm_inf(x) / 12) to avoid evaluating cdf_cvm_inf(x) twice in cdf_cvm.
        //
        // Implementation based on MAPLE code of Julian Faraway and R code of the
        // function pCvM in the package goftest (v1.1.1), permission granted
        // by Adrian Baddeley. Main difference in the implementation: the code
        // here keeps adding terms of the series until the terms are small enough.
        int k = 0;
        double tot = 0.0;
        while (true) {
            double gamma_kp1 = gamma(k + 1.0);
            double z = -psi1_term(k, x) / (M_PI * gamma_kp1);
            tot += z;
            if (std::abs(z) < 1e-7) {
                break;
            }
            ++k;
        }
        return tot;
    }

    inline double cdf_cvm(double x, int n = -1) {
        // Calculate the CDF of the Cramér-von Mises statistic for a finite sample
        // size n. If n=-1, use the asymptotic CDF (n=inf).
        //
        // See equation 1.8 in Csörgő, S. and Faraway, J. (1996) for finite samples,
        // 1.2 for the asymptotic CDF.
        //
        // For finite n, the approximation is less accurate near the support boundaries
        // [1/(12*n), n/3] and for larger x values where the CDF is close to 1.
        // The returned value is clipped to [0, 1].
        if (n == -1) {
            return cdf_cvm_inf(x);
        }
        double lower = 1.0 / (12.0 * n);
        double upper = static_cast<double>(n) / 3.0;
        if (x <= lower) {
            return 0.0;
        }
        if (x >= upper) {
            return 1.0;
        }
        double cdf = cdf_cvm_inf(x) * (1.0 + lower) + psi1_mod(x) / n;
        return std::min(std::max(cdf, 0.0), 1.0);
    }

    // Logarithm of the generalized inverse Gaussian probability density function.
    inline double geninvgauss_logpdf(double x, double p, double b) {
        if (x <= 0) {
            return -std::numeric_limits<double>::infinity();
        }

        double z = cyl_bessel_ke(p, b);
        if (std::isinf(z)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return -std::log(2.0) - std::log(z) + b + (p - 1.0) * std::log(x) - b * (x + 1.0 / x) / 2.0;
    }

    inline float geninvgauss_logpdf(float x, float p, float b) {
        return static_cast<float>(
            geninvgauss_logpdf(static_cast<double>(x), static_cast<double>(p), static_cast<double>(b))
        );
    }

    // Generalized inverse Gaussian probability density function.
    inline double geninvgauss_pdf(double x, double p, double b) { return std::exp(geninvgauss_logpdf(x, p, b)); }

    inline float geninvgauss_pdf(float x, float p, float b) {
        return static_cast<float>(
            geninvgauss_pdf(static_cast<double>(x), static_cast<double>(p), static_cast<double>(b))
        );
    }

    namespace detail {

        inline double genhyperbolic_log_norming_constant(double p, double a, double b) {
            double t1, t2, t3, t4, t5, t6;

            t1 = (a + b) * (a - b);
            t2 = p * 0.5 * std::log(t1);
            t3 = 0.5 * std::log(2 * M_PI);
            t4 = (p - 0.5) * std::log(a);
            t5 = std::sqrt(t1);
            t6 = std::log(cyl_bessel_ke(p, t5)) - t5;

            return t2 - t3 - t4 - t6;
        }

    } // namespace detail

    // Logarithm of the generalized hyperbolic probability density function.
    inline double genhyperbolic_logpdf(double x, double p, double a, double b) {
        double t1, t2, t3, t4, t5;

        t1 = detail::genhyperbolic_log_norming_constant(p, a, b);
        t2 = std::sqrt(1.0 + x * x);
        t3 = (p - 0.5) * std::log(t2);
        t4 = std::log(cyl_bessel_ke(p - 0.5, a * t2)) - a * t2;
        t5 = b * x;

        return t1 + t3 + t4 + t5;
    }

    inline float genhyperbolic_logpdf(float x, float p, float a, float b) {
        return static_cast<float>(genhyperbolic_logpdf(
            static_cast<double>(x), static_cast<double>(p), static_cast<double>(a), static_cast<double>(b)
        ));
    }

    // Generalized hyperbolic probability density function.
    inline double genhyperbolic_pdf(double x, double p, double a, double b) {
        return std::exp(genhyperbolic_logpdf(x, p, a, b));
    }

    inline float genhyperbolic_pdf(float x, float p, float a, float b) {
        return static_cast<float>(genhyperbolic_pdf(
            static_cast<double>(x), static_cast<double>(p), static_cast<double>(a), static_cast<double>(b)
        ));
    }

} // namespace cpu
} // namespace xsf
