#pragma once

#include "xsf/bessel.h"
#include "xsf/binom.h"
#include "xsf/cephes/bdtr.h"
#include "xsf/cephes/chdtr.h"
#include "xsf/cephes/fdtr.h"
#include "xsf/cephes/gdtr.h"
#include "xsf/cephes/incbet.h"
#include "xsf/cephes/incbi.h"
#include "xsf/cephes/kolmogorov.h"
#include "xsf/cephes/nbdtr.h"
#include "xsf/cephes/ndtr.h"
#include "xsf/cephes/ndtri.h"
#include "xsf/cephes/owens_t.h"
#include "xsf/cephes/pdtr.h"
#include "xsf/cephes/tukey.h"
#include "xsf/erf.h"
#include "xsf/gamma.h"
#include <numeric>

namespace xsf {

inline double bdtr(double k, int n, double p) { return cephes::bdtr(k, n, p); }

inline double bdtri(double k, int n, double y) { return cephes::bdtri(k, n, y); }

inline double bdtrc(double k, int n, double p) { return cephes::bdtrc(k, n, p); }

inline double chdtr(double df, double x) { return cephes::chdtr(df, x); }

inline float chdtr(float df, float x) { return static_cast<float>(cephes::chdtr(df, x)); }

inline double chdtrc(double df, double x) { return cephes::chdtrc(df, x); }

inline float chdtrc(float df, float x) { return static_cast<float>(cephes::chdtrc(df, x)); }

inline double chdtri(double df, double y) { return cephes::chdtri(df, y); }

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

inline double fdtr(double a, double b, double x) { return cephes::fdtr(a, b, x); }

inline double fdtrc(double a, double b, double x) { return cephes::fdtrc(a, b, x); }

inline double fdtri(double a, double b, double y) { return cephes::fdtri(a, b, y); }

inline double gdtr(double a, double b, double x) { return cephes::gdtr(a, b, x); }

inline double gdtrc(double a, double b, double x) { return cephes::gdtrc(a, b, x); }

inline double kolmogorov(double x) { return cephes::kolmogorov(x); }

inline double kolmogc(double x) { return cephes::kolmogc(x); }

inline double kolmogi(double x) { return cephes::kolmogi(x); }

inline double kolmogci(double x) { return cephes::kolmogci(x); }

inline double kolmogp(double x) { return cephes::kolmogp(x); }

inline double ndtr(double x) { return cephes::ndtr(x); }

inline float ndtr(float x) { return ndtr(static_cast<double>(x)); }

inline std::complex<double> ndtr(std::complex<double> z) { return 0.5 * erfc(-z * M_SQRT1_2); }

inline std::complex<float> ndtr(std::complex<float> z) {
    return static_cast<std::complex<float>>(ndtr(static_cast<std::complex<double>>(z)));
}

/*
 * Log of the CDF of the normal distribution for double x.
 *
 * Let F(x) be the CDF of the standard normal distribution.
 * This implementation of log(F(x)) is based on the identities
 *
 *   F(x) = erfc(-x/√2)/2
 *        = 1 - erfc(x/√2)/2
 *
 * We use the first formula for x < -1, with erfc(z) replaced
 * by erfcx(z)*exp(-z**2) to ensure high precision for large
 * negative values when we take the logarithm:
 *
 *   log F(x) = log(erfc(-x/√2)/2)
 *            = log(erfcx(-x/√2)/2)*exp(-x**2/2))
 *            = log(erfcx(-x/√2)/2) - x**2/2
 *
 * For x >= -1, we use the second formula for F(x):
 *
 *   log F(x) = log(1 - erfc(x/√2)/2)
 *            = log1p(-erfc(x/√2)/2)
 */
inline double log_ndtr(double x) {
    double t = x * M_SQRT1_2;
    if (x < -1.0) {
        return log(erfcx(-t) / 2) - t * t;
    } else {
        return log1p(-erfc(t) / 2);
    }
}

inline float log_ndtr(float x) { return log_ndtr(static_cast<double>(x)); }

/*
 * Log of the normal CDF for complex arguments.
 *
 * This is equivalent to log(ndtr(z)), but is more robust to overflow at $z\to\infty$.
 * This implementation uses $\erfc(z) = \exp(-z^2) w(iz)$ taking special care to select
 * the principal branch of the log function log( exp(-z^2) w(i z) )
 */
inline std::complex<double> log_ndtr(std::complex<double> z) {
    if (z.real() > 6) {
        // Underflow. Close to the real axis, expand the log in log(1 - ndtr(-z)).
        std::complex<double> w = -0.5 * erfc(z * M_SQRT1_2);
        if (std::abs(w) < 1e-8) {
            return w;
        }
    }

    z *= -M_SQRT1_2;
    double x = std::real(z);
    double y = std::imag(z);

    /* Compute the principal branch of $log(exp(-z^2))$, using the fact that
     * $log(e^t) = log|e^t| + i Arg(e^t)$, and that if $t = r + is$, then
     * $e^t = e^r (\cos(s) + i \sin(s))$.
     */
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2 * x * y;        // Im(-z^2)

    double im = fmod(mIm_z2, 2.0 * M_PI);
    if (im > M_PI) {
        im -= 2.0 * M_PI;
    }

    std::complex<double> val1 = std::complex<double>(mRe_z2, im);

    std::complex<double> val2 = log(xsf::wofz(complex<double>(-y, x)));
    std::complex<double> result = val1 + val2 - M_LN2;

    /* Again, select the principal branch: log(z) = log|z| + i arg(z), thus
     * the imaginary part of the result should belong to [-pi, pi].
     */
    im = imag(result);
    if (im >= M_PI) {
        im -= 2 * M_PI;
    }
    if (im < -M_PI) {
        im += 2 * M_PI;
    }

    return {result.real(), im};
}

inline std::complex<float> log_ndtr(std::complex<float> z) {
    return static_cast<std::complex<float>>(log_ndtr(static_cast<std::complex<double>>(z)));
}

inline double nbdtr(int k, int n, double p) { return cephes::nbdtr(k, n, p); }

inline double nbdtrc(int k, int n, double p) { return cephes::nbdtrc(k, n, p); }

inline double nbdtri(int k, int n, double p) { return cephes::nbdtri(k, n, p); }

inline double ndtri(double x) { return cephes::ndtri(x); }

inline float ndtri(float x) { return static_cast<float>(cephes::ndtri(x)); }

inline double owens_t(double h, double a) { return cephes::owens_t(h, a); }

inline double pdtr(double k, double m) { return cephes::pdtr(k, m); }

inline double pdtrc(double k, double m) { return cephes::pdtrc(k, m); }

inline double pdtri(int k, double y) { return cephes::pdtri(k, y); }

/*
 * Compute the exact p-value of the Cramer-von Mises two-sample test
 * for a given value s of the test statistic.
 * m and n are the sizes of the samples.
 *
 * [1] Y. Xiao, A. Gordon, and A. Yakovlev, "A C++ Program for
 *     the Cramér-Von Mises Two-Sample Test", J. Stat. Soft.,
 *     vol. 17, no. 8, pp. 1-15, Dec. 2006.
 * [2] T. W. Anderson "On the Distribution of the Two-Sample Cramer-von Mises
 *     Criterion," The Annals of Mathematical Statistics, Ann. Math. Statist.
 *     33(3), 1148-1159, (September, 1962)
 */
inline double pval_cvm_2samp_exact(double s, int64_t m, int64_t n) {
    // [1, p. 3]
    int64_t lcm = std::lcm(m, n);
    // [1, p. 4], below eq. 3
    int64_t a = lcm / m;
    int64_t b = lcm / n;
    // Combine Eq. 9 in [2] with Eq. 2 in [1] and solve for $\zeta$
    // Hint: `s` is $U$ in [2], and $T_2$ in [1] is $T$ in [2]
    int64_t mn = m * n;
    // Uses double floor division since s is double
    int64_t zeta = std::floor((lcm * lcm * (m + n) * (6.0 * s - mn * (4.0 * mn - 1))) / (6.0 * mn * mn));
    int64_t combinations = static_cast<int64_t>(xsf::binom(static_cast<double>(m + n), static_cast<double>(m)));
    // Each frequency table maps value -> frequency,
    // mirroring the 2-row numpy array where row 0 = values, row 1 = frequencies
    using FreqTable = std::map<int64_t, int64_t>;
    // the frequency table of g_{u, v}^+ defined in [1, p. 6]
    // gs[0] = {0: 1}, gs[1..m] = empty
    std::vector<FreqTable> gs(m + 1);
    gs[0][0] = 1;
    for (int64_t u = 0; u < n + 1; ++u) {
        std::vector<FreqTable> next_gs;
        FreqTable tmp;
        for (int64_t v = 0; v < m + 1; ++v) {
            // Calculate g recursively with eq. 11 in [1]. Even though it
            // doesn't look like it, this also does 12/13 (all of Algorithm 1).
            const FreqTable &g = gs[v];
            // Merge tmp and g: for common keys sum frequencies,
            // keep unique keys from both sides.
            // (equivalent to np.intersect1d + concatenate logic)
            FreqTable merged;
            for (const auto &[key, freq] : tmp) {
                merged[key] += freq;
            }
            for (const auto &[key, freq] : g) {
                merged[key] += freq;
            }
            int64_t diff = a * v - b * u;
            int64_t res = diff * diff;
            // tmp[0] += res  (shift all keys by res)
            tmp.clear();
            for (const auto &[key, freq] : merged) {
                tmp[key + res] += freq;
            }
            next_gs.push_back(tmp);
        }
        gs = std::move(next_gs);
    }
    // (equivalent to return np.float64(np.sum(freq[value >= zeta]) / combinations))
    const FreqTable &final_table = gs[m];
    int64_t sum_freq = 0;
    for (const auto &[value, freq] : final_table) {
        if (value >= zeta) {
            sum_freq += freq;
        }
    }
    return static_cast<double>(sum_freq) / static_cast<double>(combinations);
}

inline double smirnov(int n, double x) { return cephes::smirnov(n, x); }

inline double smirnovc(int n, double x) { return cephes::smirnovc(n, x); }

inline double smirnovi(int n, double x) { return cephes::smirnovi(n, x); }

inline double smirnovci(int n, double x) { return cephes::smirnovci(n, x); }

inline double smirnovp(int n, double x) { return cephes::smirnovp(n, x); }

inline double tukeylambdacdf(double x, double lmbda) { return cephes::tukeylambdacdf(x, lmbda); }

inline float tukeylambdacdf(float x, double lmbda) {
    return tukeylambdacdf(static_cast<double>(x), static_cast<double>(lmbda));
}

template <typename InputMat, typename OutputMat>
XSF_HOST_DEVICE inline void poisson_binom_pmf_all(InputMat p, OutputMat res) {
    /* Outputs the full pmf of a Poisson-Binomial distribution
     *
     * p should be an std::mdspan view of a 1d array of length n containing
     * the success probabilities for the n Bernoulli trials. res should be
     * a std::mdspan view of a 1d array of zeros of length n + 1. It is up to
     * the caller to pass valid values for p and res.
     *
     * Upon completion, res[k] will contain the probability of observing k
     * successes for k from 0 to n.
     */
    using T = typename OutputMat::value_type;
    auto n = p.extent(0);
    auto out_size = res.extent(0);

    if (out_size != n + 1) {
        set_error("poisson_binom_pmf", SF_ERROR_MEMORY, "out.shape[-1] must be p.shape[-1] + 1");
    }

    if (n == 0) {
        res[0] = T(1.0);
        return;
    }

    res[0] = T(1.0) - p[0];
    res[1] = p[0];

    for (decltype(n) i = 1; i < n; i++) {
        T p_i = p[i];
        T q_i = 1 - p_i;
        for (decltype(n) j = i + 1; j >= 1; j--) {
            T tmp = res[j - 1] * p_i;
            res[j - 1] *= q_i;
            res[j] += tmp;
        }
    }
}

} // namespace xsf
