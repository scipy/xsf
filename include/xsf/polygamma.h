/*
 * Calculation of the polygamma function for positive integer orders and real and complex inputs
 * Author: Lorenzo Peri
 */

#pragma once

#include "cephes/psi.h"
#include "cephes/zeta.h"
#include "config.h"
#include "digamma.h"
#include "error.h"
#include "gamma.h"
#include "trig.h"
#include "trigamma.h"

namespace xsf {

constexpr int MAX_ORDER_BERNOULLI_SERIES = 6;

namespace detail {

    XSF_HOST_DEVICE inline std::complex<double>
    polygamma_forward_recurrence(int n, std::complex<double> const z, std::complex<double> const psiz, int m) {
        /* Compute polygamma(n, z + m) using polygamma(n, z) using the recurrence relation
         *    polygamma(n, z + 1) = polygamma(n, z) + (-1)^n * n! / z^(n+1).
         * See https://dlmf.nist.gov/5.15#E5 */

        std::complex<double> zk, res = psiz;

        for (int k = 0; k < m; k++) {
            zk = (z + static_cast<double>(k));
            res += std::pow(zk, -(n + 1));
        }

        return res;
    }

    XSF_HOST_DEVICE inline std::complex<double>
    polygamma_backward_recurrence(int n, std::complex<double> const z, std::complex<double> const psiz, int m) {
        /* Compute polygamma(n, z + m) using polygamma(n, z) using the recurrence relation
         *    polygamma(n, z + 1) = polygamma(n, z) + (-1)^n * n! / z^(n+1).
         * See https://dlmf.nist.gov/5.15#E5 */

        std::complex<double> zk, res = psiz;

        for (int k = 1; k < m + 1; k++) {
            zk = (z - static_cast<double>(k));
            res -= std::pow(zk, -(n + 1));
        }

        return res;
    }

    XSF_HOST_DEVICE inline std::complex<double> polygamma_asymptotic_series(int n, std::complex<double> z) {
        /* Evaluate digamma using an asymptotic series. See
         * http://dlmf.nist.gov/5.15.E8
         * Higher order converge slower, so we need slightly more Bernoulli
         * numbers than the digamma (Bernoulli numbers calculated via sympy).
         */
        static constexpr std::array<double, 36> bernoulli2k = {
            1.66666666666666666670e-1,   -3.33333333333333333330e-2,  2.38095238095238095240e-2,
            -3.33333333333333333330e-2,  7.57575757575757575760e-2,   -2.53113553113553113550e-1,
            1.16666666666666666670e+0,   -7.09215686274509803920e+0,  5.49711779448621553880e+1,
            -5.29124242424242424240e+2,  6.19212318840579710140e+3,   -8.65802531135531135530e+4,
            1.42551716666666666670e+6,   -2.72982310678160919540e+7,  6.01580873900642368380e+8,
            -1.51163157670921568630e+10, 4.29614643061166666670e+11,  -1.37116552050883327720e+13,
            4.88332318973593166670e+14,  -1.92965793419400681490e+16, 8.41693047573682615000e+17,
            -4.03380718540594554130e+19, 2.11507486380819916060e+21,  -1.20866265222965259350e+23,
            7.50086674607696436690e+24,  -5.03877810148106891410e+26, 3.65287764848181233350e+28,
            -2.84987693024508822260e+30, 2.38654274996836276450e+32,  -2.13999492572253336660e+34,
            2.05009757234780975700e+36,  -2.09380059113463784090e+38, 2.27526964884635155600e+40,
            -2.62577102862395760470e+42, 3.21250821027180325180e+44,  -4.15982781667947109140e+46,
        };

        std::complex<double> res, term, fac_coeff;
        std::complex<double> const zn_inv = std::pow(z, -n);
        std::complex<double> const rzz = 1. / z / z;
        std::complex<double> zfac = zn_inv;

        if (!(std::isfinite(z.real()) && std::isfinite(z.imag()))) {
            /* Check for infinity (or nan) and return early.
             * Result of division by complex infinity is implementation dependent.
             * and has been observed to vary between C++ stdlib and CUDA stdlib.
             */
            return 0.; // Asymptotically psi(n, z) ~ 1/z, so just return zero?
        }

        res = zn_inv * (1. / n + 0.5 / z);

        fac_coeff = (0.5 * static_cast<std::complex<double>>(n + 1));

        for (int k = 1; k < bernoulli2k.size() + 1; k++) {

            zfac *= rzz;
            term = bernoulli2k[k - 1] * fac_coeff * zfac;
            res += term;

            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(res)) {
                break;
            }

            double const num_coeff = (2. * k + n) * (2. * k + n + 1.);
            double const den_coeff = (2. * k + 2) * (2. * k + 1.);
            fac_coeff *= num_coeff / den_coeff;
        }

        return -res;
    }

    XSF_HOST_DEVICE inline std::complex<double> polygamma_reflection_rhs(int n, std::complex<double> z) {
        /* The reflection formula for the polygamma of arbitrary positive order
         * reads (http://dlmf.nist.gov/5.15.E6)
         *   psi(n, 1-z) + (-1)**(n-1) * psi(n, z) = (-1)**(n-1) * pi * d^n/dz^n cot(pi z)
         * To calculate this derivative we can leverage the fact that, if y = cot(z), then
         * d/dz cot(z) = - (1 + y **2). Therefore, the derivatives of y = cot(z)) are all
         * polynomials in cot(z)
         *   P_n(y) = sum_{m=1}^{n+1} c_{n, m} y**m
         * that satisfy the recurrence relation
         *   c_{n+1, m} = -(m+1) c_{n, m+1} - (m-1) c_{n, m-1}
         * This allows fast calculation of the coefficients of P_n(y).
         * However, this requires us to loop through the calculation n times,
         * and the coefficients blow up, leading to loss of precision.
         * For larger n, there is little point using this strategy.
         */

        if (n > MAX_ORDER_BERNOULLI_SERIES) {
            // We are not going to use it for n>MAX_ORDER_BERNOULLI_SERIES.
            // Statically allocate the memory and check against misuse
            return std::numeric_limits<double>::quiet_NaN();
        }

        std::array<double, MAX_ORDER_BERNOULLI_SERIES + 2> coefficient_array{}; // all zeros at the start
        std::array<double, MAX_ORDER_BERNOULLI_SERIES + 2> next_coeffs;

        coefficient_array[1] = 1.; // P_0(y) = y

        /* The coefficients blow up very rapidly roughly O(n!). So we compute the
         * normalized RHS (-1)**(n-1) * pi * d^n/dz^n cot(pi z) / n!, and let the
         * general prefactor in the main polygamma driver handle the n! scaling
         */
        for (int k = 0; k < n; k++) {
            std::fill(next_coeffs.begin(), next_coeffs.end(), 0.);
            for (int m = 0; m < k + 3; m++) {
                double const term1 = (m + 1 <= k + 1) ? (-static_cast<double>(m + 1) * coefficient_array[m + 1]) : 0.;
                double const term2 = (m - 1 >= 0) ? (-static_cast<double>(m - 1) * coefficient_array[m - 1]) : 0.;
                next_coeffs[m] = (term1 + term2) / static_cast<double>(k + 1);
            }
            std::copy(next_coeffs.begin(), next_coeffs.end(), coefficient_array.begin());
        }

        std::complex<double> result = 0.0;
        std::complex<double> cot_piz = cospi(z) / sinpi(z);
        // Evaluate from highest degree down to 0
        for (int m = n + 1; m >= 0; m--) {
            result = result * cot_piz + coefficient_array[m];
        }

        double const sign = (n & 1) ? -1.0 : 1.0;
        return sign * std::pow(M_PI, n + 1) * result;
    }

    XSF_HOST_DEVICE inline std::complex<double> polygamma_Hurwitz_series(int n, std::complex<double> z) {
        /* Evaluate the sum for the Hurwitz zeta function:
         * zeta(n+1, z) = sum_{k=0}^\infty 1 / (z+k)^{n+1}
         * For n > 8, the terms decay exponentially fast, making direct
         * summation far more stable than the asymptotic Bernoulli series.
         */
        std::complex<double> res = 0.0;
        for (int k = 0; k < 10000; k++) {
            std::complex<double> const zk = z + static_cast<double>(k);
            std::complex<double> const term = std::pow(zk, -(n + 1));
            res += term;

            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(res)) {
                break;
            }
        }
        return -res;
    }

} // namespace detail

XSF_HOST_DEVICE inline double polygamma(int n, double z) {
    if (n < 0) {
        /* Although it is formally possible to extend the polygamma function to negative orders
         * this is not implemented yet. For more information see:
         * O. Espinosa, and V.H. Moll (2003), "A GENERALIZED POLYGAMMA FUNCTION"
         * https://arxiv.org/pdf/math.CA/0305079
         */
        set_error("polygamma", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    /* To compute the polygamma for real inputs we can leverage Cehpes via the identity
     *    psi(n, x) = (-1)**(n+1) * n! * zeta(n+1, x)
     * where n is a positive integer and zeta is the Hurwitz zeta function
     */
    if (n == 0) {
        return digamma(z);
    } else if (n == 1) {
        return trigamma(z);
    }

    /* (-1)^(n+1) */
    const double sign = n & 1 ? 1. : -1.;
    const double np1 = static_cast<double>(n + 1);

    double const factorial = gamma(np1);

    const double Hur_zeta = cephes::zeta(np1, z);

    return sign * factorial * Hur_zeta;
}

XSF_HOST_DEVICE inline float polygamma(int n, float z) {
    return static_cast<float>(polygamma(n, static_cast<double>(z)));
}

XSF_HOST_DEVICE inline std::complex<double> polygamma(int n, std::complex<double> z) {
    /*
     * Compute the polygamma function for complex arguments. The strategy is:
     * - For small orders (n <= (MAX_ORDER_BERNOULLI_SERIES = 6)):
     *      - If close to the origin, use a reflection relation to step away
     *      from the origin.
     *      - If close to the negative real axis, use the recurrence
     *      formula to move to the right halfplane.
     *      - If |z| is large (> 16), use the asymptotic series.
     *      - If |z| is small, use a recurrence relation to make |z| large
     *      enough to use the asymptotic series.
     * - For large orders (n > (MAX_ORDER_BERNOULLI_SERIES = 6)):
     *      - Use the recurrence relation to move in the positive half plane
     *      - Evaluate the Hurwitz zeta, which converges faster and to higher accuracy
     *      than the Bernoulli asymptotic series for large n.
     */
    if (n < 0) {
        /* Although it is formally possible to extend the polygamma function to negative orders
         * this is not implemented yet. For more information see:
         * O. Espinosa, and V.H. Moll (2003), "A GENERALIZED POLYGAMMA FUNCTION"
         * https://arxiv.org/pdf/math.CA/0305079
         */
        set_error("polygamma", SF_ERROR_DOMAIN, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    // Fall back to digamma and trigamma
    if (n == 0) {
        return digamma(z);
    } else if (n == 1) {
        return trigamma(z);
    }

    if (std::abs(z.imag()) < std::numeric_limits<double>::epsilon()) {
        // use the real zeta, because it evaluates to better precision
        return static_cast<std::complex<double>>(polygamma(n, z.real()));
    }

    double absz = std::abs(z);
    std::complex<double> res = 0.;
    /* The reflection formula flips the sign for odd polygammas.
     * We must keep track of that. */
    std::complex<double> sign = 1.;
    std::complex<double> const prefactor_sign = (n & 1) ? -1. : 1.;
    std::complex<double> const prefactor = prefactor_sign * gamma(static_cast<double>(n + 1));
    /* Use the asymptotic series for z away from the negative real axis
     * with abs(z) > smallabsz. */
    double const smallabsz = 16.;

    if (z.real() <= 0.0 && std::ceil(z.real()) == z) {
        // Poles
        set_error("polygamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    if (n <= MAX_ORDER_BERNOULLI_SERIES) {
        /* For small orders, the Hurwitz series converges slowly, best take the approaches
         * of the di- and trigamma with the Bernoulli series
         */
        if (z.real() < 0 && std::abs(z.imag()) < smallabsz) {
            std::complex<double> const reflection_rhs = detail::polygamma_reflection_rhs(n, z);
            res -= reflection_rhs;
            sign *= (n & 1) ? -1. : 1.;
            z = 1. - z;
            absz = std::abs(z);
        }
        if (absz < 0.5) {
            res -= sign * std::pow(z, -(n + 1));
            z += 1;
            absz = std::abs(z);
        }

        if (absz > smallabsz) {
            res += sign * detail::polygamma_asymptotic_series(n, z);
        } else if (z.real() >= 0.0) {
            double m = std::trunc(smallabsz - z.real()) + 1;
            std::complex<double> init = detail::polygamma_asymptotic_series(n, z + m);
            res += sign * detail::polygamma_backward_recurrence(n, z + m, init, m);
        } else {
            double m = std::trunc(smallabsz + z.real()) + 1;
            std::complex<double> init = detail::polygamma_asymptotic_series(n, z - m);
            res += sign * detail::polygamma_forward_recurrence(n, z - m, init, m);
        }
    } else {
        /* For large n, Hurwitz series converges converges very fast, so we can use it directly.
         * However, we need to make sure that Re[z] is negative-enough that
         * sum_{k=0}^\infty 1 / (z+k)^{n+1} converges quickly.
         * The reflection is very expensive to calculate and its coefficients blow up like n!.
         * We can however use the recurrence formula to step away towards
         * the positive real semi-plane */
        while ((absz < 0.5) || (z.real() < 0)) {
            res -= std::pow(z, -(n + 1));
            z += 1;
            absz = std::abs(z);
        }

        res += detail::polygamma_Hurwitz_series(n, z);
    }

    return prefactor * res;
}

XSF_HOST_DEVICE inline std::complex<float> polygamma(int n, std::complex<float> z) {
    return static_cast<std::complex<float>>(polygamma(n, static_cast<std::complex<double>>(z)));
}

} // namespace xsf
