/* 
 * Author: Lorenzo Peri
 */

#pragma once

#include "cephes/psi.h"
#include "cephes/zeta.h"
#include "config.h"
#include "error.h"
#include "trig.h"
#include "gamma.h"
#include "digamma.h"
#include "trigamma.h"

namespace xsf {
namespace detail {

    XSF_HOST_DEVICE inline std::complex<double>
    polygamma_forward_recurrence(int n, std::complex<double>const  z, std::complex<double> const psiz, int m){
        /* Compute polygamma(n, z + m) using polygamma(n, z) using the recurrence relation
        *    polygamma(n, z + 1) = polygamma(n, z) + (-1)^n * n! / z^(n+1).
        * See https://dlmf.nist.gov/5.15#E5 */

        std::complex<double> const np1_neg = ( -1. * static_cast<std::complex<double>>(n + 1));
        std::complex<double> zk, res = psiz;

        for (int k=0;k<m;k++){
            zk = (z + static_cast<double>(k));
            res += std::pow(zk, np1_neg);
        }

        return res;
    }

    XSF_HOST_DEVICE inline std::complex<double>
    polygamma_backward_recurrence(int n, std::complex<double>const  z, std::complex<double> const psiz, int m){
        /* Compute polygamma(n, z + m) using polygamma(n, z) using the recurrence relation
        *    polygamma(n, z + 1) = polygamma(n, z) + (-1)^n * n! / z^(n+1).
        * See https://dlmf.nist.gov/5.15#E5 */

        std::complex<double> const np1_neg = ( -1. * static_cast<std::complex<double>>(n + 1));
        std::complex<double> zk, res = psiz;

        for (int k=1;k<m+1;k++) {
            zk = (z - static_cast<double>(k));
            res -= std::pow(zk, np1_neg);
        }

        return res;
    }

    XSF_HOST_DEVICE inline std::complex<double> 
    polygamma_asymptotic_series(int n, std::complex<double> z) {
        /* Evaluate digamma using an asymptotic series. See
        * http://dlmf.nist.gov/5.15.E8
        * Higher order converge slower, so we need slightly more Bernoulli 
        * numbers than the digamma (Bernoulli numbers calculated via mpmath).
        */
        static constexpr std::array<double, 36> bernoulli2k = {
            0.16666666666666665741e+0, -0.033333333333333332871e+0, 0.023809523809523808202e+0, 
            -0.033333333333333332871e+0, 0.075757575757575759678e+0, -0.25311355311355310249e+0, 
            1.1666666666666667407e+0, -7.0921568627450977118e+0, 54.971177944862155584e+0, 
            -529.12424242424242493e+0, 6192.1231884057970092e+0, -86580.253113553117146e+0, 
            1425517.1666666667443e+0, -27298231.067816093564e+0, 601580873.90064239502e+0, 
            -15116315767.092157364e+0, 429614643061.16668701e+0, -13711655205088.332031e+0, 
            488332318973593.1875e+0, -19296579341940068.0e+0, 841693047573682560.0e+0, 
            -40338071854059454464.0e+0, 2.1150748638081992622e+21, -1.2086626522296526202e+23, 
            7.500866746076964166e+24, -5.0387781014810688499e+26, 3.6528776484818122276e+28, 
            -2.8498769302450882361e+30, 2.3865427499683627448e+32, -2.1399949257225334859e+34, 
            2.050097572347809739e+36, -2.0938005911346379301e+38, 2.2752696488463514863e+40, 
            -2.6257710286239577207e+42, 3.2125082102718031743e+44, -4.1598278166794711978e+46, 
        };

        std::complex<double> res, term, fac_coeff;
        std::complex<double> const zn_inv = std::pow(z, static_cast<double>(-n));
        std::complex<double> const rzz = 1. / z / z;
        std::complex<double> zfac = zn_inv;

        res = zn_inv * (1. / n + 0.5 / z);

        fac_coeff = (0.5 * static_cast<std::complex<double>>(n+1));

        for (int k=1;k<bernoulli2k.size()+1;k++) {

            zfac *= rzz;
            term = bernoulli2k[k-1] * fac_coeff * zfac;
            res += term;

            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(res)) {
                break;
            }

            double const num_coeff = (2.*k + n)*(2.*k + n + 1.);
            double const den_coeff = (2.*k + 2) * (2.*k + 1.);
            fac_coeff *= num_coeff / den_coeff;
        }

        return -res;
    }

    XSF_HOST_DEVICE inline std::complex<double> 
    polygamma_reflection_rhs_nleq25(int n, std::complex<double> z) {
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
         * However, this requires us to loop through the calculation n times.
         * For larger n, there is little point using this strategy. 
         * Recurrence should be used instead.
         * Moreover, allowing a maximum n allows us to statically allocate
         * the coefficient buffers and avoid any dynamical allocation.
         */
        
        if(n>25){
            // We are not going to use it for n>25.
            // Statically allocate the memory and check against misuse
            return std::numeric_limits<double>::quiet_NaN();
        }

        std::array<double, 27> coefficient_array{}; // all zeros at the start
        std::array<double, 27> next_coeffs; 

        coefficient_array[1] = 1.; // P_0(y) = y

        /* The coefficients blow up very rapidly roughly O(n!). So we compute the 
         * normalized RHS (-1)**(n-1) * pi * d^n/dz^n cot(pi z) / n!, and let the 
         * general prefactor in the main polygamma driver handle the n! scaling
         */
        for(int k=0; k<n; k++){
            std::fill(next_coeffs.begin(), next_coeffs.end(), 0.);
            for(int m=0; m<k+3; m++){
                double const term1 = (m + 1 <= k + 1) ? 
                    (-static_cast<double>(m + 1) * coefficient_array[m + 1]) :
                    0.;
                double const term2 = (m - 1 >= 0) ? 
                    (-static_cast<double>(m - 1) * coefficient_array[m - 1]) :
                    0.;
                next_coeffs[m] = (term1 + term2)/static_cast<double>(k + 1);
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
    if (n == 0){
        return digamma(z);
    } else if (n == 1){
        return trigamma(z);
    }

    /* (-1)^(n+1) */
    const double sign = n & 1 ? 1. : -1.;
    const double np1 = static_cast<double>(n + 1);

    double const factorial = gamma(np1);

    const double Hur_zeta = cephes::zeta(np1, z);

    return sign * factorial * Hur_zeta;
}

XSF_HOST_DEVICE inline float polygamma(int n, float z) { return static_cast<float>(polygamma(n, static_cast<double>(z))); }

XSF_HOST_DEVICE inline std::complex<double> polygamma(int n, std::complex<double> z) {
    /*
     * Compute the polygamma function for complex arguments. The strategy is:
     *
     * - If close to the origin, use a recurrence relation to step away
     * from the origin.
     * - If close to the negative real axis, use either the reflection 
     * or the recurrence formula to move to the right halfplane.
     * - If |z| is large (> 16), use the asymptotic series.
     * - If |z| is small, use a recurrence relation to make |z| large
     * enough to use the asymptotic series.
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
    } else if (n==1){
        return trigamma(z);
    }

    if(std::abs(z.imag()) < std::numeric_limits<double>::epsilon()){
        // use the zeta, because it evaluates to better precision
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
    double const smallabsz = 22.5;
    /* When recursively stepping away via recurrence, it is best to leave numbers away from the
     * real axis alone. So we want smallimag < smallabsz. 
     * Some testing yielded  this value of smallimag as a sensible value.
     * To increase the precision, having this value depend on n with some ansatz would probably be ideal.
     */
    double const smallimag = 7.5;

    if (z.real() <= 0.0 && std::ceil(z.real()) == z) {
        // Poles
        set_error("polygamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    if (n <= 25){ 
        // For small orders, we can use the reflection formula
        if(z.real() < 0 && std::abs(z.imag()) < smallabsz) {
            std::complex<double> const reflection_rhs = detail::polygamma_reflection_rhs_nleq25(n, z);
            res -= reflection_rhs;
            sign *= (n & 1) ? -1. : 1.;
            z = 1. - z;
            absz = std::abs(z);
        } 
        if(absz < 0.5){
            std::complex<double> const np1_neg =  -1. * static_cast<std::complex<double>>(n + 1);
            res -= sign * std::pow(z, np1_neg);
            z += 1;
            absz = std::abs(z);
        }
    } else {
        /* For large n, the reflection formula blows up an leads to loss of precision.
         * We can however use the recurrence formula to step away from poles
         * in the negative real plane */
        while((absz < 0.5) || (z.real() < 0 && std::abs(z.imag()) < smallimag)){
            std::complex<double> const np1_neg =  -1. * static_cast<std::complex<double>>(n + 1);
            res -= std::pow(z, np1_neg);
            z += 1;
            absz = std::abs(z);
        }
    }

    if (absz > smallabsz) {
        res += sign * detail::polygamma_asymptotic_series(n, z);
    } else if (z.real() >= 0.0) {
        double m = std::trunc(smallabsz - absz) + 1;
        std::complex<double> init = detail::polygamma_asymptotic_series(n, z + m);
        res += sign * detail::polygamma_backward_recurrence(n, z + m, init, m);
    } else {
        double m = std::trunc(smallabsz - absz) - 1;
        std::complex<double> init = detail::polygamma_asymptotic_series(n, z - m);
        res += sign * detail::polygamma_forward_recurrence(n, z - m, init, m);
    }

    return prefactor * res;
}

XSF_HOST_DEVICE inline std::complex<float> polygamma(int n, std::complex<float> z) {
    return static_cast<std::complex<float>>(polygamma(n, static_cast<std::complex<double>>(z)));
}

} // namespace xsf
