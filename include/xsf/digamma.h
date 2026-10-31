/* Translated from Cython into C++ by SciPy developers in 2024.
 * Original header comment appears below.
 */

/* An implementation of the digamma function for complex arguments.
 *
 * Author: Josh Wilson
 *
 * Distributed under the same license as Scipy.
 *
 * Sources:
 * [1] "The Digital Library of Mathematical Functions", dlmf.nist.gov
 *
 * [2] mpmath (version 0.19), http://mpmath.org
 */

#pragma once

#include "cephes/psi.h"
#include "cephes/zeta.h"
#include "config.h"
#include "error.h"
#include "trig.h"

namespace xsf {
namespace detail {
    // All of the following were computed with mpmath
    // Location of the positive root
    constexpr double digamma_posroot = 1.4616321449683623;
    // Value of the positive root
    constexpr double digamma_posrootval = -9.2412655217294275e-17;
    // Location of the negative root
    constexpr double digamma_negroot = -0.504083008264455409;
    // Value of the negative root
    constexpr double digamma_negrootval = 7.2897639029768949e-17;

    // zeta(n + 1, digamma_posroot) for n = 1, ..., 99, computed with mpmath
    constexpr double digamma_posroot_zeta[] = {
        0.9676722454476212,     0.4427631689835922,     0.2584997609556511,     0.16394270544240658,
        0.10782405069126241,    0.07219956125645474,    0.04880428816414313,    0.033161126474847376,
        0.022597648232218118,   0.01542476590494897,    0.010538791616612184,   0.007204534386356875,
        0.004926781395729858,   0.0033698016554393312,  0.00230512632673493,    0.001576936771430199,
        0.0010788252019162978,  0.000738070938996006,   0.0005049532658346027,  0.0003454680251063082,
        0.00023635601564027086, 0.00016170622091974827, 0.00011063372768747428, 7.569179582195078e-05,
        5.17857579522209e-05,   3.543007094765967e-05,  2.4240066118601362e-05, 1.6584242271854165e-05,
        1.1346384584663872e-05, 7.76281766846211e-06,   5.3110609208898746e-06, 3.6336507898010536e-06,
        2.4860227331295433e-06, 1.7008538854332646e-06, 1.163667536354887e-06,  7.96142543124199e-07,
        5.446941930669458e-07,  3.726616128343832e-07,  2.549626552021561e-07,  1.7443695117727793e-07,
        1.1934394829830275e-07, 8.165115189488431e-08,  5.586299683532188e-08,  3.8219600619175055e-08,
        2.6148576951961945e-08, 1.7889984864911548e-08, 1.22397314032337e-08,   8.374016297671817e-09,
        5.7292228598500124e-09, 3.919743335942345e-09,  2.681757752408118e-09,  1.8347692760029745e-09,
        1.2552879890596146e-09, 8.58826205608992e-10,   5.875802667349371e-10,  4.020028354997662e-10,
        2.750369420128259e-10,  1.881711092354027e-10,  1.287404015320643e-10,  8.807989204072243e-11,
        6.026132658886509e-11,  4.1228791249777315e-11, 2.8207364891163943e-11, 1.9298538957471046e-11,
        1.3203417168886059e-11, 9.033337980653022e-12,  6.180308781351105e-12,  4.2283612895533845e-12,
        2.8929038705870207e-12, 1.979228412939453e-12,  1.3541221159871894e-12, 9.264452212882197e-13,
        6.338429436417964e-13,  4.3365421718712694e-13, 2.966917624793437e-13,  2.0298661568213233e-13,
        1.3887667726856546e-13, 9.501479407568141e-14,  6.500595543329273e-14,  4.447490817513446e-14,
        3.0428249904217275e-14, 2.0817994465273552e-14, 1.4242977986587842e-14, 9.7445708454203e-15,
        6.666910603304518e-15,  4.561278038564776e-15,  3.1206744147403155e-15, 2.135061428063943e-15,
        1.4607378712994558e-15, 9.993881677602775e-16,  6.837480765599266e-16,  4.677976458808154e-16,
        3.2005155845210365e-16, 2.1896860954644052e-16, 1.4981102481923055e-16, 1.0249571024758304e-16,
        7.012414895254074e-17,  4.7976605600760517e-17, 3.28239945775132e-17
    };

    // zeta(n + 1, digamma_negroot) for n = 1, ..., 99, computed with mpmath
    constexpr double digamma_negroot_zeta[] = {
        8.939798558792134,      0.8093454625306499,     32.258572064915505,     2.7610116487908227,
        128.27342343090902,     14.706721735530332,     513.2704629562861,      75.37700066344615,
        2055.5357119463315,     368.5783503832184,      8234.667099203069,      1743.3646088673593,
        32997.785655255044,     8051.818155363481,      132262.95480834588,     36529.945969485576,
        530280.3912420176,      163451.7934270022,      2126602.372596875,      723318.8861804291,
        8530604.878543982,      3172131.31557447,       34228316.69580112,      13807481.99726386,
        137373435.87411496,     59721085.41129418,      551481251.8585448,      256915093.7987531,
        2214464568.0139103,     1100067148.9416316,     8894384049.698204,      4691108165.464198,
        35733149689.1997,       19932951831.31558,      143593265361.61185,     84427777205.99568,
        577169303800.4994,      356587583747.61273,     2320479554068.252,      1502246808129.2224,
        9331609972791.018,      6314216613755.217,      37535190870239.27,      26484684511361.957,
        151015931286965.84,     110878926240539.7,      607726007220707.0,      463397408554349.8,
        2446203594688412.5,     1933621531448282.8,     9848631064699276.0,     8056719714934248.0,
        3.966033418706299e+16,  3.352458218519522e+16,  1.5974703846095734e+17, 1.3932560733139752e+17,
        6.435820112041242e+17,  5.783630021481716e+17,  2.5933926695029007e+18, 2.398322687545823e+18,
        1.0452606668648067e+19, 9.935377914500428e+18,  4.2137766423157514e+19, 4.112064452831596e+19,
        1.6990556309721845e+20, 1.7004372382577053e+20, 6.8522226506873e+20,    7.026038929889101e+20,
        2.764023275935674e+21,  2.900893253646807e+21,  1.115159170357447e+22,  1.196862745715766e+22,
        4.5000311645821485e+22, 4.9347693676470485e+22, 1.816252543790399e+23,  2.0333695330919317e+23,
        7.331917212319315e+23,  8.373526739745053e+23,  2.9603156332049827e+24, 3.4463364563646236e+24,
        1.19546310752253e+25,   1.4176778047087605e+25, 4.828481533255975e+25,  5.828817746743598e+25,
        1.950561973063074e+26,  2.395400118392808e+26,  7.881016936454601e+26,  9.839682679529297e+26,
        3.1847595510935494e+27, 4.040169682322971e+27,  1.287186491919594e+28,  1.6582231942712865e+28,
        5.20325725283157e+28,   6.803305986480627e+28,  2.1036655452544412e+29, 2.7902227277743842e+29,
        8.506369437934226e+29,  1.1439506354873035e+30, 3.440143256278426e+30
    };

    // Limit of digamma(z) for non-finite complex z.
    // Precondition: !std::isfinite(z.real()) || !std::isfinite(z.imag()).
    template <typename T>
    XSF_HOST_DEVICE cxx::complex<T> digamma_nonfinite_limit(cxx::complex<T> z) {
        if (cxx::isnan(z.real()) || !cxx::isfinite(z.imag()) || z.real() <= T{0})
            return {cxx::numeric_limits<T>::quiet_NaN(), cxx::numeric_limits<T>::quiet_NaN()};
        return {cxx::numeric_limits<T>::infinity(), T{0}};
    }

    template <typename T>
    XSF_HOST_DEVICE T digamma_zeta_series(T z, double root, double rootval, const double *zeta) {
        T res = rootval;
        T coeff = -1.0;

        z = z - root;
        T term;
        for (int n = 1; n < 100; n++) {
            coeff *= -z;
            term = coeff * zeta[n - 1];
            res += term;
            if (cxx::abs(term) < cxx::numeric_limits<double>::epsilon() * cxx::abs(res)) {
                break;
            }
        }
        return res;
    }

    XSF_HOST_DEVICE inline cxx::complex<double>
    digamma_forward_recurrence(cxx::complex<double> z, cxx::complex<double> psiz, int n) {
        /* Compute digamma(z + n) using digamma(z) using the recurrence
         * relation
         *
         * digamma(z + 1) = digamma(z) + 1/z.
         *
         * See https://dlmf.nist.gov/5.5#E2 */
        cxx::complex<double> res = psiz;

        for (int k = 0; k < n; k++) {
            res += 1.0 / (z + static_cast<double>(k));
        }
        return res;
    }

    XSF_HOST_DEVICE inline cxx::complex<double>
    digamma_backward_recurrence(cxx::complex<double> z, cxx::complex<double> psiz, int n) {
        /* Compute digamma(z - n) using digamma(z) and a recurrence relation. */
        cxx::complex<double> res = psiz;

        for (int k = 1; k < n + 1; k++) {
            res -= 1.0 / (z - static_cast<double>(k));
        }
        return res;
    }

    XSF_HOST_DEVICE inline cxx::complex<double> digamma_asymptotic_series(cxx::complex<double> z) {
        /* Evaluate digamma using an asymptotic series. See
         *
         * https://dlmf.nist.gov/5.11#E2 */
        double bernoulli2k[] = {0.166666666666666667,   -0.0333333333333333333, 0.0238095238095238095,
                                -0.0333333333333333333, 0.0757575757575757576,  -0.253113553113553114,
                                1.16666666666666667,    -7.09215686274509804,   54.9711779448621554,
                                -529.124242424242424,   6192.12318840579710,    -86580.2531135531136,
                                1425517.16666666667,    -27298231.0678160920,   601580873.900642368,
                                -15116315767.0921569};
        cxx::complex<double> rzz = 1.0 / z / z;
        cxx::complex<double> zfac = 1.0;
        cxx::complex<double> term;
        cxx::complex<double> res;

        if (!(cxx::isfinite(z.real()) && cxx::isfinite(z.imag()))) {
            /* Check for infinity (or nan) and return early.
             * Result of division by complex infinity is implementation dependent.
             * and has been observed to vary between C++ stdlib and CUDA stdlib.
             */
            return cxx::log(z);
        }

        res = cxx::log(z) - 0.5 / z;

        for (int k = 1; k < 17; k++) {
            zfac *= rzz;
            term = -bernoulli2k[k - 1] * zfac / (2 * static_cast<double>(k));
            res += term;
            if (cxx::abs(term) < cxx::numeric_limits<double>::epsilon() * cxx::abs(res)) {
                break;
            }
        }
        return res;
    }

} // namespace detail

XSF_HOST_DEVICE inline double digamma(double z) {
    /* Wrap Cephes' psi to take advantage of the series expansion around
     * the smallest negative zero.
     */
    if (cxx::abs(z - detail::digamma_negroot) < 0.3) {
        return detail::digamma_zeta_series(
            z, detail::digamma_negroot, detail::digamma_negrootval, detail::digamma_negroot_zeta
        );
    }
    return cephes::psi(z);
}

XSF_HOST_DEVICE inline float digamma(float z) { return static_cast<float>(digamma(static_cast<double>(z))); }

XSF_HOST_DEVICE inline cxx::complex<double> digamma(cxx::complex<double> z) {
    /*
     * Compute the digamma function for complex arguments. The strategy
     * is:
     *
     * - Around the two zeros closest to the origin (posroot and negroot)
     * use a Taylor series with precomputed zero order coefficient.
     * - If close to the origin, use a recurrence relation to step away
     * from the origin.
     * - If close to the negative real axis, use the reflection formula
     * to move to the right halfplane.
     * - If |z| is large (> 16), use the asymptotic series.
     * - If |z| is small, use a recurrence relation to make |z| large
     * enough to use the asymptotic series.
     */

    // Non-finite z: std::abs(complex) can return NaN, false-routing dispatch.
    if (!cxx::isfinite(z.real()) || !cxx::isfinite(z.imag())) {
        return detail::digamma_nonfinite_limit(z);
    }

    double absz = cxx::abs(z);
    cxx::complex<double> res = 0;
    /* Use the asymptotic series for z away from the negative real axis
     * with abs(z) > smallabsz. */
    int smallabsz = 16;
    /* Use the reflection principle for z with z.real < 0 that are within
     * smallimag of the negative real axis.
     * int smallimag = 6  # unused below except in a comment */

    if (z.real() <= 0.0 && cxx::ceil(z.real()) == z) {
        // Poles
        set_error("digamma", SF_ERROR_SINGULAR, NULL);
        return {cxx::numeric_limits<double>::quiet_NaN(), cxx::numeric_limits<double>::quiet_NaN()};
    }
    if (cxx::abs(z - detail::digamma_negroot) < 0.3) {
        // First negative root.
        return detail::digamma_zeta_series(
            z, detail::digamma_negroot, detail::digamma_negrootval, detail::digamma_negroot_zeta
        );
    }

    if (z.real() < 0 && cxx::abs(z.imag()) < smallabsz) {
        /* Reflection formula for digamma. See
         *
         *https://dlmf.nist.gov/5.5#E4
         */
        res = -M_PI * cospi(z) / sinpi(z);
        z = 1.0 - z;
        absz = cxx::abs(z);
    }

    if (absz < 0.5) {
        /* Use one step of the recurrence relation to step away from
         * the pole. */
        res = -1.0 / z;
        z += 1.0;
        absz = cxx::abs(z);
    }

    if (cxx::abs(z - detail::digamma_posroot) < 0.5) {
        res += detail::digamma_zeta_series(
            z, detail::digamma_posroot, detail::digamma_posrootval, detail::digamma_posroot_zeta
        );
    } else if (absz > smallabsz) {
        res += detail::digamma_asymptotic_series(z);
    } else if (z.real() >= 0.0) {
        double n = cxx::trunc(smallabsz - absz) + 1;
        cxx::complex<double> init = detail::digamma_asymptotic_series(z + n);
        res += detail::digamma_backward_recurrence(z + n, init, n);
    } else {
        // z.real() < 0, absz < smallabsz, and z.imag() > smallimag
        double n = cxx::trunc(smallabsz - absz) - 1;
        cxx::complex<double> init = detail::digamma_asymptotic_series(z - n);
        res += detail::digamma_forward_recurrence(z - n, init, n);
    }
    return res;
}

XSF_HOST_DEVICE inline cxx::complex<float> digamma(cxx::complex<float> z) {
    return static_cast<cxx::complex<float>>(digamma(static_cast<cxx::complex<double>>(z)));
}

} // namespace xsf
