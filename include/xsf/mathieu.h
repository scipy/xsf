/* Kernels for computing Mathieu functions. This file does not contain complete
 * implementations of the Mathieu functions, but only kernels which can be
 * used within SciPy and CuPy to implement the Mathieu functions. One currently
 * must supply an eigenvalue and eigenvector solver for symmetric tridiagonal
 * matrices, such as dstevd from LAPACK, to construct a complete implementation.
 *
 * This code is based on code written by Stuart Brorson (github @brorson) which is
 * available at https://github.com/scipy/xsf/pull/99. It isolates out the parts
 * which are fit for xsf as a library of simple numerical kernels which can be used
 * on both CPU and GPU. Parts for which the approach used on CPU and GPU must
 * necessarily be different will be left to SciPy and CuPy respectively. Hence
 * the lack of the eigenvalue/eigenvector finding step here. The code has
 * has been refactored to reduce code duplication, but is otherwise faithful
 * to the original (with the exception that a policy template arg has been
 * added to allow either degrees or radians to be used for angular arguments).
 *
 * Stuart Brorson's implementations follow prototypes written in MATLAB and
 * maintained on GitHub at https://github.com/brorson/MathieuFcnsFourier.
 * A full write up for the algorithms is available at
 * https://github.com/brorson/ScipyMathieuPaper.
 *
 * SciPy Developers 2026
 */

#pragma once

#include "xsf/config.h"
#include "xsf/error.h"
#include "xsf/trig.h"

namespace xsf {
namespace mathieu {

    /* Used in template arguments for parity. Can stand for either
     * function parity (even Mathieu functions ce or odd Mathieu functions se)
     * or the parity of the order ``m``. */
    enum class Parity { Even, Odd };

    /* Policy to determine whether to use radians or degrees for angular
     * Mathieu functions. */
    enum class AngleUnitPolicy { Radians, Degrees };

    // Get index of characteristic value in sorted array of eigenvalues.
    template <Parity FuncParity, typename T>
    XSF_HOST_DEVICE T cv_index(T m) {
        static_assert(std::is_integral_v<T>, "m must be of integer type");

        if constexpr (FuncParity == Parity::Even) {
            return m / 2;
        }
        return (m - 1) / 2;
    }

    /* Calculates the first element along the diagonal of the recurrence
     * relation matrix as a function of the shape parameter ``q``.
     *
     * Templated on function parity (even vs odd Mathieu functions) and
     * order parity (whether the order parameter ``m`` is even or odd).
     */
    template <Parity FuncParity, Parity OrderParity>
    XSF_HOST_DEVICE double d0(double q) {
        constexpr auto Even = Parity::Even;
        constexpr auto Odd = Parity::Odd;

        if constexpr (FuncParity == Even && OrderParity == Even) {
            return 0.0;
        }
        if constexpr (FuncParity == Even && OrderParity == Odd) {
            return 1.0 + q;
        }
        if constexpr (FuncParity == Odd && OrderParity == Even) {
            return 4.0;
        }
        // FuncParity == Odd && OrderParity == Odd
        return 1.0 - q;
    }

    /* Calculates the first element along the off-diagonal of the recurrence
     * relation matrix as a function of the shape parameter ``q``.
     *
     * Templated on function parity (even vs odd Mathieu functions) and
     * order parity (whether the order parameter ``m`` is even or odd).
     */
    template <Parity FuncParity, Parity OrderParity>
    XSF_HOST_DEVICE double e0(double q) {
        constexpr auto Even = Parity::Even;

        if constexpr (FuncParity == Even && OrderParity == Even) {
            return M_SQRT2 * q;
        }
        return q;
    }

    /* Calculates the square root of the diagonal entries of the recurrence
     * relation matrix as a function of the index into the diagonal. */
    template <Parity FuncParity, Parity OrderParity, typename U>
    XSF_HOST_DEVICE double sqrt_di(U j) {
        static_assert(std::is_integral_v<U>, "j must be of integer type");
        constexpr auto Even = Parity::Even;
        constexpr auto Odd = Parity::Odd;
        if constexpr (OrderParity == Odd) {
            return 2.0 * static_cast<double>(j) + 1.0;
        }
        if constexpr (FuncParity == Even) {
            return 2.0 * static_cast<double>(j);
        }
        // FuncParity == Odd && OrderParity == Even
        return 2.0 * (static_cast<double>(j) + 1.0);
    }

    /*-----------------------------------------------
      This creates the recurrence relation matrix for
      the even-even (ce_2n), even-odd (ce_2n+1),
      (se_2n) -- sometimes called se_2n+2, and odd-odd (se_2n+1)
      Mathieu functions, depending on the values of the template
      arguments FuncParity and OrderParity.

      Template arguments:

      FuncParity = parity of Mathieu function, Parity::Even for the
      even Mathieu function and Parity::Odd for the odd
      Mathieu function.

      OrderParity = parity of the order m of a Mathieu function.

      Inputs:
      q = shape parameter.

      Outputs:
      D = Diagonal elements in recurrence matrix. mdspan view of a 1d array of
      size n.
      E = Off-diagonal elements in recurrence matrix. mdspan view of 1d array of
      size n - 1.
      -------------------------------------------------*/

    template <Parity FuncParity, Parity OrderParity, typename OutputMat>
    XSF_HOST_DEVICE void make_matrix(typename OutputMat::value_type q, OutputMat D, OutputMat E) {
        auto n = D.extent(0);
        if (n == 0) {
            return;
        }

        D(0) = d0<FuncParity, OrderParity>(q);

        for (decltype(n) j = 1; j < n; j++) {
            auto tj = sqrt_di<FuncParity, OrderParity>(j);
            D(j) = tj * tj;
        }

        if (n == 1) {
            return;
        }

        // Make off-diagonal entries
        E(0) = e0<FuncParity, OrderParity>(q);
        for (decltype(n) j = 1; j < (n - 1); j++) {
            E(j) = q;
        }
    }

    // Determine size of recursion matrix to use.
    XSF_HOST_DEVICE inline int get_partial_sum_N(int m, double q) {
        int N;

        // This is sort of ad-hoc ...
        double abs_q = std::abs(q);
        if (abs_q > 1.0) {
            double qq = std::log10(abs_q); // I need to use size of q to compute N.
            N = m + 25 + 10 * qq;
        } else {
            N = m + 25;
        }
        return N;
    }

    /* Sum the Fourier series for computing angular Mathieu functions.
     *
     * The Fourier coefficients ``X`` are the values of the eigenvector
     * associated to the characteristic eigenvalue. xsf does not currently
     * supply a kernel for computing this eigenvector.
     *
     * Template arguments:
     *
     * FuncParity = parity of Mathieu function, Parity::Even for the
     * even Mathieu function and Parity::Odd for the odd
     * Mathieu function.
     *
     * OrderParity = parity of the order m of a Mathieu function.
     *
     * AngleUnits = AngleUnitPolicy::Radians or AngleUnitPolicy::Degrees.
     *
     * Inputs:
     * X = A 1d mdspan view of the Fourier coefficients.
     * v = Angular argument, in either radians or degrees depending on the
     * value of the AngleUnits template argument.
     *
     * Ouputs:
     * out = Value of angular Mathieu function for given Fourier coefficients at angle v.
     * out_diff = Derivative of angular Mathieu function for given Fourier coefficients at angle v.
     */
    template <
        Parity FuncParity, Parity OrderParity, AngleUnitPolicy AngleUnits = AngleUnitPolicy::Radians, typename InputMat>
    XSF_HOST_DEVICE void sum_fourier_series(
        InputMat X, double v, typename InputMat::value_type &out, typename InputMat::value_type &out_diff
    ) {
        auto N = X.extent(0);
        // Local scope variables used in summing the Fourier series.
        double tt, td, xep{0.0}, xem{0.0}, xedp{0.0}, xedm{0.0};

        // Sum from smallest to largest coeff.
        for (decltype(N) kp1 = N; kp1 > 0; kp1--) {
            auto k = kp1 - 1;
            double r = sqrt_di<FuncParity, OrderParity>(k);
            double phi = r * v;
            double x_cos, x_sin;
            if constexpr (AngleUnits == AngleUnitPolicy::Degrees) {
                x_cos = xsf::cospi(phi / 180.0);
                x_sin = xsf::sinpi(phi / 180.0);
            } else {
                x_cos = std::cos(phi);
                x_sin = std::sin(phi);
            }
            if constexpr (FuncParity == Parity::Even) {
                tt = X(k) * x_cos;
                td = -r * X(k) * x_sin;
            } else {
                tt = X(k) * x_sin;
                td = r * X(k) * x_cos;
            }
            if (tt < 0) {
                xem = xem + tt; // Neg running sum
            } else {
                xep = xep + tt; // Pos running sum
            }

            if (td < 0) {
                xedm = xedm + td;
            } else {
                xedp = xedp + td;
            }
        } // for (k=(N-1) ...

        // This makes sure the fcn has the right overall sign for q<0.
        // Someday combine this with the above sums into the same for loop.
        double x = 0.0;
        for (decltype(N) l = 0; l < N; l++) {
            x += X(l);
        }
        auto sgn = std::copysign(1.0, x);
        out = sgn * (xep + xem);
        out_diff = sgn * (xedp + xedm);
    }

} // namespace mathieu
} // namespace xsf
