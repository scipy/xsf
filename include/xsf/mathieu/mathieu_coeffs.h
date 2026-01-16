#ifndef MATHIEU_COEFFS_H
#define MATHIEU_COEFFS_H

#include "../config.h"
#include "../error.h"
#include "make_matrix_tridiagonal.h"
#include "matrix_utils.h"
#include <vector>

#define SQRT2 1.414213562373095

/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the functions
 * which return the Mathieu A and B coefficients used in the Fourier
 * series computing the Mathieu fcns.
 *
 * Stuart Brorson, Summer 2025.
 *
 */

/* DSTEVD_ prototype */
#ifdef __cplusplus
extern "C" {
#endif
// This Lapack routine computes the eigenvalues/vectors of a symmetric tridiagonal matrix.
// I need to do a forward declaration here.
void dstevd_(
    const char *jobz,  // 'N' = eigenvalues only, 'V' = eigenvalues + eigenvectors
    const int *n,      // Matrix dimension
    double *d,         // Diagonal elements (input), eigenvalues (output)
    double *e,         // Off-diagonal elements (input), destroyed on output
    double *z,         // Eigenvectors (output if jobz='V')
    const int *ldz,    // Leading dimension of z (>= n if jobz='V', >= 1 if 'N')
    double *work,      // Workspace array
    const int *lwork,  // Size of work array
    int *iwork,        // Integer workspace array
    const int *liwork, // Size of iwork array
    int *info          // Status: 0 = success, <0 = illegal argument, >0 = convergence failure
);
#ifdef __cplusplus
}
#endif

namespace xsf {
namespace mathieu {

    //------------------------------------------------------
    int mathieu_coeffs_ee(int N, double q, int m, double *AA) {
        // Returns Fourier coeffs for the mth order ce_2n Mathieu fcn.
        // Allowed value of m = 0, 2, 4, 6, ...
        // Inputs:
        // N = size of recursion matrix to use.
        // q = frequency parameter
        // m = order of Mathieu fcn desired.
        // Output:
        // AA = length N vector preallocated to hold coeffs.
        // Returns SF_ERROR_OK if all goes well.

        int retcode = SF_ERROR_OK;

        // Bail out if m is not even.
        if (m % 2 != 0)
            return SF_ERROR_ARG;

        // Allocate diag part of recursion matrix
        std::vector<double> D(N);

        // Allocate off-diag part of recursion matrix.
        std::vector<double> E(N - 1);

        // Allocate return matrix for eigenvectors.
        std::vector<double> Z(N * N);

        // Do EVD
        retcode = make_matrix_ee(N, q, D.data(), E.data());
        if (retcode != 0) {
            return SF_ERROR_NO_RESULT;
        }

        char jobz = 'V'; // Eigenvalues & eigenvecs

        /* Allocate the optimal workspace */
        int lwork = 1 + 4 * N + N * N;
        int liwork = 3 + 5 * N;

        /* Allocate worksapce */
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);

        int info;

        /* Solve eigenproblem */
        dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

        // Check return code from dstevd and bail if it's not 0.
        if (info != 0) {
            return SF_ERROR_NO_RESULT;
        }

        // Sort D vector from lowest to highest
        // quickSort(D, 0, N-1);
        // print_matrix(D, N, 1);

        // Undo sqrt(2) in make_matrix by normalizing elet in first col by sqrt(2).
        int idx;
        int col = m / 2;
        idx = MATRIX_IDX(N, 0, col);
        AA[0] = Z[idx] / SQRT2;
        // Transfer remaining elets in correct col to coeff vector.
        // Lapack is column major
        for (int j = 1; j < N; j++) {
            idx = MATRIX_IDX(N, j, col);
            AA[j] = Z[idx];
        }
        return retcode;
    }

    //------------------------------------------------------
    int mathieu_coeffs_eo(int N, double q, int m, double *AA) {
        // Returns Fourier coeffs for the mth order ce_2n+1 Mathieu fcn.
        // Allowed value of m = 1, 3, 5, 7 ...

        int retcode = SF_ERROR_OK;

        // Bail out if m is not odd.
        if (m % 2 != 1)
            return SF_ERROR_ARG;

        // Allocate diag part of recursion matrix
        std::vector<double> D(N);

        // Allocate off-diag part of recursion matrix.
        std::vector<double> E(N - 1);

        // Allocate matrix for eigenvectors
        std::vector<double> Z(N * N);

        // Do EVD
        retcode = make_matrix_eo(N, q, D.data(), E.data());
        if (retcode != 0) {
            return SF_ERROR_NO_RESULT;
        }

        char jobz = 'V'; // Eigenvalues & eigenvecs

        /* Allocate the optimal workspace */
        int lwork = 1 + 4 * N + N * N;
        int liwork = 3 + 5 * N;

        /* Allocate worksapce */
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);

        int info;

        /* Solve eigenproblem */
        dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

        // Check return code from dstevd and bail if it's not 0.
        if (info != 0) {
            return SF_ERROR_NO_RESULT;
        }

        // Sort D vector from lowest to highest
        // quickSort(D, 0, N-1);
        // print_matrix(D, N, 1);

        // Transfer correct col to coeff vector.
        int idx;
        int col = (m - 1) / 2;
        // Transfer elets in correct col to coeff vector.
        for (int j = 0; j < N; j++) {
            idx = MATRIX_IDX(N, j, col);
            AA[j] = Z[idx];
        }
        return retcode;
    }

    //------------------------------------------------------
    int mathieu_coeffs_oe(int N, double q, int m, double *AA) {
        // Returns Fourier coeffs for the mth order se_2n+2 Mathieu fcn.
        // Allowed value of m = 0, 2, 4, 6, ...
        // Inputs:
        // N = size of recursion matrix to use.
        // q = frequency parameter
        // m = order of Mathieu fcn desired.
        // Output:
        // AA = length N vector preallocated to hold coeffs.
        // Returns 0 if all goes well.  Must put check on calloc
        // here.

        int retcode = SF_ERROR_OK;

        // Bail out if m is not even or >= 2.
        // if ((m % 2 != 0) || (m < 2))
        //    return SF_ERROR_ARG;

        // Bail out if m is not even.
        if (m % 2 != 0)
            return SF_ERROR_ARG;

        // Allocate diag part of recursion matrix
        std::vector<double> D(N);

        // Allocate off-diag part of recursion matrix.
        std::vector<double> E(N - 1);

        // Allocate matrix for eigenvectors
        std::vector<double> Z(N * N);

        // Do EVD
        retcode = make_matrix_oe(N, q, D.data(), E.data());
        if (retcode != 0) {
            return SF_ERROR_NO_RESULT;
        }

        char jobz = 'V'; // Eigenvalues & eigenvecs

        /* Allocate the optimal workspace */
        int lwork = 1 + 4 * N + N * N;
        int liwork = 3 + 5 * N;

        /* Allocate worksapce */
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);

        int info;

        /* Solve eigenproblem */
        dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

        // Bail out if dstevd doesn't return 0.
        if (info != 0) {
            return SF_ERROR_NO_RESULT;
        }

        // Sort D vector from lowest to highest
        // quickSort(D, 0, N-1);
        // print_matrix(D, N, 1);

        // Transfer remaining elets in correct col to coeff vector.
        int idx;
        int col = m / 2;
        for (int j = 0; j < N; j++) {
            idx = MATRIX_IDX(N, j, col);
            AA[j] = Z[idx];
        }
        return retcode;
    }

    //------------------------------------------------------
    int mathieu_coeffs_oo(int N, double q, int m, double *AA) {
        // Returns Fourier coeffs for the mth order se_2n+1 Mathieu fcn.
        // Allowed value of m = 1, 3, 5, 7 ...

        int retcode = SF_ERROR_OK;

        // Bail out if m is not odd.
        if (m % 2 != 1)
            return SF_ERROR_ARG;

        // Allocate diag part of recursion matrix
        std::vector<double> D(N);

        // Allocate off-diag part of recursion matrix.
        std::vector<double> E(N - 1);

        // Allocate matrix for eigenvectors
        std::vector<double> Z(N * N);

        // Do EVD
        retcode = make_matrix_oo(N, q, D.data(), E.data());
        if (retcode != 0) {
            return SF_ERROR_NO_RESULT;
        }

        char jobz = 'V'; // Eigenvalues & eigenvecs

        /* Allocate the optimal workspace */
        int lwork = 1 + 4 * N + N * N;
        int liwork = 3 + 5 * N;

        /* Allocate worksapce */
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);

        int info;

        /* Solve eigenproblem */
        dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

        // Bail out if dstevd didn't return 0;
        if (info != 0) {
            return SF_ERROR_NO_RESULT;
        }

        // Sort D vector from lowest to highest
        // quickSort(D, 0, N-1);
        // print_matrix(D, N, 1);

        // Transfer correct col to coeff vector.
        int idx;
        int col = (m - 1) / 2;
        // Transfer elets in correct col to coeff vector.
        for (int j = 0; j < N; j++) {
            idx = MATRIX_IDX(N, j, col);
            AA[j] = Z[idx];
        }
        return retcode;
    }

} // namespace mathieu
} // namespace xsf

#endif // #ifndef MATHIEU_COEFFS_H
