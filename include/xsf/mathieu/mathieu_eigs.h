#ifndef MATHIEU_EIGS_H
#define MATHIEU_EIGS_H

#include "../config.h"
#include "../error.h"
#include "make_matrix_tridiagonal.h"
#include "matrix_utils.h"
#include <vector>

/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the functions
 * which return the Mathieu eigenvalues (characteristic values) a and
 * b as a function of parameter q.
 *
 * Stuart Brorson, summer 2025.
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
    // This is the Mathieu characteristic value (eigenvalue)
    // a for even fcns.
    int mathieu_a(int m, double q, double *a) {
        // printf("--> mathieu_a, m = %d, q = %e\n", m, q);

        int N = m + 25; // Sets size of recursion matrix.  I make it larger than the order
                        // for accuracy's sake.
        int retcode = SF_ERROR_OK;

        if (m > 500) {
            // Don't support absurdly larger orders for now.
            *a = std::numeric_limits<double>::quiet_NaN();
            return SF_ERROR_DOMAIN;
        }

        // Allocate diag part of recursion matrix
        std::vector<double> D(N);

        // Allocate off-diag part of recursion matrix.
        std::vector<double> E(N - 1);

        // Allocate vector for eigenvectors
        std::vector<double> Z(N);

        // Do EVD
        if (m % 2 == 0) {
            // Even order m
            retcode = make_matrix_ee(N, q, D.data(), E.data());
            if (retcode != SF_ERROR_OK) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER; // Not sure what went wrong.
            }

            char jobz = 'N'; // Eigenvalues only.

            /* Allocate the optimal workspace */
            int lwork = 1 + 2 * N;
            int liwork = 1;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            int info;

            /* Solve eigenproblem */
            dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

            // Check if dstevd was successful
            if (info != 0) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = m / 2;
            *a = D[idx];

        } else {
            // Odd order m
            retcode = make_matrix_eo(N, q, D.data(), E.data());
            if (retcode != SF_ERROR_OK) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER;
            }

            char jobz = 'N'; // Eigenvalues only.

            /* Allocate the optimal workspace */
            int lwork = 1 + 2 * N;
            int liwork = 1;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            int info;

            /* Solve eigenproblem */
            dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

            // Check if dstevd was successful
            if (info != 0) {
                *a = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = (m - 1) / 2;
            *a = D[idx];
        }

        // printf("<-- mathieu_a\n");
        return retcode;
    }

    //---------------------------------------------------------
    int mathieu_b(int m, double q, double *b) {
        // This computes the Mathieu characteristic value (eigenvalue)
        // for odd fcns.
        // printf("--> mathieu_b, m = %d, q = %e\n", m, q);
        int N = m + 25; // Sets size of recursion matrix
        int retcode = SF_ERROR_OK;

        if (m > 500) {
            // Don't support absurdly larger orders for now.
            *b = std::numeric_limits<double>::quiet_NaN();
            return SF_ERROR_DOMAIN;
        }

        // Allocate diag part of recursion matrix
        std::vector<double> D(N);

        // Allocate off-diag part of recursion matrix.
        std::vector<double> E(N - 1);

        // Allocate vector for eigenvectors
        std::vector<double> Z(N);

        // Do EVD
        if (m % 2 == 0) {
            // Even order m
            retcode = make_matrix_oe(N, q, D.data(), E.data());
            if (retcode != SF_ERROR_OK) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER;
            }

            char jobz = 'N'; // Eigenvalues only.

            /* Allocate the optimal workspace */
            int lwork = 1 + 2 * N;
            int liwork = 1;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            int info;

            /* Solve eigenproblem */
            dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

            if (retcode != 0) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = (m - 2) / 2;
            *b = D[idx];

        } else {
            // Odd order m
            retcode = make_matrix_oo(N, q, D.data(), E.data());
            if (retcode != SF_ERROR_OK) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_OTHER;
            }

            char jobz = 'N'; // Eigenvalues only.

            /* Allocate the optimal workspace */
            int lwork = 1 + 2 * N;
            int liwork = 1;

            /* Allocate worksapce */
            std::vector<double> work(lwork);
            std::vector<int> iwork(liwork);

            int info;

            /* Solve eigenproblem */
            dstevd_(&jobz, &N, D.data(), E.data(), Z.data(), &N, work.data(), &lwork, iwork.data(), &liwork, &info);

            if (retcode != 0) {
                *b = std::numeric_limits<double>::quiet_NaN();
                return SF_ERROR_NO_RESULT;
            }

            // Sort ww vector from lowest to highest
            // quickSort(ww, 0, N-1);
            // print_matrix(ww, N, 1);

            // Now figure out which one to return.
            int idx = (m - 1) / 2;
            *b = D[idx];
        }

        // printf("<-- mathieu_b\n");
        return retcode;
    }

} // namespace mathieu
} // namespace xsf

#endif // #ifndef MATHIEU_EIGS_H
