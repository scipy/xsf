#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include "matrix_utils.h"
#include <math.h>
#include <stdio.h>

// These fcns are meant to make it easier to deal with
// matrices in C.  We use col major format since that's
// what underlies Lapack.

// returns +/-1 depending upon sign of x
#define SIGN(x) (((x) > 0) - ((x) < 0))

// Macros to extract matrix index and element.
// Matrix is NxN, i = row idx, j = col idx.
// MATRIX_IDX is where col major format is enforced.
#define MATRIX_IDX(N, I, J) (((N) * (I)) + (J))
#define MATRIX_ELEMENT(A, m, n, i, j) A[MATRIX_IDX(n, i, j)]

// Min and max macros for scalars.
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

//===========================================================
// This file holds utility functions for dealing with vectors
// and matrices.  The idea is to be able to reuse common matrix
// operations.  I will name the utils analogously to their names
// in Matlab.
// Note that C matrices are row-major.

namespace xsf {
namespace mathieu {

    //-----------------------------------------------------
    void print_matrix(const double *A, int m, int n) {
        // prints matrix as 2-dimensional tablei -- this is how we
        // usually think of matrices.
        int i, j;
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                printf("% 10.4e  ", MATRIX_ELEMENT(A, m, n, i, j));
            }
            printf("\n");
        }
    }

    //-----------------------------------------------------
    // Stuff to sort a vector.
    // Function to swap two elements
    void swap(double *a, double *b) {
        double temp = *a;
        *a = *b;
        *b = temp;
    }

    // Partition function for quicksort
    int partition(double *arr, int low, int high) {
        double pivot = arr[high]; // Choose last element as pivot
        int i = (low - 1);        // Index of smaller element

        for (int j = low; j <= high - 1; j++) {
            // If current element is smaller than or equal to pivot
            if (arr[j] <= pivot) {
                i++;
                swap(&arr[i], &arr[j]);
            }
        }
        swap(&arr[i + 1], &arr[high]);
        return (i + 1);
    }

    // Quicksort function
    void quickSort(double *arr, int low, int high) {
        if (low < high) {
            // Partition the array and get pivot index
            int pivotIndex = partition(arr, low, high);

            // Recursively sort elements before and after partition
            quickSort(arr, low, pivotIndex - 1);
            quickSort(arr, pivotIndex + 1, high);
        }
    }

} // namespace mathieu
} // namespace xsf

#endif // #ifndef MATRIX_UTILS_H
