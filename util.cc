#include <stdio.h>
#include <stdlib.h>
#include "matrix_io.h"
#include "util.h"

void print_matrix_head(Matrix mat) {
    printf("ROWS: %d, COLS: %d, NNZ: %d\n", mat.M, mat.N, mat.nz);
}

void print_csr_matrix(Matrix mat) {
    print_matrix_head(mat);
    for(int i = 0; i < mat.M; i++) {
        // access ith column
        for(int j=mat.I[i]; j < mat.I[i+1]; j++) {
            double a = mat.val[j];
            printf("%f ", a);
        }
        printf("\n");
    }
}

void print_csc_matrix(Matrix mat) {
    print_matrix_head(mat);
    for(int i = 0; i < mat.N; i++) {
        // access ith row
        for(int j=mat.J[i]; j < mat.J[i+1]; j++) {
            double a = mat.val[j];
            printf("%f ", a);
        }
        printf("\n");
    }
}

void print_coo_matrix(Matrix mat) {
    print_matrix_head(mat);
    for (int i=0; i<mat.nz; i++)
        fprintf(stdout, "(%d, %d, %f)\n", mat.I[i], mat.J[i], mat.val[i]);
}