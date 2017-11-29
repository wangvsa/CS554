#include <stdio.h>
#include <stdlib.h>
#include "matrix_io.h"

void print_csr_matrix(Matrix mat) {
    printf("ROWS: %d, COLS: %d, NNZ: %d\n", mat.M, mat.N, mat.nz);
    for(int i = 0; i < mat.M; i++) {
        // access ith row
        for(int j=0; j < mat.I[i+1]-mat.I[i]; j++) {
            double a = mat.val[mat.I[i] + j];
            printf("%f ", a);
        }
        printf("\n");
    }
}

void print_csc_matrix(Matrix mat) {
    printf("ROWS: %d, COLS: %d, NNZ: %d\n", mat.M, mat.N, mat.nz);
    for(int i = 0; i < mat.N; i++) {
        // access ith column
        for(int j=0; j < mat.J[i+1]-mat.J[i]; j++) {
            double a = mat.val[mat.J[i] + j];
            printf("%f ", a);
        }
        printf("\n");
    }
}

void print_coo_matrix(Matrix mat) {
    printf("ROWS: %d, COLS: %d, NNZ: %d\n", mat.M, mat.N, mat.nz);
    for (int i=0; i<mat.nz; i++)
        fprintf(stdout, "(%d, %d, %f)\n", mat.I[i], mat.J[i], mat.val[i]);
}
