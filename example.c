#include <stdio.h>
#include <stdlib.h>
#include "matrix_io.h"

void print_csr_matrix(Matrix mat) {
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
    for (int i=0; i<mat.nz; i++)
        fprintf(stdout, "(%d, %d, %f)\n", mat.I[i], mat.J[i], mat.val[i]);
}


/** C = AA
    all matrices are in csr format */
void self_multiplication(int M, int N, int *I, int *J, double *val) {
}

/** C = AB
    all matrices are in coo format
    use Expansion, Sort, Compress(ESC) algorithm */
void esc_multiplication(Matrix A, Matrix B) {
}


int main(int argc, char *argv[]) {
    Matrix A;   // Matrix is defined in matrix_io.h

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    
    read_mm_matrix_coo(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_coo_matrix(A);

    read_mm_matrix_csr(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_csr_matrix(A);

    read_mm_matrix_csc(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_csc_matrix(A);
    return 0;
}
