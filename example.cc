#include <stdio.h>
#include <stdlib.h>
#include "matrix_io.h"
#include "util.h"
#include "esc.h"

int main(int argc, char *argv[]) {
    Matrix A, B;   // Matrix is defined in matrix_io.h

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }

    /*
    read_mm_matrix_coo(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_coo_matrix(A);

    read_mm_matrix_coo(argv[2], &(B.M), &(B.N), &(B.nz), &(B.I), &(B.J), &(B.val));
    print_coo_matrix(B);

    esc_multiplication(A, B);
    */


    read_mm_matrix_csc(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_matrix_head(A);
    //print_csc_matrix(A);

    read_mm_matrix_csr(argv[2], &(B.M), &(B.N), &(B.nz), &(B.I), &(B.J), &(B.val));
    print_matrix_head(A);
    //print_csr_matrix(B);

    esc(A, B);

    return 0;
}
