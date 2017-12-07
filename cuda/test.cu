#include <iostream>
#include <stdio.h>
#include <chrono>
#include "../matrix_io.h"
#include "../util.h"
#include "timer.h"

using namespace std;
using Clock=std::chrono::high_resolution_clock;

// Define key(i, j), convert coordinate(i, j) to a size_t value
inline size_t key(int i, int j) {return (size_t) i << 32 | (unsigned int) j;}
inline int get_first(size_t C) { return C>>32; }
inline int get_second(size_t C) { return C & 0xFFFFFFFF; }


// Every thread in the same block will work on this method with different scalar
// So mat should be stored in shared memory
__device__
void scale_csr_row(Matrix mat, float scalar, int A_row, int A_col) {
    // access ith row of B, i.e. ith col of A
    for(int i=mat.I[A_col]; i < mat.I[A_col+1]; i++) {
        int B_col = mat.J[i];
        double val = mat.val[i] * scalar;
        //size_t p = key(A_row, B_col);
        /*
        if((*C).find(p) == (*C).end()) {
            (*C)[p] = val;
        } else {
            (*C)[p] += val;
        }
        */
    }
}

__global__
void cuda_multiplication(Matrix A, Matrix B) {
    // go through all columns of csc matrix A
    // each block is responsible for processing one column of A
    for (int i = 0; i < A.N; i++) {
        // get all values in ith column of A
        for(int j = A.J[i]; j < A.J[i+1]; j++) {
            int row = A.I[j];
            float scalar = A.val[j];
            scale_csr_row(B, scalar, row, i);
        }
    }
}


int main(int argc, char *argv[]) {
    Matrix A, B;
    if (argc < 3) {
        printf("Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }

    read_mm_matrix_csc(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_matrix_head(A);
    read_mm_matrix_csr(argv[2], &(B.M), &(B.N), &(B.nz), &(B.I), &(B.J), &(B.val));
    print_matrix_head(B);

    timer t;
    cuda_multiplication<<<10, 128>>>(A, B);
    printf("time:%f milliseconds\n", t.milliseconds_elapsed());

    return 0;
}
