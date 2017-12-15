#include <iostream>
#include <stdio.h>
#include <chrono>
#include "../matrix_io.h"
#include "../util.h"
#include "esc.h"

using namespace std;
using Clock=std::chrono::high_resolution_clock;


int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s [martix-market-filename] [matrix-market-filename]\n", argv[0]);
        exit(1);
    }


    // Read matrix from matrix-market file
    Matrix A, B;
    read_mm_matrix_csc(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_matrix_head(A);
    read_mm_matrix_csr(argv[2], &(B.M), &(B.N), &(B.nz), &(B.I), &(B.J), &(B.val));
    print_matrix_head(B);

    int conflict_size = 100;     // size of collision array
    if(argc == 4)
        conflict_size = atoi(argv[3]);
    esc(A, B, conflict_size);
}
