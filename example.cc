#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include <vector>
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
    Matrix tmpC;

    tmpC.I = (int*)malloc(A.nz*B.nz*sizeof(int));
    tmpC.J = (int*)malloc(A.nz*B.nz*sizeof(int));
    tmpC.val = (double*)malloc(A.nz*B.nz*sizeof(double));

    // expression
    int k = 0;
    for (int i=0;i<A.nz;i++) {
        for(int j=0;j<B.nz;j++) {
            if (A.J[i] == B.I[j]) {
                tmpC.I[k] = A.I[i];
                tmpC.J[k] = B.J[j];
                tmpC.val[k] = A.val[i] * B.val[j];
                k++;
            }
        }
    }
    tmpC.nz = k;
    printf("tmpC:\n");
    print_coo_matrix(tmpC);

    // sort
    std::vector<std::tuple<int, int, double>> vals;
    for(int i = 0; i < tmpC.nz; i++)
        vals.push_back(std::make_tuple(tmpC.I[i], tmpC.J[i], tmpC.val[i]));
    std::sort(vals.begin(), vals.end(), 
        [](tuple<int, int, double> const &t1, tuple<int, int, double> const &t2) {
            return (get<0>(t1) < get<0>(t2) && get<1>(t1) < get<1(t2) ); // or use a custom compare function
    });


    Matrix C;
    C.I = (int*)malloc(tmpC.nz*sizeof(int));
    C.J = (int*)malloc(tmpC.nz*sizeof(int));
    C.val = (double*)malloc(tmpC.nz*sizeof(double));
    // compression
    k = 0;
    for(int i = 0; i < tmpC.nz; i++) {
        C.I[k] = tmpC.I[i];
        C.J[k] = tmpC.J[i];
        C.val[k] = tmpC.val[i];
        for(int j = i+1; j < C.nz; j++) {
            if(tmpC.I[i] == tmpC.I[j] && tmpC.J[i] == tmpC.J[j]) {
                k++;
            } else {
                break;
            }
        }
    }
    C.nz = k;

}


int main(int argc, char *argv[]) {
    Matrix A, B;   // Matrix is defined in matrix_io.h

    if (argc < 3) {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    
    read_mm_matrix_coo(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_coo_matrix(A);

    read_mm_matrix_coo(argv[2], &(B.M), &(B.N), &(B.nz), &(B.I), &(B.J), &(B.val));
    print_coo_matrix(B);

    esc_multiplication(A, B);

    //read_mm_matrix_csr(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    //print_csr_matrix(A);

    //read_mm_matrix_csc(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    //print_csc_matrix(A);
    return 0;
}