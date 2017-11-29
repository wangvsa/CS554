#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include <vector>
#include <algorithm>
#include "matrix_io.h"
#include "util.h"
using namespace std;

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
    vector<tuple<int, int, double>> vals;
    for(int i = 0; i < tmpC.nz; i++)
        vals.push_back(make_tuple(tmpC.I[i], tmpC.J[i], tmpC.val[i]));
    sort(vals.begin(), vals.end(), 
        [](tuple<int, int, double> const &t1, tuple<int, int, double> const &t2) -> bool {
            return (get<0>(t1)<get<0>(t2)) || (get<1>(t1)<get<1>(t2)); // or use a custom compare function
    });
	vector<tuple<int,int,double>>::iterator it;
	printf("Sort:\n");
	for(it = vals.begin(); it != vals.end(); it++)    {
        printf("%d %d %f\n", get<0>(*it),  get<1>(*it), get<2>(*it));
	}

    Matrix C;
    C.I = (int*)malloc(tmpC.nz*sizeof(int));
    C.J = (int*)malloc(tmpC.nz*sizeof(int));
    C.val = (double*)malloc(tmpC.nz*sizeof(double));
    // compression
    k = 0;
    for(int i = 0; i < tmpC.nz; i++) {
        C.I[k] = get<0>(vals[i]);
        C.J[k] = get<1>(vals[i]);
        C.val[k] = get<2>(vals[i]);
        for(int j = i+1; j < tmpC.nz; j++) {
            if( C.I[k] == get<0>(vals[j]) && C.J[k] == get<1>(vals[j]) ) {
                C.val[k] += get<2>(vals[j]);
            } else {
                i = j - 1;
                break;
            }
        }
        k++;
    }
    C.nz = k;
    printf("C:\n");
    print_coo_matrix(C);

}


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

    read_mm_matrix_csr(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_csr_matrix(A);

    read_mm_matrix_csc(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_csc_matrix(A);

    return 0;
}
