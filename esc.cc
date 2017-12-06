#include <iostream>
#include <stdio.h>
#include <map>
#include "util.h"
#include "matrix_io.h"
#include "esc.h"
using namespace std;

// Scale ith row of matrix B
void scale_csr_row(Matrix mat, float scalar, int A_row, int A_col, map<pair<int, int>, float> *C) {
    // access ith row of B, i.e. ith col of A
    for(int i=mat.I[A_col]; i < mat.I[A_col+1]; i++) {
        int B_col = mat.J[i];
        double val = mat.val[i] * scalar;
        printf("(%d, %d): %f ", A_row, B_col, val);

        pair<int, int> p = pair<int, int>(A_row, B_col);
        if((*C).find(p) == (*C).end())
        	(*C)[p] = val;
        else
        	(*C)[p] += val;
    }
    printf("\n");
}

// Convert matrix stored in map to COO Matrix
Matrix map_to_coo_matrix(map<pair<int,int>, float> m) {
	Matrix mat;
	int size = m.size();
	mat.I = (int*)malloc(sizeof(int)*size);
	mat.J = (int*)malloc(sizeof(int)*size);
	mat.val = (float*)malloc(sizeof(float)*size);
	mat.M = mat.N =  mat.nz = size;

	int i = 0;
	map<pair<int, int>, float>::iterator it;
	for (it = m.begin(); it != m.end(); it++,i++) {
		printf("(%d, %d): %f\n", it->first.first, it->first.second, it->second);
		mat.I[i] = it->first.first;
		mat.J[i] = it->first.second;
		mat.val[i] = it->second;
	}
	return mat;
}


/**
 *  Matrix A: CSC Matrix
 *  Matrix B: CSR Matrix
 */
void esc(Matrix A, Matrix B) {
	map<pair<int, int>, float> C;

	for (int i = 0; i < A.N; i++) {
		// ith column of A
		for(int j = A.J[i]; j < A.J[i+1]; j++) {
			int row = A.I[j];
			float scalar = A.val[j];
			scale_csr_row(B, scalar, row, i, &C);
		}
	}

	printf("C:\n");
	Matrix res = map_to_coo_matrix(C);
	print_coo_matrix(res);
}

