#include <iostream>
#include <stdio.h>
#include <map>
#include <chrono>
#include "util.h"
#include "matrix_io.h"
#include "esc.h"

using namespace std;
using Clock=std::chrono::high_resolution_clock;

// Scale ith row of matrix B
void scale_csr_row(Matrix mat, float scalar, int A_row, int A_col, map<pair<int, int>, float> *C) {
    // access ith row of B, i.e. ith col of A
    for(int i=mat.I[A_col]; i < mat.I[A_col+1]; i++) {
        int B_col = mat.J[i];
        double val = mat.val[i] * scalar;

        pair<int, int> p = pair<int, int>(A_row, B_col);
        if((*C).find(p) == (*C).end())
        	(*C)[p] = val;
        else
        	(*C)[p] += val;
    }
}

// Convert matrix stored in map to COO Matrix
Matrix map_to_coo_matrix(map<pair<int,int>, float> m, int ROWS, int COLS) {
	Matrix mat;
	int size = m.size();
	mat.M = ROWS;
	mat.N = COLS;
	mat.nz = size;
	mat.I = (int*)malloc(sizeof(int)*size);
	mat.J = (int*)malloc(sizeof(int)*size);
	mat.val = (float*)malloc(sizeof(float)*size);

	int i = 0;
	map<pair<int, int>, float>::iterator it;
	for (it = m.begin(); it != m.end(); it++,i++) {
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
	auto start = Clock::now();
	map<pair<int, int>, float> C;

	for (int i = 0; i < A.N; i++) {
		// ith column of A
		for(int j = A.J[i]; j < A.J[i+1]; j++) {
			int row = A.I[j];
			float scalar = A.val[j];
			scale_csr_row(B, scalar, row, i, &C);
		}
	}
	auto end = Clock::now();

	Matrix res = map_to_coo_matrix(C, A.M, B.N);
	print_matrix_head(res);
	//print_coo_matrix(res);

	std::chrono::duration<double, milli> fp_ms = end - start;
	cout<<"time: "<<fp_ms.count()<<" milliseconds"<<endl;
}

