#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <vector>
#include <omp.h>
#include <mkl.h>
#include "util.h"
#include "matrix_io.h"
#include "esc.h"

using namespace std;

// Define key(i, j), convert coordinate(i, j) to a size_t value
inline int key(int i, int j, int COLS) {return i * COLS + j;}
inline int get_first(size_t C) { return C>>32; }
inline int get_second(size_t C) { return C & 0xFFFFFFFF; }


// Scale ith row of matrix B
inline void scale_csr_row(Matrix mat, float scalar, int A_row, int A_col, unordered_map<int, float> *C) {
    // access ith row of B, i.e. ith col of A
    for(int i=mat.I[A_col]; i < mat.I[A_col+1]; i++) {
        int B_col = mat.J[i];
        double val = mat.val[i] * scalar;
        size_t p = key(A_row, B_col, mat.N);
        // The operation of map takes most of computation time
        /*
        if((*C).find(p) == (*C).end()) {
            (*C)[p] = val;
        } else {
            (*C)[p] += val;
        }
        */
        (*C)[p] = val;
    }
}

// Convert matrix stored in map to COO Matrix
Matrix map_to_coo_matrix(unordered_map<size_t, float> m, int ROWS, int COLS) {
	Matrix mat;
	int size = m.size();
	mat.M = ROWS;
	mat.N = COLS;
	mat.nz = size;
	mat.I = (int*)malloc(sizeof(int)*size);
	mat.J = (int*)malloc(sizeof(int)*size);
	mat.val = (float*)malloc(sizeof(float)*size);

	int i = 0;
	unordered_map<size_t, float>::iterator it;
	for (it = m.begin(); it != m.end(); it++,i++) {
		mat.I[i] = get_first(it->first);
		mat.J[i] = get_second(it->first);
		mat.val[i] = it->second;
	}
	return mat;
}


/**
 *  Matrix A: CSC Matrix
 *  Matrix B: CSR Matrix
 *  Use unordered_map to avoid sort and compression
 */
void esc(Matrix A, Matrix B) {
    double t1 = dsecnd();

    int const T = 8;
    // each thread will operate on its own map
    // and we'll combine all maps at the end
	unordered_map<int, float> maps[T];

    omp_set_num_threads(T);
    #pragma omp parallel for
	for (int i = 0; i < A.N; i++) {
        int tid = omp_get_thread_num();
		// ith column of A
		for(int j = A.J[i]; j < A.J[i+1]; j++) {
			int row = A.I[j];
			float scalar = A.val[j];
			scale_csr_row(B, scalar, row, i, &(maps[tid]));
		}
	}

    double t2 = dsecnd();

	//Matrix res = map_to_coo_matrix(C, A.M, B.N);
	//print_matrix_head(res);
	//print_coo_matrix(res);

	cout<<"time: "<<t2-t1<<" seconds"<<endl;
}

