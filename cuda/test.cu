#include <iostream>
#include <stdio.h>
#include <chrono>
#include "../matrix_io.h"
#include "../util.h"
#include "timer.h"

using namespace std;
using Clock=std::chrono::high_resolution_clock;

#define HASH_TABLE_SIZE 32


// Define key(i, j), convert coordinate(i, j) to a size_t value
__device__ inline int key(int i, int j, int COLS) { return i * COLS + j;}
inline int get_i(int index, int COLS) { return index / COLS; }
inline int get_j(int index, int COLS) { return index % COLS; }


// Every thread in the same block will run this method with different scalar
// So mat better be stored in shared memory
__device__
inline void scale_csr_row(Matrix mat, float scalar, int A_row, int A_col, int *dev_C_key, float *dev_C_val) {
    // access ith row of B, i.e. ith col of A
    for(int i=mat.I[A_col]; i < mat.I[A_col+1]; i++) {
        int B_col = mat.J[i];
        float val = mat.val[i] * scalar;

        // p is a simple hash code for (A_row, B_col) pair
        int p = key(A_row, B_col, mat.M);
        int index = p & (HASH_TABLE_SIZE-1);

        dev_C_key[index] = p;
        atomicAdd(&(dev_C_val[index]), val);

        //printf("p:%d, index:%d, (%d %d %f %f)\n", p, index, A_row, B_col, val, dev_C_val[index]);
    }
}

__global__
void cuda_multiplication(Matrix A, Matrix B, int *dev_C_key, float *dev_C_val) {
    /*
     *  go through all columns of csc matrix A
     *  we set blockDim.x == A.N
     *  so each block is responsible for processing one column of A
     */

    //for (int i = 0; i < A.N; i++) {
    int i = blockIdx.x;
        // get all values in ith column of A
        // each thread in the same block get a portion of tasks
        int load = (A.J[i+1] - A.J[i]) / blockDim.x;
        if((A.J[i+1] - A.J[i]) % blockDim.x !=0) load = load + 1;
        for(int j = 0; j < load; j++) {
            int k = A.J[i] + threadIdx.x * load + j;
            if(k < A.J[i+1]) {
                int row = A.I[k];
                float scalar = A.val[k];
                scale_csr_row(B, scalar, row, i, dev_C_key, dev_C_val);
            }
        }
    //}
}


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


    // Allocate device memory and copy matrices from the host to the device
    timer t1;
    Matrix dev_A, dev_B;
    dev_A.M = A.M; dev_A.N = A.N; dev_A.nz = A.nz;
    cudaMalloc((void**)&dev_A.I, sizeof(int)*(A.nz));
    cudaMalloc((void**)&dev_A.J, sizeof(int)*(A.N+1));
    cudaMalloc((void**)&dev_A.val, sizeof(float)*(A.nz));
    cudaMemcpy(dev_A.I, A.I, A.nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_A.J, A.J, (A.N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_A.val, A.val, A.nz*sizeof(float), cudaMemcpyHostToDevice);
    dev_B.M = B.M; dev_B.N = B.N; dev_B.nz = B.nz;
    cudaMalloc((void**)&dev_B.I, sizeof(int)*(B.M+1));
    cudaMalloc((void**)&dev_B.J, sizeof(int)*(B.nz));
    cudaMalloc((void**)&dev_B.val, sizeof(float)*(B.nz));
    cudaMemcpy(dev_B.I, B.I, (B.M+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B.J, B.J, B.nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B.val, B.val, B.nz*sizeof(float), cudaMemcpyHostToDevice);

    // Allocate device memory for has table
    int *dev_C_key;
    float *dev_C_val;
    cudaMalloc((void**)&dev_C_key, sizeof(int)*HASH_TABLE_SIZE);
    cudaMemset(dev_C_key, -1, sizeof(int)*HASH_TABLE_SIZE);
    cudaMalloc((void**)&dev_C_val, sizeof(float)*HASH_TABLE_SIZE);
    cudaMemset(dev_C_val, 0, sizeof(float)*HASH_TABLE_SIZE);
    printf("time for allocating memory:%f milliseconds\n", t1.milliseconds_elapsed());

    timer t2;
    cuda_multiplication<<<A.N, 32>>>(dev_A, dev_B, dev_C_key, dev_C_val);
    printf("time for SpMM: %f milliseconds\n", t2.milliseconds_elapsed());

    // Copy back the hash table
    int *C_key = new int[HASH_TABLE_SIZE];
    float *C_val = new float[HASH_TABLE_SIZE];
    cudaMemcpy(C_key, dev_C_key, HASH_TABLE_SIZE*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(C_val, dev_C_val, HASH_TABLE_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
    for(int i = 0; i < HASH_TABLE_SIZE; i++) {
        int index = C_key[i];
        if(index != -1)
            printf("%d: (%d %d %f)\n", index, get_i(index, A.N), get_j(index, A.N), C_val[i]);
    }



    cudaFree(dev_A.I);
    cudaFree(dev_A.J);
    cudaFree(dev_A.val);
    cudaFree(dev_B.I);
    cudaFree(dev_B.J);
    cudaFree(dev_B.val);
    cudaFree(dev_C_key);
    cudaFree(dev_C_val);
    delete C_key;
    delete C_val;
    return 0;
}
