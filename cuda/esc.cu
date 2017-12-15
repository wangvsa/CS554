#include <iostream>
#include <stdio.h>
#include <chrono>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "../matrix_io.h"
#include "../util.h"
#include "esc.h"
#include "timer.h"

using namespace std;
using Clock=std::chrono::high_resolution_clock;


#define HASH_TABLE_SIZE 268435456
#define GRID_SIZE 10000                // the number of blocks
#define BLOCK_SIZE 32                  // the number of threads in each block
int const INVALID_KEY = 0;
int BLOCK_CONFLICT_SIZE = 100;        // the size of confict array for each block


// For storing collision results
typedef struct ConflictArray_t {
    int *counter;     // each block has its own conflict counter
    unsigned int *key;
    float *val;
} ConflictArray;


// Hash function to convert coordinate(i, j) to a single integer
__device__ inline unsigned int hashcode(int i, int j, int COLS) { return i * COLS + j; }
inline int get_i(unsigned int key, int COLS) { return (key-1) / COLS; }
inline int get_j(unsigned int key, int COLS) { return (key-1) % COLS; }


// Every thread in the same block will run this method with different scalar
// So mat better be stored in shared memory
__device__
inline void scale_csr_row(Matrix mat, float scalar, int A_row, int A_col, unsigned int *dev_C_key, float *dev_C_val, ConflictArray conflict_array, int BLOCK_CONFLICT_SIZE) {
    // access ith row of B, i.e. ith col of A
    for(int i=mat.I[A_col]; i < mat.I[A_col+1]; i++) {
        int B_col = mat.J[i];
        float val = mat.val[i] * scalar;

        // p is a simple hash code for (A_row, B_col) pair
        unsigned int h = hashcode(A_row, B_col, mat.M) + 1;      // leave zero as the INVALID_KEY
        unsigned int index = h & (HASH_TABLE_SIZE-1);

        unsigned int old_h = atomicCAS(&(dev_C_key[index]), INVALID_KEY, h);
        if(old_h == INVALID_KEY) {
            atomicAdd(&(dev_C_val[index]), val);
        } else {
            if(old_h == h) {
                atomicAdd(&(dev_C_val[index]), val);
            } else {
                int counter = atomicAdd(&(conflict_array.counter[blockIdx.x]), 1);
                atomicExch(&(conflict_array.key[blockIdx.x*BLOCK_CONFLICT_SIZE+counter]), h);
                atomicExch(&(conflict_array.val[blockIdx.x*BLOCK_CONFLICT_SIZE+counter]), val);
            }
        }
        //printf("h:%d, index:%d, (%d %d %f %f)\n", h, index, A_row, B_col, val, dev_C_val[index]);
    }
}

__global__
void cuda_multiplication(Matrix A, Matrix B, unsigned int *dev_C_key, float *dev_C_val, ConflictArray conflict_array, int BLOCK_CONFLICT_SIZE) {
    /*
     *  go through all columns of csc matrix A
     *  each block is responsible for processing one column of A
     */
    int block_load = A.N / gridDim.x;
    if(A.N % gridDim.x != 0) block_load = block_load + 1;

    for (int i = blockIdx.x*block_load; i < (blockIdx.x+1)*block_load; i++) {
        if(i >= A.N) return;
        // get all values in ith column of A
        // each thread in the same block get a portion of tasks
        int load = (A.J[i+1] - A.J[i]) / blockDim.x;
        if((A.J[i+1] - A.J[i]) % blockDim.x !=0) load = load + 1;
        for(int j = 0; j < load; j++) {
            int k = A.J[i] + threadIdx.x * load + j;
            if(k < A.J[i+1]) {
                int row = A.I[k];
                float scalar = A.val[k];
                scale_csr_row(B, scalar, row, i, dev_C_key, dev_C_val, conflict_array, BLOCK_CONFLICT_SIZE);
            }
        }
    }

    // Debug, print conflict for each block
    /*
    if(blockIdx.x == 0 && threadIdx.x == 0) {
        unsigned int total = 0;
        for(int i = 0; i < gridDim.x; i++) {
            printf("conflict of block %d: %d\n", i, conflict_array.counter[i]);
            total += conflict_array.counter[i];
        }
        printf("total conflict: %d\n", total);
    }
    */
}

inline void set_coo_element(Matrix mat, int index, int i, int j, float val) {
    mat.I[index] = i;
    mat.J[index] = j;
    mat.val[index] = val;
}
void hashtable_to_coo(int COLS, unsigned int *key1, float *val1, thrust::host_vector<unsigned int> key2, thrust::host_vector<float> val2) {
    int nnz = 0;
    for(int i = 0; i < HASH_TABLE_SIZE; i++)
        if(key1[i] != INVALID_KEY) nnz++;
    printf("nnz of C: %d\n", nnz);
    for(int i = 0; i<key2.size(); i++)
        if(key2[i] != INVALID_KEY) nnz++;
    printf("nnz of C: %d\n", nnz);

    Matrix C;
    C.nz = nnz;
    C.I = new int[nnz];
    C.J = new int[nnz];
    C.val = new float[nnz];
    int t = 0;
    for(int i = 0; i < HASH_TABLE_SIZE; i++) {
        unsigned int index = key1[i];
        if(index != INVALID_KEY)
            set_coo_element(C, t++, get_i(index, COLS), get_j(index, COLS), val1[i]);
    }
    for(int i = 0; i < key2.size(); i++) {
        unsigned int index = key2[i];
        if(index != INVALID_KEY)
            set_coo_element(C, t++, get_i(index, COLS), get_j(index, COLS), val2[i]);
    }
}

void esc(Matrix A, Matrix B, int conflict_size) {
    BLOCK_CONFLICT_SIZE = conflict_size;
    printf("\n====================================================\n");
    printf("blocks: %d, threads: %d\nhashtable size: %d, conflict array size: %d\n", GRID_SIZE, BLOCK_SIZE, HASH_TABLE_SIZE, BLOCK_CONFLICT_SIZE);
    printf("====================================================\n\n");
    esc(A, B);
}

// Matrix A: csc format; Matrix B: csr format
void esc(Matrix A, Matrix B) {
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
    unsigned int *dev_C_key;
    float *dev_C_val;
    cudaMalloc((void**)&dev_C_key, sizeof(unsigned int)*HASH_TABLE_SIZE);
    cudaMalloc((void**)&dev_C_val, sizeof(float)*HASH_TABLE_SIZE);
    cudaMemset(dev_C_val, 0, sizeof(float)*HASH_TABLE_SIZE);

    // Allocate device memory for collision results
    ConflictArray conflict_array;
    cudaMalloc((void**)&conflict_array.counter, sizeof(int)*GRID_SIZE);
    cudaMalloc((void**)&conflict_array.key, sizeof(unsigned int)*GRID_SIZE * BLOCK_CONFLICT_SIZE);  // each block has 1M for collision results
    cudaMalloc((void**)&conflict_array.val, sizeof(float)*GRID_SIZE * BLOCK_CONFLICT_SIZE);
    printf("time for allocating memory: %f milliseconds\n", t1.milliseconds_elapsed());

    timer t2;
    cuda_multiplication<<<GRID_SIZE, BLOCK_SIZE>>>(dev_A, dev_B, dev_C_key, dev_C_val, conflict_array, BLOCK_CONFLICT_SIZE);
    printf("time for SpMM: %f milliseconds\n", t2.milliseconds_elapsed());

    // Copy back the hash table
    unsigned int *C_key = new unsigned int[HASH_TABLE_SIZE];
    float *C_val = new float[HASH_TABLE_SIZE];
    cudaMemcpy(C_key, dev_C_key, HASH_TABLE_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaMemcpy(C_val, dev_C_val, HASH_TABLE_SIZE*sizeof(float), cudaMemcpyDeviceToHost);

    timer t3;
    thrust::device_ptr<unsigned int> key_ptr(conflict_array.key);
    thrust::device_ptr<float> val_ptr(conflict_array.val);
    thrust::device_vector<unsigned int> output_keys(GRID_SIZE*BLOCK_CONFLICT_SIZE);
    thrust::device_vector<float> output_vals(GRID_SIZE*BLOCK_CONFLICT_SIZE);
    thrust::sort_by_key(key_ptr, key_ptr+GRID_SIZE*BLOCK_CONFLICT_SIZE, val_ptr);
    thrust::reduce_by_key(key_ptr, key_ptr+GRID_SIZE*BLOCK_CONFLICT_SIZE, val_ptr, output_keys.begin(), output_vals.begin(), thrust::equal_to<unsigned int>(), thrust::plus<float>());
    printf("time for reduction: %f milliseconds\n", t3.milliseconds_elapsed());

    thrust::host_vector<unsigned int> key2 = output_keys;
    thrust::host_vector<float> val2 = output_vals;
    hashtable_to_coo(A.N, C_key, C_val, key2, val2);

    cudaFree(dev_A.I);
    cudaFree(dev_A.J);
    cudaFree(dev_A.val);
    cudaFree(dev_B.I);
    cudaFree(dev_B.J);
    cudaFree(dev_B.val);
    cudaFree(dev_C_key);
    cudaFree(dev_C_val);
    cudaFree(conflict_array.counter);
    cudaFree(conflict_array.key);
    cudaFree(conflict_array.val);
    delete C_key;
    delete C_val;
}
