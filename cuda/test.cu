#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <algorithm>
#include <cstdlib>
#include "../matrix_io.h"
#include "../util.h"

typedef struct DeviceMatrix_t {
    int M, N, nz;       // number of rows, number of cloumns, and number of non-zeros
    thrust::device_vector<int> I;
    thrust::device_vector<int> J;
    thrust::device_vector<double> val;
    DeviceMatrix_t(int _M, int _N, int _nz, int *_I, int *_J, double *_val) {
        M = _M;
        N = _N;
        nz = _nz;
        thrust::copy(_I, _I+M+1, I.begin());
        thrust::copy(_J, _J+N+1, J.begin());
        thrust::copy(_val, _val+nz, val.begin());
    }
} DeviceMatrix;


void test()
{
    // generate random data serially
    thrust::host_vector<int> h_vec(100);
    std::generate(h_vec.begin(), h_vec.end(), rand);
    // transfer to device and compute sum
    thrust::device_vector<int> d_vec = h_vec;
    int x = thrust::reduce(d_vec.begin(), d_vec.end(), 0, thrust::plus<int>());
    printf("x:%d", x);
}

int main(int argc, char *argv[]) {
    Matrix A, B;
    if (argc < 3) {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    read_mm_matrix_csc(argv[1], &(A.M), &(A.N), &(A.nz), &(A.I), &(A.J), &(A.val));
    print_csc_matrix(A);

    read_mm_matrix_csr(argv[2], &(B.M), &(B.N), &(B.nz), &(B.I), &(B.J), &(B.val));
    print_csr_matrix(B);
    return 0;
}
