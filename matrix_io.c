/* 
*   Matrix Market I/O 
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include "matrix_io.h"
#include "mmio.h"

void sort(int *col_idx, double *a, int start, int end);
void coo2csr_in(int n, int nz, int *i_idx, int *j_idx, double *a);
void csr2csc(int n, int m, int nz, int *row_start, int *col_idx, double *a, int *row_idx, int *col_start, double *csc_a);

void read_mm_matrix_coo(char *fname, int *M, int *N, int *nz, int **I, int **J, double **val) {

    int ret_code;
    MM_typecode matcode;

    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) exit(1);

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode) ) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }


    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) !=0)
        exit(1);


    /* reseve memory for matrices */
    *I = (int *) malloc(*nz * sizeof(int));
    *J = (int *) malloc(*nz * sizeof(int));
    *val = (double *) malloc(*nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (int i=0; i<*nz; i++) {
        fscanf(f, "%d %d %lg\n", &(*I)[i], &(*J)[i], &(*val)[i]);
        (*I)[i]--;  /* adjust from 1-based to 0-based */
        (*J)[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, *M, *N, *nz);
}


void read_mm_matrix_csr(char *fname, int *M, int *N, int *nz, int **I, int **J, double **val) {
    read_mm_matrix_coo(fname, M, N, nz, I, J, val);
    coo2csr_in(*M, *nz, *I, *J, *val);
}

void read_mm_matrix_csc(char *fname, int *M, int *N, int *nz, int **I, int **J, double **val) {
    int *csr_I, *csr_J;
    double *csr_val;
    read_mm_matrix_csr(fname, M, N, nz, &csr_I, &csr_J, &csr_val);
    csr2csc(*M, *N, *nz, csr_I, csr_J, csr_val, *I, *J, *val);
}


void sort(int *col_idx, double *a, int start, int end) {
    int i, j, it;
    double dt;
    for (i=end-1; i>start; i--)
        for(j=start; j<i; j++)
            if (col_idx[j] > col_idx[j+1]) {
                if (a) {
                    dt=a[j];
                    a[j]=a[j+1];
                    a[j+1]=dt;
                }
                it=col_idx[j];
                col_idx[j]=col_idx[j+1];
                col_idx[j+1]=it;
            }
}


/* converts COO format to CSR format, in-place 
   n is the number of rows. */
void coo2csr_in(int n, int nz, int *i_idx, int *j_idx, double *a) {
    int *row_start;
    int i, j;
    int init, i_next, j_next, i_pos;
    double dt, a_next;

    row_start = (int *)malloc((n+1)*sizeof(int));
    if (!row_start) {
        printf ("coo2csr_in: cannot allocate temporary memory\n");
        exit (1);
    }
    for (i=0; i<=n; i++) row_start[i] = 0;

    /* determine row lengths */
    for (i=0; i<nz; i++) row_start[i_idx[i]+1]++;

    for (i=0; i<n; i++) row_start[i+1] += row_start[i];

    for (init=0; init<nz;) {
        dt = a[init];
        i = i_idx[init];
        j = j_idx[init];
        i_idx[init] = -1;
        while (1) {
            i_pos = row_start[i];
            a_next = a[i_pos];
            i_next = i_idx[i_pos];
            j_next = j_idx[i_pos];

            a[i_pos] = dt;
            j_idx[i_pos] = j;
            i_idx[i_pos] = -1;
            row_start[i]++;
            if (i_next < 0) break;
            dt = a_next;
            i = i_next;
            j = j_next;
        }
        init++;
        while ((i_idx[init] < 0) && (init < nz))  init++;
    }


    /* shift back row_start */
    for (i=0; i<n; i++) i_idx[i+1] = row_start[i];
    i_idx[0] = 0;


    for (i=0; i<n; i++){
        sort(j_idx, a, i_idx[i], i_idx[i+1]);
    }

}


/*
   converts CSR format to CSC format, not in-place,
   if a == NULL, only pattern is reorganized.
   the size of matrix is m x n.
   */
void csr2csc(int m, int n, int nz, int *row_start, int *col_idx, double *a, int *row_idx, int *col_start, double *csc_a) {
    int i, j, k, l;
    int *ptr;

    for (i=0; i<=m; i++) col_start[i] = 0;

    /* determine column lengths */
    for (i=0; i<nz; i++) col_start[col_idx[i]+1]++;

    for (i=0; i<m; i++) col_start[i+1] += col_start[i];


    /* go through the structure once more. Fill in output matrix. */

    for (i=0, ptr=row_start; i<n; i++, ptr++)
        for (j=*ptr; j<*(ptr+1); j++){
            k = col_idx[j];
            l = col_start[k]++;
            row_idx[l] = i;
            if (a) csc_a[l] = a[j];
        }

    /* shift back col_start */
    for (i=m; i>0; i--) col_start[i] = col_start[i-1];

    col_start[0] = 0;
}
