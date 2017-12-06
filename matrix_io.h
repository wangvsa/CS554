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
*      float precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading floats, otherwise errors will occur.
*/


#ifndef MATRIX_IO_H
#define MATRIX_IO_H

typedef struct type_matrix {
	int M, N, nz;	// number of rows, columns, and non-zeros
	int *I, *J;		// row pointers, column indices
	float *val;	// all non-zero values
} Matrix;

/* reads matrix market format (coordinate) and returns csr format 
    int: fname, filename
    out: M, N, nz, I, J, val*/
void read_mm_matrix_csr(char *fname, int *M, int *N, int *nz, int **I, int **J, float **val);

/* reads matrix market format (coordinate) and returns coo format
    int: fname, filename
    out: M, N, nz, I, J, val*/
void read_mm_matrix_coo(char *fname, int *M, int *N, int *nz, int **I, int **J, float **val);

/* reads matrix market format (coordinate) and returns coo format
    int: fname, filename
    out: M, N, nz, I, J, val*/
void read_mm_matrix_csc(char *fname, int *M, int *N, int *nz, int **I, int **J, float **val);

#endif