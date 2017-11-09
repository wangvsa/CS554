void read_mm_matrix (char *fn, int *m, int *n, int *nz, int **i_idx, int **j_idx, double **a);
void write_csr (char *fn, int m, int n, int nz, int *row_start, int *col_idx, double *a);
