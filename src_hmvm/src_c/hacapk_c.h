#ifndef _HACAPK_C_H
#define _HACAPK_C_H

// H matrix
typedef struct submatrix{
  int ltmtx, kt, ndl, ndt, nstrtl, nstrtt;
  double *a1, *a2, *a2t;
}submatrix;
typedef struct matrix{
  int nd;
  int nlf;
  int ktmax;
  int *nleavs;
  submatrix *submat;
  int b1, b2, b3;
  int e1, e2, e3;
}matrix;

typedef struct matrix2{
  int nd;
  int nlf;
  int ktmax;
  int *nleavs;
  int *ltmtx, *kt, *ndl, *ndt, *nstrtl, *nstrtt;
  int *a1, *a2; // head
  int len;
  double *rowmat, *rowmat_t; // data
  // 近似行列と密行列を分離してみる（確認用）
  int napprox, ndense;
  int *approx, *dense;
}matrix2;

// setup matrix
int loadmatrix(const char *fname, matrix*, matrix2*);
int dummymatrix(matrix *mat, matrix2 *mat2, int ndl, int ndt, int n);

// seq
void hmvm_seq(matrix mat, matrix2 mat2, double *b);
void hmvm_seq_bench(matrix mat, matrix2 mat2, double *b);
// omp
void hmvm_omp(matrix mat, matrix2 mat2, double *b);
void hmvm_omp_bench(matrix mat, matrix2 mat2, double *b);

// mkl
#ifdef __INTEL_COMPILER
void hmvm_blas_p(matrix mat, matrix2 mat2, double *b);
void hmvm_blas_p_bench(matrix mat, matrix2 mat2, double *b);
void hmvm_blas_s(matrix mat, matrix2 mat2, double *b);
void hmvm_blas_s_bench(matrix mat, matrix2 mat2, double *b);
void hmvm_cblas_p(matrix mat, matrix2 mat2, double *b);
void hmvm_cblas_p_bench(matrix mat, matrix2 mat2, double *b);
void hmvm_cblas_s(matrix mat, matrix2 mat2, double *b);
void hmvm_cblas_s_bench(matrix mat, matrix2 mat2, double *b);
#endif

// cuda
#ifdef _USE_CUDA
#ifdef __cplusplus
extern "C"{
#endif
void hmvm_cuda1(matrix2 mat2, double *b, int kernel);
#ifdef __cplusplus
}
#endif
#endif

#endif
