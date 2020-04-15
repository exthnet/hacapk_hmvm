#ifndef _HACAPK_C_HPP
#define _HACAPK_C_HPP

// H matrix
template <class T>
struct submatrix{
  int ltmtx, kt, ndl, ndt, nstrtl, nstrtt;
  T *a1, *a2, *a2t;
};

template <class T>
struct matrix{
  int nd;
  int nlf;
  int ktmax;
  int *nleavs;
  submatrix<T> *submat;
  int b1, b2, b3;
  int e1, e2, e3;
};

// H matrix for GPU
template <class T>
struct matrix2{
  int nd;
  int nlf;
  int ktmax;
  int len;
  int *nleavs;
  int *ltmtx, *kt, *ndl, *ndt, *nstrtl, *nstrtt;
  int *a1, *a2; // head
  T *rowmat, *rowmat_t; // data
  // 近似行列と密行列を分離してみる（確認用）
  int napprox, ndense;
  int *approx, *dense;
};

// setup matrix
int loadHmatrix(const char *fname, matrix<double>*, matrix2<double>*);
int loadHmatrix2(const char *fname, matrix2<double>*);
int convertD2F(matrix<float> *matF, matrix2<float> *matF2, const matrix<double> *matD, const matrix2<double> *matD2);

// mkl
#ifdef __INTEL_COMPILER
template <class T>void hmvm_blas_p(matrix<T> *mat, matrix2<T> *mat2, T *b);
template <class T>void hmvm_blas_p_bench(matrix<T> *mat, matrix2<T> *mat2, T *b);
template <class T>void hmvm_blas_s(matrix<T> *mat, matrix2<T> *mat2, T *b);
template <class T>void hmvm_blas_s_bench(matrix<T> *mat, matrix2<T> *mat2, T *b);
template <class T>void hmvm_cblas_p(matrix<T> *mat, matrix2<T> *mat2, T *b);
template <class T>void hmvm_cblas_p_bench(matrix<T> *mat, matrix2<T> *mat2, T *b);
template <class T>void hmvm_cblas_s(matrix<T> *mat, matrix2<T> *mat2, T *b);
template <class T>void hmvm_cblas_s_bench(matrix<T> *mat, matrix2<T> *mat2, T *b);
#endif

// cuda
#ifdef _USE_CUDA
//template <class T>void hmvm_cuda1(matrix2<T> *mat2, T *b, int kernel);
#endif

#endif

