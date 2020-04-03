#ifndef HMVM_OMP_H
#define HMVM_OMP_H

template <class T>
void hmvm_omp(matrix<T> mat, matrix2<T> mat2, T *b, int dump_result);
template <class T>
void hmvm_omp_bench(matrix<T> mat, matrix2<T> mat2, T *b);

#endif
