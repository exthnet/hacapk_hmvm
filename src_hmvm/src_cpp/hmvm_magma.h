#ifndef _HMVM_MAGMA_H
#define _HMVM_MAGMA_H

template <class T> void hmvm_magma(matrix2<T> *mat2, T *b, int kernel, int dump_result);
template <class T> void hmvm_magma_batched(matrix2<T> *mat2, T *b, int kernel, int dump_result);
template <class T> void hmvm_magma_batched2(matrix2<T> *mat2, T *b, int kernel, int dump_result);

#endif
