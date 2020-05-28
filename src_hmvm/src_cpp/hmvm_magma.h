#ifndef _HMVM_MAGMA_H
#define _HMVM_MAGMA_H

template <class T> void hmvm_magma(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench);
template <class T> void hmvm_magma_batched1(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench);
template <class T> void hmvm_magma_batched2(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench);

#endif
