#ifndef _HMVM_CUDA_H
#define _HMVM_CUDA_H

template <class T> void hmvm_cuda0(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench);
template <class T> void hmvm_cuda1(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench);
template <class T> void hmvm_cuda2(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench);
template <class T> void hmvm_cuda3(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench);
//template <class T> void hmvm_cuda4(matrix2<T> *mat2, T *b, int kernel, int dump_result);

#endif
