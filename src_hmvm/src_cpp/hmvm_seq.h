#ifndef HMVM_SEQ_H
#define HMVM_SEQ_H

template <class T>
void hmvm_seq(matrix<T> mat, matrix2<T> mat2, T *b, int dump_result);
template <class T>
void hmvm_seq_bench(matrix<T> mat, matrix2<T> mat2, T *b);

#endif
