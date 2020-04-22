/* -*- C++ -*- */
#ifndef _HMVM_CUDA_KERNELS_H
#define _HMVM_CUDA_KERNELS_H

template <class T>
__global__ void hmvm_cuda_seq
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense);

template <class T>
__global__ void hmvm_cuda_block
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense);

template <class T, int div, int atomic>
__global__ void hmvm_cuda_hybrid1
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense);

template <class T>
void hmvm_cuda_hybrid2_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, T *v, T *b, int nd, char *fname, int bench);

template <class T>
void hmvm_cuda_hybrid3_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, T *v, T *b, int nd, char *fname, int bench);

#endif
