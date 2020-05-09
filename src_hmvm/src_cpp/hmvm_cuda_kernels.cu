// -*- c++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>
#include <cuda_runtime_api.h>
#include <cooperative_groups.h>

#include "hacapk.h"

// できれば重複計算していないかチェックもしたい

#define CHECK_DO(act,msg) {ret=act; if(ret!=cudaSuccess){printf("%s failed\n",msg);exit(-1);};}






template
void hmvm_cuda_hybrid4_proxy<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat, float *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int mul, int atomic, float *v, float *b, int nd, char *fname, int bench);
template
void hmvm_cuda_hybrid4_proxy<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat, double *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int mul, int atomic, double *v, double *b, int nd, char *fname, int bench);
// ######## ######## ######## ######## ######## ######## ######## ########
#endif
