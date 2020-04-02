// -*- C++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <cuda_runtime_api.h>

#include "hacapk_c.h"

__global__ void hmvm_cudaD_kernel000000
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat, int ndense, int *dense);

__global__ void hmvm_cudaD
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);
__global__ void hmvm_cudaD_block
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);

template <int div>
__global__ void hmvm_cudaD_kernel00dd00
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat, int ndense, int *dense);

#define CHECK_DO(act,msg) {ret=act; if(ret!=cudaSuccess){printf("%s failed\n",msg);exit(-1);};}

void hmvm_cuda1(matrix2 mat2, double *b, int kernel)
{
  const int L=5, M=5;
  FILE *F;
  matrix2 d_sm;
  int i, l, nd = mat2.nd;
  double d1, d2, dtimes[L+M], dmin, dmax, davg1, davg2;
  double *v=NULL, *tmp=NULL, *zero;
  double *d_b, *d_v, *d_zaut, *d_zbut;
  int ip;
  int len, offset=0;
  cudaError_t ret;
  printf("hmvm_cuda1: begin\n");
  v=(double*)malloc(sizeof(double)*mat2.nd);
  tmp=(double*)malloc(sizeof(double)*mat2.ktmax);
  zero = (double*)malloc(sizeof(double)*mat2.ktmax);
  for(i=0;i<nd;i++){
	v[i] = 0.0;
  }
  for(i=0;i<mat2.ktmax;i++){
	zero[i] = 0.0;
  }
  CHECK_DO(cudaMalloc((void**)&d_zaut, sizeof(double)*mat2.nd),"cudaMalloc z_aut");
  CHECK_DO(cudaMalloc((void**)&d_zbut, sizeof(double)*mat2.ktmax),"cudaMalloc zbut");
  CHECK_DO(cudaMalloc((void**)&d_b, sizeof(double)*mat2.nd),"cudaMalloc d_b");
  CHECK_DO(cudaMalloc((void**)&d_v, sizeof(double)*mat2.nd),"cudaMallod d_v");
  //for(i=0;i<mat2.nd;i++){d_b[i]=NULL;d_v[i]=NULL;}

  len = mat2.len;
  printf("total length = %d\n", len);
  // host alloc
  // device alloc
  d_sm.nd = mat2.nd;
  d_sm.nlf = mat2.nlf;
  d_sm.ktmax = mat2.ktmax;
  cudaMalloc((void**)&d_sm.ltmtx, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.ndl, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.ndt, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.nstrtl, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.nstrtt, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.kt, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.a1, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.a2, sizeof(int)*mat2.nlf);
  cudaMalloc((void**)&d_sm.rowmat,sizeof(double)*mat2.len);
  cudaMalloc((void**)&d_sm.rowmat_t,sizeof(double)*mat2.len);
  // memcpy
  CHECK_DO(cudaMemcpy(d_sm.ltmtx, mat2.ltmtx, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.ltmtx");
  CHECK_DO(cudaMemcpy(d_sm.ndt, mat2.ndt, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.ndt");
  CHECK_DO(cudaMemcpy(d_sm.ndl, mat2.ndl, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.ndl");
  CHECK_DO(cudaMemcpy(d_sm.nstrtl, mat2.nstrtl, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.nstrtl");
  CHECK_DO(cudaMemcpy(d_sm.nstrtt, mat2.nstrtt, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.nstrtt");
  CHECK_DO(cudaMemcpy(d_sm.kt, mat2.kt, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.kt");
  CHECK_DO(cudaMemcpy(d_sm.a1, mat2.a1, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.a1");
  CHECK_DO(cudaMemcpy(d_sm.a2, mat2.a2, sizeof(int)*mat2.nlf, cudaMemcpyHostToDevice),"d_sm.a2");
  CHECK_DO(cudaMemcpy(d_sm.rowmat, mat2.rowmat, sizeof(double)*mat2.len, cudaMemcpyHostToDevice),"d_sm.rowmat");
  CHECK_DO(cudaMemcpy(d_sm.rowmat_t, mat2.rowmat_t, sizeof(double)*mat2.len, cudaMemcpyHostToDevice),"d_sm.rowmat_t");

  // 分離
  printf("begin splitting\n");
  mat2.approx = (int*)malloc(sizeof(int)*mat2.nlf);
  mat2.dense  = (int*)malloc(sizeof(int)*mat2.nlf);
  mat2.napprox = mat2.ndense = 0;
  for(ip=0; ip<mat2.nlf; ip++){
	if(mat2.ltmtx[ip]==1){
	  mat2.approx[mat2.napprox++] = ip;
	}else{
	  mat2.dense[mat2.ndense++] = ip;
	}
  }
  cudaMalloc((void**)&d_sm.approx,sizeof(int)*mat2.napprox);
  cudaMalloc((void**)&d_sm.dense,sizeof(int)*mat2.ndense);
  CHECK_DO(cudaMemcpy(d_sm.approx, mat2.approx, sizeof(int)*mat2.napprox, cudaMemcpyHostToDevice),"d_sm.approx");
  CHECK_DO(cudaMemcpy(d_sm.dense, mat2.dense, sizeof(int)*mat2.ndense, cudaMemcpyHostToDevice),"d_sm.dense");
  d_sm.napprox = mat2.napprox;
  d_sm.ndense = mat2.ndense;
  printf("end splitting\n");

#if 0
  hmvm_cudaD_kernel00dd00<16><<<mat2.ndense,32>>>						\
  (d_v, d_b, mat2.nlf, mat2.ktmax, mat2.ltmtx, mat2.ndt, mat2.ndl, mat2.nstrtl, mat2.nstrtt, mat2.kt, mat2.a1, mat2.a2, mat2.rowmat, mat2.ndense, mat2.dense);\


	hmvm_cudaD_kernel000000<<<mat2.ndense,32>>>						\
	  (d_v, d_b, mat2.nlf, mat2.ktmax, mat2.ltmtx, mat2.ndt, mat2.ndl, mat2.nstrtl, mat2.nstrtt, mat2.kt, mat2.a1, mat2.a2, mat2.rowmat, mat2.ndense, mat2.dense); \

  //FUNCNAME(d_v, d_b, d_sm, a1, a2);									\

#endif

#define BENCH(FUNCNAME,B,T,S)											\
  printf("nd = %d\n", nd);												\
  cudaMemcpy(d_v, v, sizeof(double)*nd, cudaMemcpyHostToDevice);		\
  cudaMemcpy(d_b, b, sizeof(double)*nd, cudaMemcpyHostToDevice);		\
  FUNCNAME<<<B,T,S>>> \
  (d_v, d_b, d_sm.nlf, d_sm.ktmax, d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, \
   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);					\
  cudaDeviceSynchronize();												\
  cudaMemcpy(v, d_v, sizeof(double)*nd, cudaMemcpyDeviceToHost);		\
  printf("write to %s\n", fname);												\
  F = fopen(fname, "w");												\
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);							\
  fclose(F);

#if 0
  hmvm_cudaD<<<1,1,d_sm.ktmax*sizeof(double)>>>							\
  (d_v, d_b, d_sm.nlf, d_sm.ktmax, d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, \
   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);					\

#endif

#if 0
  for(l=0;l<M+L;l++){													\
	for(i=0;i<nd;i++)v[i] = 0.0;										\
	d1 = omp_get_wtime();												\
	cudaDeviceSynchronize();											\
	d2 = omp_get_wtime();												\
	dtimes[l] = d2-d1;													\
  }																		\
  dmin = 9999.99;         dmax = 0.0;									\
  davg1 = 0.0;            davg2 = 0.0;									\
  for(i=0;i<M+L;i++)davg1 += dtimes[i];									\
  for(i=M;i<M+L;i++){													\
	if(dmin>dtimes[i])dmin=dtimes[i];									\
	if(dmax<dtimes[i])dmax=dtimes[i];									\
	davg2 += dtimes[i];													\
  }																		\
  davg1 /= (M+L);         davg2 /= L;
#endif

  if(kernel==0)
  {
	int a1, a2;
	char name[8], fname[32];
	a1=a2=0;
	snprintf(name,8,"_%d_%d",a1,a2);
	snprintf(fname,32,"hmvm_cuda1%s.txt",name);
	printf("fname = %s\n", fname);
	BENCH(hmvm_cudaD,1,1,d_sm.ktmax*sizeof(double));
	printf("TIME %d hmvm_cuda1%s min %e max %e avg1 %e avg2 %e |", M+L, name, dmin, dmax, davg1, davg2);
	for(i=0;i<M+L;i++)printf(" %e", dtimes[i]);
	printf("\n");
  }
  if(kernel==1)
  {
	int a1, a2;
	char name[8], fname[32];
	a1=a2=0;
	snprintf(name,8,"_%d_%d",a1,a2);
	snprintf(fname,32,"hmvm_cuda1blk%s.txt",name);
	printf("fname = %s\n", fname);
	BENCH(hmvm_cudaD_block,d_sm.napprox+d_sm.ndense,1,d_sm.ktmax*sizeof(double));
	printf("TIME %d hmvm_cuda1%s min %e max %e avg1 %e avg2 %e |", M+L, name, dmin, dmax, davg1, davg2);
	for(i=0;i<M+L;i++)printf(" %e", dtimes[i]);
	printf("\n");
  }

  // ######## ######## ######## ######## ######## ######## ######## ########
  cudaFree(d_sm.ltmtx);
  cudaFree(d_sm.ndl);
  cudaFree(d_sm.ndt);
  cudaFree(d_sm.nstrtl);
  cudaFree(d_sm.nstrtt);
  cudaFree(d_sm.kt);
  cudaFree(d_sm.a1);
  cudaFree(d_sm.a2);
  cudaFree(d_sm.rowmat);
  cudaFree(d_sm.rowmat_t);
  cudaFree(d_zaut);
  cudaFree(d_zbut);
  cudaFree(d_b);
  cudaFree(d_v);

  free(v); free(tmp); free(zero);
  printf("hmvm_cuda1: end\n");
}
