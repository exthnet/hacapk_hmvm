// -*- C++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>
#include <cuda_runtime_api.h>

#include "hacapk.h"
#include "hmvm_cuda_kernels.h"

#define CHECK_DO(act,msg) {ret=act; if(ret!=cudaSuccess){printf("%s failed\n",msg);exit(-1);};}


template<class T>
void hmvm_cuda1(matrix2<T> *mat2, T *b, int kernel, int dump_result)
{
  const int L=10, M=5;
  FILE *F;
  matrix2<T> d_sm;
  int i, l, nd = mat2->nd, ktmax = mat2->ktmax, nlf = mat2->nlf;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  T *v=NULL, *tmp=NULL, *zero;
  T *d_b, *d_v;//, *d_zaut, *d_zbut;
  int ip;
  int len;
  cudaError_t ret;
  printf("hmvm_cuda1_%s: begin\n", typeid(T).name());
  v    = new T[nd];    //(double*)malloc(sizeof(double)*mat2->nd);
  tmp  = new T[ktmax]; //(double*)malloc(sizeof(double)*mat2->ktmax);
  zero = new T[ktmax]; //(double*)malloc(sizeof(double)*mat2->ktmax);
  for(i=0;i<nd;i++){
	v[i] = (T)0.0;
  }
  for(i=0;i<ktmax;i++){
	zero[i] = (T)0.0;
  }
  //CHECK_DO(cudaMalloc((void**)&d_zaut, sizeof(T)*mat2->nd),"cudaMalloc z_aut");
  //CHECK_DO(cudaMalloc((void**)&d_zbut, sizeof(T)*mat2->ktmax),"cudaMalloc zbut");
  CHECK_DO(cudaMalloc((void**)&d_b, sizeof(T)*nd),"cudaMalloc d_b");
  CHECK_DO(cudaMalloc((void**)&d_v, sizeof(T)*nd),"cudaMalloc d_v");
  //for(i=0;i<mat2->nd;i++){d_b[i]=NULL;d_v[i]=NULL;}

  printf("nd = %d\n", nd);												\
  len = mat2->len;
  printf("total length = %d\n", len);
  // host alloc
  // device alloc
  d_sm.nd    = nd;
  d_sm.nlf   = nlf;
  d_sm.ktmax = ktmax;
  cudaMalloc((void**)&d_sm.ltmtx, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.ndl, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.ndt, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.nstrtl, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.nstrtt, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.kt, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.a1, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.a2, sizeof(int)*nlf);
  cudaMalloc((void**)&d_sm.rowmat,sizeof(T)*mat2->len);
  cudaMalloc((void**)&d_sm.rowmat_t,sizeof(T)*mat2->len);
  // memcpy
  CHECK_DO(cudaMemcpy(d_sm.ltmtx, mat2->ltmtx, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.ltmtx");
  CHECK_DO(cudaMemcpy(d_sm.ndt, mat2->ndt, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.ndt");
  CHECK_DO(cudaMemcpy(d_sm.ndl, mat2->ndl, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.ndl");
  CHECK_DO(cudaMemcpy(d_sm.nstrtl, mat2->nstrtl, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.nstrtl");
  CHECK_DO(cudaMemcpy(d_sm.nstrtt, mat2->nstrtt, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.nstrtt");
  CHECK_DO(cudaMemcpy(d_sm.kt, mat2->kt, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.kt");
  CHECK_DO(cudaMemcpy(d_sm.a1, mat2->a1, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.a1");
  CHECK_DO(cudaMemcpy(d_sm.a2, mat2->a2, sizeof(int)*nlf, cudaMemcpyHostToDevice),"d_sm.a2");
  CHECK_DO(cudaMemcpy(d_sm.rowmat, mat2->rowmat, sizeof(T)*mat2->len, cudaMemcpyHostToDevice),"d_sm.rowmat");
  CHECK_DO(cudaMemcpy(d_sm.rowmat_t, mat2->rowmat_t, sizeof(T)*mat2->len, cudaMemcpyHostToDevice),"d_sm.rowmat_t");

  // 分離
  printf("begin splitting\n");
  mat2->approx = (int*)malloc(sizeof(int)*nlf);
  mat2->dense  = (int*)malloc(sizeof(int)*nlf);
  mat2->napprox = mat2->ndense = 0;
  for(ip=0; ip<nlf; ip++){
	if(mat2->ltmtx[ip]==1){
	  mat2->approx[mat2->napprox++] = ip;
	}else{
	  mat2->dense[mat2->ndense++] = ip;
	}
  }
  cudaMalloc((void**)&d_sm.approx,sizeof(int)*mat2->napprox);
  cudaMalloc((void**)&d_sm.dense,sizeof(int)*mat2->ndense);
  CHECK_DO(cudaMemcpy(d_sm.approx, mat2->approx, sizeof(int)*mat2->napprox, cudaMemcpyHostToDevice),"d_sm.approx");
  CHECK_DO(cudaMemcpy(d_sm.dense, mat2->dense, sizeof(int)*mat2->ndense, cudaMemcpyHostToDevice),"d_sm.dense");
  d_sm.napprox = mat2->napprox;
  d_sm.ndense  = mat2->ndense;
  printf("end splitting (napprox=%d, ndense=%d)\n", mat2->napprox, mat2->ndense);

#define EXEC(FUNCNAME,BLOCKS,THREADS,S)									\
  for(i=0;i<nd;i++)v[i] = (T)0.0;										\
  CHECK_DO(cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v"); \
  CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b"); \
  cudaDeviceSynchronize();												\
  FUNCNAME<<<BLOCKS,THREADS,S>>>										\
	(d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,			\
	 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,							\
	 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
  cudaDeviceSynchronize();												\
  CHECK_DO(cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v"); \
  printf("write to %s\n", fname);										\
  F = fopen(fname, "w");												\
  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);							\
  fclose(F);

#define BENCH(FUNCNAME,BLOCKS,THREADS,S)
#if 0
#define BENCH(FUNCNAME,BLOCKS,THREADS,S)								\
  for(l=0;l<L;l++){														\
	for(i=0;i<nd;i++)v[i] = (T)0.0;										\
	CHECK_DO(cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v"); \
	CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b"); \
	cudaDeviceSynchronize();											\
	d1 = omp_get_wtime();												\
	FUNCNAME<<<BLOCKS,THREADS,S>>>										\
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,		\
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,							\
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
	cudaDeviceSynchronize();											\
	d2 = omp_get_wtime();												\
	dtimes[l] = d2-d1;													\
  }																		\
  dmin = 9999.99;														\
  dmax = 0.0;															\
  davg = 0.0;															\
  for(i=M;i<L;i++){														\
	davg += dtimes[i];													\
	if(dmin>dtimes[i])dmin=dtimes[i];									\
	if(dmax<dtimes[i])dmax=dtimes[i];									\
  }																		\
  davg /= (L-M);
#endif

#define EXEC2(FUNCNAME,BLOCKS,THREADS,S,DIV,OPT)							\
  for(i=0;i<nd;i++)v[i] = (T)0.0;										\
  CHECK_DO(cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v"); \
  CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b"); \
  cudaDeviceSynchronize();												\
  FUNCNAME<T,DIV,OPT><<<BLOCKS,THREADS,S>>>							\
	(d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,			\
	 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,							\
	 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
  cudaDeviceSynchronize();												\
  CHECK_DO(cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v"); \
  printf("write to %s\n", fname);										\
  F = fopen(fname, "w");												\
  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);							\
  fclose(F);

#define BENCH2(FUNCNAME,BLOCKS,THREADS,S,DIV,OPT)
#if 0
#define BENCH2(FUNCNAME,BLOCKS,THREADS,S,DIV,OPT)							\
  for(l=0;l<L;l++){														\
	for(i=0;i<nd;i++)v[i] = (T)0.0;										\
	CHECK_DO(cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v"); \
	CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b"); \
	cudaDeviceSynchronize();											\
	d1 = omp_get_wtime();												\
	FUNCNAME<T,DIV,OPT><<<BLOCKS,THREADS,S>>>						\
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,		\
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,							\
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
	cudaDeviceSynchronize();											\
	d2 = omp_get_wtime();												\
	dtimes[l] = d2-d1;													\
  }																		\
  dmin = 9999.99;														\
  dmax = 0.0;															\
  davg = 0.0;															\
  for(i=M;i<L;i++){														\
	davg += dtimes[i];													\
	if(dmin>dtimes[i])dmin=dtimes[i];									\
	if(dmax<dtimes[i])dmax=dtimes[i];									\
  }																		\
  davg /= (L-M);
#endif

  // sequential
  if(kernel==0)
  {
	int a1, a2;
	char name[8], fname[0xff];
	a1=a2=0;
	snprintf(name,8,"_%d_%d",a1,a2);
	snprintf(fname,32,"result_cuda1_seq%s_%s.txt", name, typeid(T).name());
	printf("fname = %s\n", fname);
	EXEC(hmvm_cuda_seq,1,1,d_sm.ktmax*sizeof(T));
	BENCH(hmvm_cuda_seq,1,1,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_seq%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }

  // block parallel
  // whole hmvm calculation in 1 GPU kernel
  if(kernel==1)
  {
	int a1, a2;
	char name[8], fname[32];
	a1=a2=0;
	snprintf(name,8,"_%d_%d",a1,a2);
	snprintf(fname,32,"result_cuda1_block%s_%s.txt", name, typeid(T).name());
	printf("fname = %s\n", fname);
	EXEC(hmvm_cuda_block,d_sm.napprox+d_sm.ndense,1,d_sm.ktmax*sizeof(T));
	BENCH(hmvm_cuda_block,d_sm.napprox+d_sm.ndense,1,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_block%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }

  // block & thread parallel
  // whole hmvm calculation in 1 GPU kernel
  // under development
  if(kernel==10)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div1_shuffl_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,1,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,1,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_%s min %e max %e avg %e\n", L-M, name, dmin, dmax, davg);
  }
  if(kernel==11)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div2_shuffl_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,2,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,2,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_%s min %e max %e avg %e\n", L-M, name, dmin, dmax, davg);
  }
  if(kernel==12)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div4_shuffl_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,4,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,4,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_%s min %e max %e avg %e\n", L-M, name, dmin, dmax, davg);
  }
  if(kernel==13)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div8_shuffl_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,8,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,8,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_%s min %e max %e avg %e\n", L-M, name, dmin, dmax, davg);
  }
  if(kernel==14)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div16_shuffl_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,16,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,16,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_%s min %e max %e avg %e\n", L-M, name, dmin, dmax, davg);
  }
  if(kernel==15)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div32_shuffl_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,32,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,32,0>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_%s min %e max %e avg %e\n", L-M, name, dmin, dmax, davg);
  }

  if(kernel==20)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div1_atomic_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,1,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,1,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_hybrid1b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  if(kernel==21)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div2_atomic_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,2,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,2,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_hybrid1b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  if(kernel==22)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div4_atomic_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,4,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,4,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_hybrid1b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  if(kernel==23)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div8_atomic_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,8,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,8,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_hybrid1b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  if(kernel==24)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div16_atomic_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,16,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,16,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_hybrid1b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  if(kernel==25)
  {
	char name[32], fname[64];
	snprintf(name,32,"hybrid1_div32_atomic_%d_%s", kernel, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	EXEC((hmvm_cuda_hybrid1<T,32,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	BENCH((hmvm_cuda_hybrid1<T,32,1>),d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1_hybrid1b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }

#if 1
  if(kernel>=1000&&kernel<1192){
	int subkernel = kernel-1000;
	int DIV, MUL, ATOMIC;
	ATOMIC = subkernel/96;
	MUL = (subkernel%96)/6 + 1;
	DIV = pow(2,subkernel%6);
	char name[32], fname[64];
	snprintf(name,32,"hybrid2_div%d_mul%d_atomic%d_%s", DIV, MUL, ATOMIC, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	printf("DIV = %d, MUL = %d, ATOMIC = %d\n", DIV, MUL, ATOMIC);
	// EXEC
	{
	  hmvm_cuda_hybrid2_proxy<T>(d_v, d_b, d_sm.nlf, d_sm.ktmax,
								 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
								 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
								 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
								 d_sm.napprox+d_sm.ndense,32*MUL,d_sm.ktmax*sizeof(T),DIV,MUL,ATOMIC, v, b, nd, fname, 0);
	}
	// BENCH
	{
	  hmvm_cuda_hybrid2_proxy<T>(d_v, d_b, d_sm.nlf, d_sm.ktmax,
								 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
								 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
								 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
								 d_sm.napprox+d_sm.ndense,32*MUL,d_sm.ktmax*sizeof(T),DIV,MUL,ATOMIC, v, b, nd, fname, 1);
	}
  }
#endif

#if 1
  if(kernel>=2000&&kernel<2192){
	int subkernel = kernel-2000;
	int DIV, MUL, ATOMIC;
	ATOMIC = subkernel/96;
	MUL = (subkernel%96)/6 + 1;
	DIV = pow(2,subkernel%6);
	ATOMIC = 1; MUL = 1; DIV = 1; // test
	char name[32], fname[64];
	snprintf(name,32,"hybrid3_div%d_mul%d_atomic%d_%s", DIV, MUL, ATOMIC, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	printf("DIV = %d, MUL = %d, ATOMIC = %d\n", DIV, MUL, ATOMIC);
	// EXEC
	{
	  hmvm_cuda_hybrid3_proxy<T>(d_v, d_b, d_sm.nlf, d_sm.ktmax,
								 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
								 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
								 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
								 (d_sm.napprox+MUL-1)/MUL+(d_sm.ndense+MUL-1)/MUL,32*MUL,d_sm.ktmax*sizeof(T),DIV,MUL,ATOMIC, v, b, nd, fname, 0);
	}
	// BENCH
	if(0){
	  hmvm_cuda_hybrid3_proxy<T>(d_v, d_b, d_sm.nlf, d_sm.ktmax,
								 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
								 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
								 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
								 (d_sm.napprox+MUL-1)/MUL+(d_sm.ndense+MUL-1)/MUL,32*MUL,d_sm.ktmax*sizeof(T),DIV,MUL,ATOMIC, v, b, nd, fname, 1);
	}
  }
#endif
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
  //cudaFree(d_zaut);
  //cudaFree(d_zbut);
  cudaFree(d_b);
  cudaFree(d_v);

  delete [] v; delete [] tmp; delete [] zero;
  //free(v); free(tmp); free(zero);
  printf("hmvm_cuda1: end\n");
}


// ######## ######## ######## ######## ######## ######## ######## ########
// template関数の実体化のための宣言
// ######## ######## ######## ######## ######## ######## ######## ########
template void hmvm_cuda1<float>(matrix2<float>  *mat2, float *b, int kernel, int dump_result);
template void hmvm_cuda1<double>(matrix2<double> *mat2, double *b, int kernel, int dump_result);
