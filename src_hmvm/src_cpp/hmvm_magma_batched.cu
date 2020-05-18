// -*- C++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>
#include <cuda_runtime_api.h>
#include <magma_v2.h>

#include "hacapk.h"

#define CHECK_DO(act,msg) {ret=act; if(ret!=cudaSuccess){printf("%s failed\n",msg);exit(-1);};}

// ######## ######## ######## ######## ######## ######## ######## ########
/*
  magma batched simple blas
 */

void  hmvm_magma_batched_calc
(
 int n1a,
 double **d_As1a, double **d_Xs1a, double **d_Ys1a,
 int *d_Ms1a, int *d_Ns1a,
 int *d_lddas1a, int *d_incxs1a, int *d_incys1a,
 int n1b,
 double **d_As1b, double **d_Xs1b, double **d_Ys1b,
 int *d_Ms1b, int *d_Ns1b,
 int *d_lddas1b, int *d_incxs1b, int *d_incys1b,
 int n2,
 double **d_As2, double **d_Xs2, double **d_Ys2,
 int *d_Ms2, int *d_Ns2,
 int *d_lddas2, int *d_incxs2, int *d_incys2,
 magma_queue_t queue
 )
{
  double done=1.0, dzero=0.0;

  magmablas_dgemv_vbatched
	(MagmaTrans, d_Ms1a, d_Ns1a,
	 done, d_As1a, d_lddas1a,
	       d_Xs1a, d_incxs1a,
	dzero, d_Ys1a, d_incys1a, n1a, queue);

  magmablas_dgemv_vbatched_atomic
	(MagmaNoTrans, d_Ms1b, d_Ns1b,
	 done, d_As1b, d_lddas1b,
	       d_Xs1b, d_incxs1b,
	       d_Ys1b, d_incys1b, n1b, queue);
  magmablas_dgemv_vbatched_atomic
	(MagmaTrans, d_Ms2, d_Ns2,
	 done, d_As2, d_lddas2,
	       d_Xs2, d_incxs2,
	       d_Ys2, d_incys2, n2, queue);
}

void  hmvm_magma_batched_calc
(
 int n1a,
 float **d_As1a, float **d_Xs1a, float **d_Ys1a,
 int *d_Ms1a, int *d_Ns1a,
 int *d_lddas1a, int *d_incxs1a, int *d_incys1a,
 int n1b,
 float **d_As1b, float **d_Xs1b, float **d_Ys1b,
 int *d_Ms1b, int *d_Ns1b,
 int *d_lddas1b, int *d_incxs1b, int *d_incys1b,
 int n2,
 float **d_As2, float **d_Xs2, float **d_Ys2,
 int *d_Ms2, int *d_Ns2,
 int *d_lddas2, int *d_incxs2, int *d_incys2,
 magma_queue_t queue
 )
{
  float done=1.0f, dzero=0.0f;

  magmablas_sgemv_vbatched
	(MagmaTrans, d_Ms1a, d_Ns1a,
	 done, d_As1a, d_lddas1a,
	       d_Xs1a, d_incxs1a,
	dzero, d_Ys1a, d_incys1a, n1a, queue);

  magmablas_sgemv_vbatched_atomic
	(MagmaNoTrans, d_Ms1b, d_Ns1b,
	 done, d_As1b, d_lddas1b,
	       d_Xs1b, d_incxs1b,
	       d_Ys1b, d_incys1b, n1b, queue);
  magmablas_sgemv_vbatched_atomic
	(MagmaTrans, d_Ms2, d_Ns2,
	 done, d_As2, d_lddas2,
	       d_Xs2, d_incxs2,
	       d_Ys2, d_incys2, n2, queue);
}

void  hmvm_magma_batched_calc_nocheck
(
 int n1a,
 double **d_As1a, double **d_Xs1a, double **d_Ys1a,
 int *d_Ms1a, int *d_Ns1a,
 int *d_lddas1a, int *d_incxs1a, int *d_incys1a,
 int n1b,
 double **d_As1b, double **d_Xs1b, double **d_Ys1b,
 int *d_Ms1b, int *d_Ns1b,
 int *d_lddas1b, int *d_incxs1b, int *d_incys1b,
 int n2,
 double **d_As2, double **d_Xs2, double **d_Ys2,
 int *d_Ms2, int *d_Ns2,
 int *d_lddas2, int *d_incxs2, int *d_incys2,
 magma_queue_t queue
 )
{
  double done=1.0, dzero=0.0;

  magmablas_dgemv_vbatched_nocheck
	(MagmaTrans, d_Ms1a, d_Ns1a,
	 done, d_As1a, d_lddas1a,
	       d_Xs1a, d_incxs1a,
	dzero, d_Ys1a, d_incys1a, n1a, queue);

  magmablas_dgemv_vbatched_nocheck_atomic
	(MagmaNoTrans, d_Ms1b, d_Ns1b,
	 done, d_As1b, d_lddas1b,
	       d_Xs1b, d_incxs1b,
	       d_Ys1b, d_incys1b, n1b, queue);
  magmablas_dgemv_vbatched_nocheck_atomic
	(MagmaTrans, d_Ms2, d_Ns2,
	 done, d_As2, d_lddas2,
	       d_Xs2, d_incxs2,
	       d_Ys2, d_incys2, n2, queue);
}

void  hmvm_magma_batched_calc_nocheck
(
 int n1a,
 float **d_As1a, float **d_Xs1a, float **d_Ys1a,
 int *d_Ms1a, int *d_Ns1a,
 int *d_lddas1a, int *d_incxs1a, int *d_incys1a,
 int n1b,
 float **d_As1b, float **d_Xs1b, float **d_Ys1b,
 int *d_Ms1b, int *d_Ns1b,
 int *d_lddas1b, int *d_incxs1b, int *d_incys1b,
 int n2,
 float **d_As2, float **d_Xs2, float **d_Ys2,
 int *d_Ms2, int *d_Ns2,
 int *d_lddas2, int *d_incxs2, int *d_incys2,
 magma_queue_t queue
 )
{
  float done=1.0f, dzero=0.0f;

  magmablas_sgemv_vbatched_nocheck
	(MagmaTrans, d_Ms1a, d_Ns1a,
	 done, d_As1a, d_lddas1a,
	       d_Xs1a, d_incxs1a,
	dzero, d_Ys1a, d_incys1a, n1a, queue);

  magmablas_sgemv_vbatched_nocheck_atomic
	(MagmaNoTrans, d_Ms1b, d_Ns1b,
	 done, d_As1b, d_lddas1b,
	       d_Xs1b, d_incxs1b,
	       d_Ys1b, d_incys1b, n1b, queue);
  magmablas_sgemv_vbatched_nocheck_atomic
	(MagmaTrans, d_Ms2, d_Ns2,
	 done, d_As2, d_lddas2,
	       d_Xs2, d_incxs2,
	       d_Ys2, d_incys2, n2, queue);
}

template <class T>
void hmvm_magma_batched_proxy
(
 int nbatch1a,
 T **d_As1a, T **d_Xs1a, T **d_Ys1a,
 int *d_Ms1a, int *d_Ns1a,
 int *d_lddas1a, int *d_incxs1a, int *d_incys1a,
 int nbatch1b,
 T **d_As1b, T **d_Xs1b, T **d_Ys1b,
 int *d_Ms1b, int *d_Ns1b,
 int *d_lddas1b, int *d_incxs1b, int *d_incys1b,
 int nbatch2,
 T **d_As2, T **d_Xs2, T **d_Ys2,
 int *d_Ms2, int *d_Ns2,
 int *d_lddas2, int *d_incxs2, int *d_incys2,
 magma_queue_t queue, int opt,
 T *d_v, int nd, char *fname, int bench
){
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  T *d_tmp, *h_zero;
  T *v;
  v = new T[nd];
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	//CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();

	if(opt==0){
	  hmvm_magma_batched_calc
		(
		 nbatch1a, d_As1a, d_Xs1a, d_Ys1a, d_Ms1a, d_Ns1a, d_lddas1a, d_incxs1a, d_incys1a,
		 nbatch1b, d_As1b, d_Xs1b, d_Ys1b, d_Ms1b, d_Ns1b, d_lddas1b, d_incxs1b, d_incys1b,
		 nbatch2, d_As2, d_Xs2, d_Ys2, d_Ms2, d_Ns2, d_lddas2, d_incxs2, d_incys2,
		 queue
		 );
	}else{
	  hmvm_magma_batched_calc_nocheck
		(
		 nbatch1a, d_As1a, d_Xs1a, d_Ys1a, d_Ms1a, d_Ns1a, d_lddas1a, d_incxs1a, d_incys1a,
		 nbatch1b, d_As1b, d_Xs1b, d_Ys1b, d_Ms1b, d_Ns1b, d_lddas1b, d_incxs1b, d_incys1b,
		 nbatch2, d_As2, d_Xs2, d_Ys2, d_Ms2, d_Ns2, d_lddas2, d_incxs2, d_incys2,
		 queue
		 );
	}

	cudaDeviceSynchronize();
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }
  if(bench==0){
	CHECK_DO(cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v");
	printf("write to %s\n", fname);
	F = fopen(fname, "w");
	for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	fclose(F);
  }else{
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=M;i<L;i++){
	  davg += dtimes[i];
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	}
	davg /= (L-M);
	printf("TIME %d hmvm_magma_batched%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}

// ######## ######## ######## ######## ######## ######## ######## ########

// ######## ######## ######## ######## ######## ######## ######## ########

template<class T>
void hmvm_magma_batched(matrix2<T> *mat2, T *b, int kernel, int dump_result)
{
  matrix2<T> d_sm;
  int i, l, nd = mat2->nd, ktmax = mat2->ktmax, nlf = mat2->nlf;
  T *v=NULL;
  T *d_b, *d_v;
  int ip;
  int len;
  cudaError_t ret;
  printf("hmvm_magma_batched_%s: begin\n", typeid(T).name()); fflush(stdout);
  v    = new T[nd];
  for(i=0;i<nd;i++){
	v[i] = (T)0.0;
  }
  CHECK_DO(cudaMalloc((void**)&d_b, sizeof(T)*nd),"cudaMalloc d_b");
  CHECK_DO(cudaMalloc((void**)&d_v, sizeof(T)*nd),"cudaMalloc d_v");
  CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");

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

#if 0
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
#endif

  int n, kt, ndl, ndt, nstrtl, nstrtt, ltmtx, head;
  magma_device_t dev;
  magma_queue_t queue;
  magma_queue_create(dev, &queue);
  int nbatch1a, nbatch1b, nbatch2;
  T **tmpvec;
  tmpvec = (T**)malloc(sizeof(T*)*nlf);
  for(ip=0;ip<nlf;ip++){
	CHECK_DO(cudaMalloc((void**)&tmpvec[ip],sizeof(T)*mat2->kt[ip]),"cudaMalloc ");
	//cudaMemcpy(tmpvec[ip], tmpzero, sizeof(T)*mat2->kt[ip], cudaMemcpyHostToDevice);
  }

  // host
  T **h_As1a, **h_Xs1a, **h_Ys1a;
  int *h_Ms1a, *h_Ns1a;
  int *h_lddas1a, *h_incxs1a, *h_incys1a;
  T **h_As1b, **h_Xs1b, **h_Ys1b;
  int *h_Ms1b, *h_Ns1b;
  int *h_lddas1b, *h_incxs1b, *h_incys1b;
  T **h_As2, **h_Xs2, **h_Ys2;
  int *h_Ms2, *h_Ns2;
  int *h_lddas2, *h_incxs2, *h_incys2;
  // device
  T **d_As1a, **d_Xs1a, **d_Ys1a;
  int *d_Ms1a, *d_Ns1a;
  int *d_lddas1a, *d_incxs1a, *d_incys1a;
  T **d_As1b, **d_Xs1b, **d_Ys1b;
  int *d_Ms1b, *d_Ns1b;
  int *d_lddas1b, *d_incxs1b, *d_incys1b;
  T **d_As2, **d_Xs2, **d_Ys2;
  int *d_Ms2, *d_Ns2;
  int *d_lddas2, *d_incxs2, *d_incys2;
  // allocation
  // 1a approx. 1
  h_As1a = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Xs1a = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Ys1a = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Ms1a = (int*)malloc(sizeof(int)*(nlf+1));
  h_Ns1a = (int*)malloc(sizeof(int)*(nlf+1));
  h_lddas1a = (int*)malloc(sizeof(int)*(nlf+1));
  h_incxs1a = (int*)malloc(sizeof(int)*(nlf+1));
  h_incys1a = (int*)malloc(sizeof(int)*(nlf+1));
  CHECK_DO(cudaMalloc((void**)&d_As1a, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Xs1a, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ys1a, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ms1a, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ns1a, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_lddas1a, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_incxs1a, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_incys1a, sizeof(int)*(nlf+1)),"cudaMalloc ");
  // 1b. approx. 2
  h_As1b = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Xs1b = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Ys1b = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Ms1b = (int*)malloc(sizeof(int)*(nlf+1));
  h_Ns1b = (int*)malloc(sizeof(int)*(nlf+1));
  h_lddas1b = (int*)malloc(sizeof(int)*(nlf+1));
  h_incxs1b = (int*)malloc(sizeof(int)*(nlf+1));
  h_incys1b = (int*)malloc(sizeof(int)*(nlf+1));
  CHECK_DO(cudaMalloc((void**)&d_As1b, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Xs1b, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ys1b, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ms1b, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ns1b, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_lddas1b, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_incxs1b, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_incys1b, sizeof(int)*(nlf+1)),"cudaMalloc ");
  // 2. dense
  h_As2 = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Xs2 = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Ys2 = (T**)malloc(sizeof(T*)*(nlf+1));
  h_Ms2 = (int*)malloc(sizeof(int)*(nlf+1));
  h_Ns2 = (int*)malloc(sizeof(int)*(nlf+1));
  h_lddas2 = (int*)malloc(sizeof(int)*(nlf+1));
  h_incxs2 = (int*)malloc(sizeof(int)*(nlf+1));
  h_incys2 = (int*)malloc(sizeof(int)*(nlf+1));
  CHECK_DO(cudaMalloc((void**)&d_As2, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Xs2, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ys2, sizeof(T*)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ms2, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_Ns2, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_lddas2, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_incxs2, sizeof(int)*(nlf+1)),"cudaMalloc ");
  CHECK_DO(cudaMalloc((void**)&d_incys2, sizeof(int)*(nlf+1)),"cudaMalloc ");
  // construction
  // 1a. approx. 1
  n=0;
  for(ip=0; ip<nlf; ip++){
	ndl    = mat2->ndl[ip];
	ndt    = mat2->ndt[ip];
	nstrtl = mat2->nstrtl[ip];
	nstrtt = mat2->nstrtt[ip];
	ltmtx  = mat2->ltmtx[ip];
	if(ltmtx==1){
	  kt     = mat2->kt[ip];
	  head = mat2->a1[ip];
	  h_As1a[n] = &d_sm.rowmat[head];
	  h_Xs1a[n] = &d_b[nstrtt-1];
	  h_Ys1a[n] = tmpvec[ip];
	  h_Ms1a[n] = ndt;
	  h_Ns1a[n] = kt;
	  h_lddas1a[n] = ndt;
	  h_incxs1a[n] = 1;
	  h_incys1a[n] = 1;
	  n++;
	  //cublasDgemv(handle,CUBLAS_OP_T, ndt,kt, &done, &d_sm.rowmat[head], ndt,&d_zu[nstrtt-1],1,&dzero,d_zbut,1);
	} else if(ltmtx==2){
	}
  }
  nbatch1a = n;
  cudaMemcpy(d_As1a, h_As1a, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Xs1a, h_Xs1a, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ys1a, h_Ys1a, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ms1a, h_Ms1a, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ns1a, h_Ns1a, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lddas1a, h_lddas1a, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_incxs1a, h_incxs1a, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_incys1a, h_incys1a, sizeof(int)*n, cudaMemcpyHostToDevice);
  // 1b. approx. 2
  n = 0;
  for(ip=0; ip<nlf; ip++){
	ndl    = mat2->ndl[ip];
	ndt    = mat2->ndt[ip];
	nstrtl = mat2->nstrtl[ip];
	nstrtt = mat2->nstrtt[ip];
	ltmtx  = mat2->ltmtx[ip];
	if(ltmtx==1){
	  kt = mat2->kt[ip];
	  head = mat2->a2[ip];
	  h_As1b[n] = &d_sm.rowmat[head];
	  h_Xs1b[n] = tmpvec[ip];
	  h_Ys1b[n] = &d_v[nstrtl-1];
	  h_Ms1b[n] = ndl;
	  h_Ns1b[n] = kt;
	  h_lddas1b[n] = ndl;
	  h_incxs1b[n] = 1;
	  h_incys1b[n] = 1;
	  n++;
	  //cublasDgemv(handle,CUBLAS_OP_N, ndl,kt, &done, &d_sm.rowmat[head], ndl,d_zbut,1,&done,&d_zaut[nstrtl-1],1);
	} else if(ltmtx==2){
	}
  }
  nbatch1b = n;
  cudaMemcpy(d_As1b, h_As1b, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Xs1b, h_Xs1b, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ys1b, h_Ys1b, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ms1b, h_Ms1b, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ns1b, h_Ns1b, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lddas1b, h_lddas1b, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_incxs1b, h_incxs1b, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_incys1b, h_incys1b, sizeof(int)*n, cudaMemcpyHostToDevice);
  // 2. dense
  n = 0;
  for(ip=0; ip<nlf; ip++){
	ndl    = mat2->ndl[ip];
	ndt    = mat2->ndt[ip];
	nstrtl = mat2->nstrtl[ip];
	nstrtt = mat2->nstrtt[ip];
	ltmtx  = mat2->ltmtx[ip];
	if(ltmtx==1){
	} else if(ltmtx==2){
	  head = mat2->a1[ip];
	  h_As2[n] = &d_sm.rowmat[head];
	  h_Xs2[n] = &d_b[nstrtt-1];
	  h_Ys2[n] = &d_v[nstrtl-1];
	  h_Ms2[n] = ndt;
	  h_Ns2[n] = ndl;
	  h_lddas2[n] = ndt;
	  h_incxs2[n] = 1;
	  h_incys2[n] = 1;
	  n++;
	  //cublasDgemv(handle,CUBLAS_OP_T, ndt,ndl, &done, &d_sm.rowmat[head], ndt,&d_zu[nstrtt-1],1,&done,&d_zaut[nstrtl-1],1);
	}
  }
  nbatch2 = n;
  cudaMemcpy(d_As2, h_As2, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Xs2, h_Xs2, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ys2, h_Ys2, sizeof(T*)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ms2, h_Ms2, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ns2, h_Ns2, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lddas2, h_lddas2, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_incxs2, h_incxs2, sizeof(int)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_incys2, h_incys2, sizeof(int)*n, cudaMemcpyHostToDevice);

  /*
	magma batched blas
	0: default
	1: nocheck
  */
  {
	char *names[]={"default","nocheck"};
	char name[0xff], fname[0xff];
	snprintf(name,0xff,"magma_batched_%d_%s_%s", kernel, names[kernel], typeid(T).name());
	snprintf(fname,0xff,"result_%s.txt", name);
	printf("fname = %s\n", fname);

	// EXEC
	hmvm_magma_batched_proxy<T>
	  (
	   nbatch1a, d_As1a, d_Xs1a, d_Ys1a, d_Ms1a, d_Ns1a, d_lddas1a, d_incxs1a, d_incys1a,
	   nbatch1b, d_As1b, d_Xs1b, d_Ys1b, d_Ms1b, d_Ns1b, d_lddas1b, d_incxs1b, d_incys1b,
	   nbatch2, d_As2, d_Xs2, d_Ys2, d_Ms2, d_Ns2, d_lddas2, d_incxs2, d_incys2,
	   queue, kernel,
	   d_v, nd, fname, 0);
	// BENCH
	hmvm_magma_batched_proxy<T>
	  (
	   nbatch1a, d_As1a, d_Xs1a, d_Ys1a, d_Ms1a, d_Ns1a, d_lddas1a, d_incxs1a, d_incys1a,
	   nbatch1b, d_As1b, d_Xs1b, d_Ys1b, d_Ms1b, d_Ns1b, d_lddas1b, d_incxs1b, d_incys1b,
	   nbatch2, d_As2, d_Xs2, d_Ys2, d_Ms2, d_Ns2, d_lddas2, d_incxs2, d_incys2,
	   queue, kernel,
	   d_v, nd, fname, 5);
  }
  magma_queue_destroy(queue);

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
  cudaFree(d_b);
  cudaFree(d_v);

  delete [] v;
  printf("hmvm_magma_batched: end\n");
}


// ######## ######## ######## ######## ######## ######## ######## ########
// template関数の実体化のための宣言
// ######## ######## ######## ######## ######## ######## ######## ########
template void hmvm_magma_batched<float>(matrix2<float>  *mat2, float *b, int kernel, int dump_result);
template void hmvm_magma_batched<double>(matrix2<double> *mat2, double *b, int kernel, int dump_result);
