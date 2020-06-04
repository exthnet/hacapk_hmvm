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
  magma simple blas
 */

void  hmvm_magma_calc
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, double *rowmat, double *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 double *d_zbut, matrix2<double> *mat)
{
  int ip;
  int ndl,ndt,nstrtl,nstrtt,kt,ltmtx;
  double dzero = 0.0;
  double done = 1.0;
  int head;
  magma_device_t dev;
  magma_queue_t queue;
  magma_queue_create(dev, &queue);

  for(ip=0; ip<nlf; ip++){
    ndl    = _ndl[ip];
    ndt    = _ndt[ip];
    nstrtl = _nstrtl[ip];
    nstrtt = _nstrtt[ip];
    ltmtx  = _ltmtx[ip];
    if(ltmtx==1){
      kt = _kt[ip];
	  //cudaMemcpy(&d_zbut, &h_zero, sizeof(double)*kt, cudaMemcpyHostToDevice);
	  head = a1[ip];
	  magma_dgemv(MagmaTrans, ndt,kt, done, &rowmat[head], ndt,&d_zu[nstrtt-1],1,dzero,d_zbut,1, queue);
	  head = a2[ip];
	  magma_dgemv(MagmaNoTrans, ndl,kt, done, &rowmat[head], ndl,d_zbut,1,done,&d_zaut[nstrtl-1],1, queue);
    } else if(ltmtx==2){
	  head = a1[ip];
	  magma_dgemv(MagmaTrans, ndt,ndl, done, &rowmat[head], ndt,&d_zu[nstrtt-1],1,done,&d_zaut[nstrtl-1],1, queue);
    }
  }
  //magma_daxpy(nd,done,d_zaut,1,d_zau,1, queue);

  magma_queue_destroy(queue);
}

void  hmvm_magma_calc
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, float *rowmat, float *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 float *d_zbut, matrix2<float> *mat)
{
  int ip;
  int ndl,ndt,nstrtl,nstrtt,kt,ltmtx;
  double dzero = 0.0;
  double done = 1.0;
  int head;
  magma_device_t dev;
  magma_queue_t queue;
  magma_queue_create(dev, &queue);

  for(ip=0; ip<nlf; ip++){
    ndl    = _ndl[ip];
    ndt    = _ndt[ip];
    nstrtl = _nstrtl[ip];
    nstrtt = _nstrtt[ip];
    ltmtx  = _ltmtx[ip];
    if(ltmtx==1){
      kt = _kt[ip];
	  //cudaMemcpy(&d_zbut, &h_zero, sizeof(double)*kt, cudaMemcpyHostToDevice);
	  head = a1[ip];
	  magma_sgemv(MagmaTrans, ndt,kt, done, &rowmat[head], ndt,&d_zu[nstrtt-1],1,dzero,d_zbut,1, queue);
	  head = a2[ip];
	  magma_sgemv(MagmaNoTrans, ndl,kt, done, &rowmat[head], ndl,d_zbut,1,done,&d_zaut[nstrtl-1],1, queue);
    } else if(ltmtx==2){
	  head = a1[ip];
	  magma_sgemv(MagmaTrans, ndt,ndl, done, &rowmat[head], ndt,&d_zu[nstrtt-1],1,done,&d_zaut[nstrtl-1],1, queue);
    }
  }
  //magma_daxpy(nd,done,d_zaut,1,d_zau,1, queue);

  magma_queue_destroy(queue);
}

template <class T>
void hmvm_magma_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 T *v, T *b, int nd, char *fname, int bench,
 matrix2<T> *mat2)
{
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  T *d_tmp, *h_zero;
  cudaMalloc((void**)&d_tmp, sizeof(T)*mat2->ktmax);
  h_zero = new T[mat2->ktmax];
  for(i=0; i<mat2->ktmax; i++)h_zero[i]=(T)0.0;
  cudaMemcpy(d_tmp, h_zero, sizeof(T)*mat2->ktmax, cudaMemcpyHostToDevice);
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	//CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();

	hmvm_magma_calc
	  (d_zaut, d_zu, nlf, ktmax,
	   ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense, d_tmp, mat2);

	cudaDeviceSynchronize();
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }
  if(bench==0){
	CHECK_DO(cudaMemcpy(v, d_zaut, sizeof(T)*nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v");
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
	printf("TIME %d hmvm_magma_%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
}

// ######## ######## ######## ######## ######## ######## ######## ########

// ######## ######## ######## ######## ######## ######## ######## ########

template<class T>
void hmvm_magma(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench)
{
  matrix2<T> d_sm;
  int i, nd = mat2->nd, ktmax = mat2->ktmax, nlf = mat2->nlf;
  T *v=NULL;
  T *d_b, *d_v;
  int len;
  cudaError_t ret;
  printf("hmvm_magma_%s: begin\n", typeid(T).name()); fflush(stdout);
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

  /*
	magma blas
  */
  {
	char name[0xff], fname[0xff];
	snprintf(name,0xff,"magma_%s", typeid(T).name());
	snprintf(fname,0xff,"result_%s.txt", name);
	printf("fname = %s\n", fname);
	// EXEC
	if(dump_result)
	  hmvm_magma_proxy<T>
		(d_v, d_b, mat2->nlf, mat2->ktmax,
		 mat2->ltmtx, mat2->ndt, mat2->ndl, mat2->nstrtl, mat2->nstrtt, mat2->kt,
		 mat2->a1, mat2->a2, d_sm.rowmat, d_sm.rowmat_t,
		 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
		 v, b, nd, fname, 0,
		 mat2);
	// BENCH
	if(nbench>0)
	  hmvm_magma_proxy<T>
		(d_v, d_b, mat2->nlf, mat2->ktmax,
		 mat2->ltmtx, mat2->ndt, mat2->ndl, mat2->nstrtl, mat2->nstrtt, mat2->kt,
		 mat2->a1, mat2->a2, d_sm.rowmat, d_sm.rowmat_t,
		 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
		 v, b, nd, fname, nbench,
		 mat2);
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
  cudaFree(d_b);
  cudaFree(d_v);

  delete [] v;
  printf("hmvm_magma: end\n");
}


// ######## ######## ######## ######## ######## ######## ######## ########
// template関数の実体化のための宣言
// ######## ######## ######## ######## ######## ######## ######## ########
template void hmvm_magma<float>(matrix2<float>  *mat2, float *b, int kernel, int dump_result, int nbench);
template void hmvm_magma<double>(matrix2<double> *mat2, double *b, int kernel, int dump_result, int nbench);
