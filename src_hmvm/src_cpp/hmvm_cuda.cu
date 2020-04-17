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




#define myAtomicAdd atomicAdd



/*
  <<<ndense,32>>>
  1 GEMV by 1 TB
  1 line by 1/div WARP
*/

__global__ void hmvm_cuda_hybrid0
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense, int div)
{
}
// nlf block, 32 thread
__global__ void hmvm_cuda_hybrid0
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense, int div)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid1 : begin\n");
#endif
  int gid   = blockIdx.x;
  int tid   = threadIdx.x;
  int bid   = threadIdx.x/(32/div);
  int blen  = div;
  int xid   = threadIdx.x%(32/div);
  int xlen  = (32/div);
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  double tmp = 0.0;
  extern __shared__ __align__(sizeof(double)) unsigned char my_smem[];
  double *tmp2 = reinterpret_cast<double *>(my_smem);

  if(gid<napprox){
#if 0
	// approx
	ip = approx[gid];
	ndl = _ndl[ip];
	ndt = _ndt[ip];
	nstrtl = _nstrtl[ip];
	nstrtt = _nstrtt[ip];
	ltmtx = _ltmtx[ip];
	kt = _kt[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt, ltmtx);
#endif
	head = a1[ip];
	for(il=bid; il<kt; il+=blen){
	  tmp2[il] = 0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp2[il] += rowmat[head+itl]*d_zu[itt];
	  }
	}
	head = a2[ip];
	for(il=bid; il<kt; il+=blen){
	  for(it=xid; it<ndl; it+=xlen){
		ill=it+nstrtl-1;
		itl=it+il*ndl;
		myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
	  }
	}
#endif
  }else{
#if 1
	// dense
	ip = dense[gid-napprox];
	ndl = _ndl[ip];
	ndt = _ndt[ip];
	nstrtl = _nstrtl[ip];
	nstrtt = _nstrtt[ip];
	ltmtx = _ltmtx[ip];
#if 1//_DEBUG_LEVEL >= 3
	if(ip==dense[0])printf("dense %d %d: %d %d %d %d %d\n", ip, gid, ndl, ndt, nstrtl, nstrtt, ltmtx);
#endif
	head = a1[ip];
	__syncthreads();
	//if(ip==0&&xid==0)printf("x %d %f + %f\n", ill, d_zaut[ill], tmp);
	__syncthreads();
	//	if(blockIdx.x==0&&threadIdx.x==0)printf("! %d %e + %f\n", 0, d_zaut[0], tmp);
	if(ip==dense[0])printf("! %d %e + %f\n", 0, d_zaut[0], tmp);
	//if(ip<ndense) // 判定いるのか？
	for(il=bid; il<ndl; il+=blen){
	  // 何故かだめ
#if 1
	  tmp = 0.0;
	  ill=il+nstrtl-1;
	  if(ip==dense[0])printf("x %d %e + %f\n", ill, d_zaut[ill], tmp);
	  __syncthreads();
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  if(ip==dense[0]&&il==bid)printf("tmp %d %f\n", threadIdx.x, tmp);
	  //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down(tmp, offset);
	  __syncthreads();
	  __syncthreads();
	  for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
	  __syncthreads();
	  if(xid==0){
		if(ip==dense[0]&&il==bid)printf("1->1 %d %f + %f\n", ill, d_zaut[ill], tmp);
		atomicAdd(&d_zaut[ill], tmp);
		if(ip==dense[0]&&il==bid)printf("1->2 %d %f + %f\n", ill, d_zaut[ill], tmp);
	  }
	  __syncthreads();
#endif
	  // ただしくうごく
#if 0
	  tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down(tmp, offset);
	  __syncthreads();
	  if(ip==dense[0]&&il==bid)printf("2->1 %d %e + %e\n", ill, d_zaut[ill], tmp);
	  __syncthreads();
	  atomicAdd(&d_zaut[ill], tmp);
	  __syncthreads();
	  if(ip==dense[0]&&il==bid)printf("2->2 %d %e + %e\n", ill, d_zaut[ill], tmp);
#endif
	}
#endif
  }

#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid1 : end\n");
#endif
}

__global__ void hmvm_cuda_test
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
}
__global__ void hmvm_cuda_test
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD : begin\n");
#endif
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  //extern __shared__ double tmp2[];
  extern __shared__ __align__(sizeof(double)) unsigned char my_smem[];
  double *tmp2 = reinterpret_cast<double *>(my_smem);
  int i;

  // approx
  for(i=0; i<napprox; i++){
#if 1
	ip = approx[i];
	ndl = _ndl[ip];
	ndt = _ndt[ip];
	nstrtl = _nstrtl[ip];
	nstrtt = _nstrtt[ip];
	ltmtx = _ltmtx[ip];
	kt = _kt[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt, ltmtx);
#endif
	head = a1[ip];
	for(il=0; il<kt; il++){
	  tmp2[il] = 0.0;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp2[il] += rowmat[head+itl]*d_zu[itt];
	  }
	}
	head = a2[ip];
	for(il=0; il<kt; il++){
	  for(it=0; it<ndl; it++){
		ill=it+nstrtl-1;
		itl=it+il*ndl;
		myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
	  }
	}
#endif
  }
  // dense
  for(i=0; i<ndense; i++){
#if 1
	ip = dense[i];
	ndl = _ndl[ip];
	ndt = _ndt[ip];
	nstrtl = _nstrtl[ip];
	nstrtt = _nstrtt[ip];
	ltmtx = _ltmtx[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt, ltmtx);
#endif
	head = a1[ip];
	for(il=0; il<ndl; il++){
	  double tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  myAtomicAdd(&d_zaut[ill], tmp);
	  printf("myAtomicAdd %d %e\n", ill, tmp);
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD : end\n");
#endif
}


template<class T>
void hmvm_cuda1(matrix2<T> *mat2, T *b, int kernel, int dump_result)
{
  const int L=10, M=5;
  FILE *F;
  matrix2<T> d_sm;
  int i, l, nd = mat2->nd;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  T *v=NULL, *tmp=NULL, *zero;
  T *d_b, *d_v;//, *d_zaut, *d_zbut;
  int ip;
  int len;
  cudaError_t ret;
  printf("hmvm_cuda1_%s: begin\n", typeid(T).name());
  v    = new T[mat2->nd];    //(double*)malloc(sizeof(double)*mat2->nd);
  tmp  = new T[mat2->ktmax]; //(double*)malloc(sizeof(double)*mat2->ktmax);
  zero = new T[mat2->ktmax]; //(double*)malloc(sizeof(double)*mat2->ktmax);
  for(i=0;i<nd;i++){
	v[i] = (T)0.0;
  }
  for(i=0;i<mat2->ktmax;i++){
	zero[i] = (T)0.0;
  }
  //CHECK_DO(cudaMalloc((void**)&d_zaut, sizeof(T)*mat2->nd),"cudaMalloc z_aut");
  //CHECK_DO(cudaMalloc((void**)&d_zbut, sizeof(T)*mat2->ktmax),"cudaMalloc zbut");
  CHECK_DO(cudaMalloc((void**)&d_b, sizeof(T)*mat2->nd),"cudaMalloc d_b");
  CHECK_DO(cudaMalloc((void**)&d_v, sizeof(T)*mat2->nd),"cudaMallod d_v");
  //for(i=0;i<mat2->nd;i++){d_b[i]=NULL;d_v[i]=NULL;}

  len = mat2->len;
  printf("total length = %d\n", len);
  // host alloc
  // device alloc
  d_sm.nd    = mat2->nd;
  d_sm.nlf   = mat2->nlf;
  d_sm.ktmax = mat2->ktmax;
  cudaMalloc((void**)&d_sm.ltmtx, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.ndl, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.ndt, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.nstrtl, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.nstrtt, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.kt, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.a1, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.a2, sizeof(int)*mat2->nlf);
  cudaMalloc((void**)&d_sm.rowmat,sizeof(T)*mat2->len);
  cudaMalloc((void**)&d_sm.rowmat_t,sizeof(T)*mat2->len);
  // memcpy
  CHECK_DO(cudaMemcpy(d_sm.ltmtx, mat2->ltmtx, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.ltmtx");
  CHECK_DO(cudaMemcpy(d_sm.ndt, mat2->ndt, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.ndt");
  CHECK_DO(cudaMemcpy(d_sm.ndl, mat2->ndl, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.ndl");
  CHECK_DO(cudaMemcpy(d_sm.nstrtl, mat2->nstrtl, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.nstrtl");
  CHECK_DO(cudaMemcpy(d_sm.nstrtt, mat2->nstrtt, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.nstrtt");
  CHECK_DO(cudaMemcpy(d_sm.kt, mat2->kt, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.kt");
  CHECK_DO(cudaMemcpy(d_sm.a1, mat2->a1, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.a1");
  CHECK_DO(cudaMemcpy(d_sm.a2, mat2->a2, sizeof(int)*mat2->nlf, cudaMemcpyHostToDevice),"d_sm.a2");
  CHECK_DO(cudaMemcpy(d_sm.rowmat, mat2->rowmat, sizeof(T)*mat2->len, cudaMemcpyHostToDevice),"d_sm.rowmat");
  CHECK_DO(cudaMemcpy(d_sm.rowmat_t, mat2->rowmat_t, sizeof(T)*mat2->len, cudaMemcpyHostToDevice),"d_sm.rowmat_t");

  // 分離
  printf("begin splitting\n");
  mat2->approx = (int*)malloc(sizeof(int)*mat2->nlf);
  mat2->dense  = (int*)malloc(sizeof(int)*mat2->nlf);
  mat2->napprox = mat2->ndense = 0;
  for(ip=0; ip<mat2->nlf; ip++){
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

#define EXEC(FUNCNAME,BLOCKS,THREADS,S)											\
  printf("nd = %d\n", nd);												\
  for(i=0;i<nd;i++)v[i] = (T)0.0;										\
  cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice);				\
  cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice);				\
  FUNCNAME<<<BLOCKS,THREADS,S>>>													\
	(d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,			\
	 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,							\
	 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
  cudaDeviceSynchronize();												\
  cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost);				\
  printf("write to %s\n", fname);										\
  F = fopen(fname, "w");												\
  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);							\
  fclose(F);

#define BENCH(FUNCNAME,B,T,S)											\
  printf("nd = %d\n", nd);												\
  for(l=0;l<L;l++){														\
	for(i=0;i<nd;i++)v[i] = (T)0.0;										\
	cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice);			\
	cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice);			\
	cudaDeviceSynchronize();											\
	d1 = omp_get_wtime();												\
	FUNCNAME<<<B,T,S>>>													\
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
  davg /= (L-5);

  // sequential
  if(kernel==0)
  {
	int a1, a2;
	char name[8], fname[0xff];
	a1=a2=0;
	snprintf(name,8,"_%d_%d",a1,a2);
	snprintf(fname,32,"result_cuda1%s_%s.txt", name, typeid(T).name());
	printf("fname = %s\n", fname);
	/*
	printf("nd = %d\n", nd);
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice);
	hmvm_cuda_seq<<<1,1,d_sm.ktmax*sizeof(T)>>>
	  //hmvm_cuda_test<<<1,1,d_sm.ktmax*sizeof(T)>>>
	  //hmvm_cuda_test<<<1,1,d_sm.ktmax*sizeof(double)>>>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);
	cudaDeviceSynchronize();
	cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost);
	printf("write to %s\n", fname);
	F = fopen(fname, "w");
	for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	fclose(F);
	*/
	EXEC(hmvm_cuda_seq,1,1,d_sm.ktmax*sizeof(T));
	//BENCH(hmvm_cuda_seq<T>,1,1,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }

  // block parallel
  // whole hmvm calculation in 1 GPU kernel
  if(kernel==1)
  {
	int a1, a2;
	char name[8], fname[32];
	a1=a2=0;
	snprintf(name,8,"_%d_%d",a1,a2);
	snprintf(fname,32,"result_cuda1blk%s_%s.txt", name, typeid(T).name());
	printf("fname = %s\n", fname);
	/*
	printf("nd = %d\n", nd);
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice);
	hmvm_cuda_block<<<d_sm.napprox+d_sm.ndense,1,d_sm.ktmax*sizeof(T)>>>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);
	cudaDeviceSynchronize();
	cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost);
	printf("write to %s\n", fname);
	F = fopen(fname, "w");
	for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	fclose(F);
	*/
	EXEC(hmvm_cuda_block,d_sm.napprox+d_sm.ndense,1,d_sm.ktmax*sizeof(T));
	//BENCH(hmvm_cuda_block<T>,d_sm.napprox+d_sm.ndense,1,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1blk%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }

#define EXEC_1(FUNCNAME,BK,TH,S)											\
  printf("nd = %d\n", nd);												\
  for(i=0;i<nd;i++)v[i] = (T)0.0;										\
  cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice);				\
  cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice);				\
  FUNCNAME<T,1><<<BK,TH,S>>>												\
	(d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,			\
	 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,							\
	 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
  cudaDeviceSynchronize();												\
  cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost);				\
  printf("write to %s\n", fname);										\
  F = fopen(fname, "w");												\
  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);							\
  fclose(F);

  // block & thread parallel
  // whole hmvm calculation in 1 GPU kernel
  // under development
  if(kernel==2)
  {
	int a1, a2;
	char name[8], fname[32];
	a1=a2=0;
	snprintf(name,8,"_%d_%d",a1,a2);
	snprintf(fname,32,"result_cuda1hyb%s_%s.txt", name, typeid(T).name());
	printf("fname = %s\n", fname);

	printf("nd = %d\n", nd);
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice); // 初期化したはずなのに参照するとおかしくない？
	cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice);
	printf("launch %d blocks\n", d_sm.napprox+d_sm.ndense);
	//hmvm_cuda_hybrid0<<<d_sm.napprox+1,32,d_sm.ktmax*sizeof(T)>>>
	hmvm_cuda_hybrid0<<<d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T)>>>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense, 1);
	/*
	hmvm_cuda_hybrid1<T,1><<<d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T)>>>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);
	*/
	cudaDeviceSynchronize();
	cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost);
	printf("write to %s\n", fname);
	F = fopen(fname, "w");
	for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	fclose(F);
	//EXEC_1(hmvm_cuda_hybrid1,d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T));
	/*
	proxy_hmvm_cuda_hybrid1<T,1><<<d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T)>>>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);
	*/
	//BENCH(hmvm_cuda_hybrid<T>,d_sm.napprox+d_sm.ndense,1,d_sm.ktmax*sizeof(T));
	printf("TIME %d hmvm_cuda1hyb%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
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
