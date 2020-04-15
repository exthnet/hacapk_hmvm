// -*- c++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <cuda_runtime_api.h>

#include "hacapk.h"

#if __CUDA_ARCH__ < 600
__device__ double myAtomicAdd(double* address, double val)
{
  unsigned long long int* address_as_ull =
	(unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
	assumed = old;
	old = atomicCAS(address_as_ull, assumed,
					__double_as_longlong(val +
										 __longlong_as_double(assumed)));

	// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

  return __longlong_as_double(old);
}
__device__ float myAtomicAdd(float* address, float val)
{
  atomicAdd(address, val);
}
#else
#define myAtomicAdd atomicAdd
#endif

// ######## ######## ######## ######## ######## ######## ######## ########
// 0: 完全逐次

template <class T>
__global__ void hmvm_cuda_seq
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD : begin\n");
#endif
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  T tmp;
  //extern __shared__ T tmp2[];
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  int i;

  // approx
  for(i=0; i<napprox; i++){
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
  }
  // dense
  for(i=0; i<ndense; i++){
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
	  tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  myAtomicAdd(&d_zaut[ill], tmp);
	}
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD : end\n");
#endif
}

// ######## ######## ######## ######## ######## ######## ######## ########

// block並列化カーネル
template <class T>
__global__ void hmvm_cuda_block
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD_block : begin\n");
#endif
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  T tmp;
  //extern __shared__ T tmp2[];
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  int i;

  if(blockIdx.x<napprox){
	// approx
	//for(i=0; i<napprox; i++){
	ip = approx[blockIdx.x];
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
  }else{
	// dense
	//for(i=0; i<ndense; i++){
	ip = dense[blockIdx.x - napprox];
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
	  tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  myAtomicAdd(&d_zaut[ill], tmp);
	}
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD_block : end\n");
#endif
}

// ######## ######## ######## ######## ######## ######## ######## ########

// まだデバッグ中

#if 1
/*
  <<<ndense,32>>>
  1 GEMV by 1 TB
  1 line by 1/div WARP
*/

// nlf block, 32 thread
template <class T, int div>
__global__ void hmvm_cuda_hybrid1
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
  return;
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
  T tmp;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);

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
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt, ltmtx);
#endif
	head = a1[ip];
	for(il=bid; il<ndl; il+=blen){
	  tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down(tmp, offset);
	  for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(tmp, offset, warpSize);
	  if(xid==0){
		myAtomicAdd(&d_zaut[ill], tmp);
	  }
	}
#endif
  }

#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid1 : end\n");
#endif
}
#endif


// ######## ######## ######## ######## ######## ######## ######## ########
// template関数の実体化のための宣言
// ######## ######## ######## ######## ######## ######## ######## ########
template __global__ void hmvm_cuda_seq<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_seq<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);

template __global__ void hmvm_cuda_block<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_block<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);

template __global__ void hmvm_cuda_hybrid1<double,1>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_hybrid1<double,2>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_hybrid1<double,4>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_hybrid1<float,1>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_hybrid1<float,2>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_hybrid1<float,4>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense);
