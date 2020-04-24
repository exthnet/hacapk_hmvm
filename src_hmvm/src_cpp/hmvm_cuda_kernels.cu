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

namespace cg = cooperative_groups;

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

#if 1
// ######## ######## ######## ######## ######## ######## ######## ########
/*
  完全逐次
  バリエーション：なし
*/
template <class T>
__global__ void hmvm_cuda_seq
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_seq : begin\n");
#endif
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  //extern __shared__ T tmp2[];
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  int i;

  // approx
  for(i=0; i<napprox; i++){
#ifndef _SKIP_APPROX
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
#ifndef _SKIP_DENSE
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
	  T tmp = (T)0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  myAtomicAdd(&d_zaut[ill], tmp);
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_seq : end\n");
#endif
}

template __global__ void hmvm_cuda_seq<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_seq<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);
#endif
// ######## ######## ######## ######## ######## ######## ######## ########

// ######## ######## ######## ######## ######## ######## ######## ########
/*
  block並列化
  ThreadBlockごとに1つの部分行列積(mat-mat-vecまたはmat-vec)を行う
  ThreadBlock内部は逐次
  バリエーション：なし
 */
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
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  int i;

  if(blockIdx.x<napprox){
#ifndef _SKIP_APPROX
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
#endif
  }else{
#ifndef _SKIP_DENSE
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
	  T tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  myAtomicAdd(&d_zaut[ill], tmp);
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD_block : end\n");
#endif
}

template __global__ void hmvm_cuda_block<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense);
template __global__ void hmvm_cuda_block<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense);
// ######## ######## ######## ######## ######## ######## ######## ########

// ######## ######## ######## ######## ######## ######## ######## ########
#if 1
/*
  hybrid1
  TBあたりスレッド数は32に固定
  1TBが1つのPMVを担当
  PMV内の1行を1/div WARPが担当
  <<<napprox+ndense,32>>>
  1 GEMV by 1 TB(=1WARP)
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - a2t：a2を転置版で計算するか否か
  - a2interchange：a2のループを入れ替えるか否か
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか
*/
template <class T, int div, int atomic>
__global__ void hmvm_cuda_hybrid1
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid1 : begin\n");
#endif
  int gid   = blockIdx.x;
  //int tid   = threadIdx.x;
  int bid   = threadIdx.x/(32/div);
  int blen  = div;
  int xid   = threadIdx.x%(32/div);
  int xlen  = (32/div);
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  T tmp = 0.0;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  cg::thread_block_tile<32/div> g = cg::tiled_partition<32/div>(cg::this_thread_block());

  if(gid<napprox){
#ifndef _SKIP_APPROX
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
	  if(xid==0)tmp2[il] = 0.0;
	  tmp = 0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
	  for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
	  if(xid==0)tmp2[il] = tmp;
	}
	head = a2[ip];
	for(il=bid; il<kt; il+=blen){
	  for(it=xid; it<ndl; it+=xlen){
		ill=it+nstrtl-1;
		itl=it+il*ndl;
		myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
	  }
	}
#endif // approx
  }else{
#ifndef _SKIP_DENSE
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
	  if(atomic==0){
		//for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
		for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
		if(xid==0){
		  myAtomicAdd(&d_zaut[ill], tmp);
		}
	  }else{
		myAtomicAdd(&d_zaut[ill], tmp);
	  }
	}
#endif // dense
  }

#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid1 : end\n");
#endif
}

template <class T>
void hmvm_cuda_hybrid1_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int atomic, T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  const int L=10, M=5;
  FILE *F;
  int i, l, lmax;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  cudaError_t ret;
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	switch(atomic){
	case 0:
	  switch(div){
	  case  1: hmvm_cuda_hybrid1<T, 1,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case  2: hmvm_cuda_hybrid1<T, 2,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case  4: hmvm_cuda_hybrid1<T, 4,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case  8: hmvm_cuda_hybrid1<T, 8,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case 16: hmvm_cuda_hybrid1<T,16,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case 32: hmvm_cuda_hybrid1<T,32,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  }
	  break;
	case 1:
	  switch(div){
	  case  1: hmvm_cuda_hybrid1<T, 1,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case  2: hmvm_cuda_hybrid1<T, 2,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case  4: hmvm_cuda_hybrid1<T, 4,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case  8: hmvm_cuda_hybrid1<T, 8,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case 16: hmvm_cuda_hybrid1<T,16,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  case 32: hmvm_cuda_hybrid1<T,32,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense); break;
	  }
	  break;
	}
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
	printf("TIME %d hmvm_cuda1_hybrid1b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
#endif
}

template
void hmvm_cuda_hybrid1_proxy<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int atomic, float *v, float *b, int nd, char *fname, int bench);
template
void hmvm_cuda_hybrid1_proxy<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int atomic, double *v, double *b, int nd, char *fname, int bench);
// ######## ######## ######## ######## ######## ######## ######## ########
#endif

/*
template __global__ void hmvm_cuda_hybrid2<float,DIV,MUL,ATOMIC>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
int napprox, int *approx, int ndense, int *dense);

for a in 0 1
do
for m in `seq 1 16`
do
for d in 1 2 4 8 16 32
do
echo "template __global__ void hmvm_cuda_hybrid2<float,${d},${m},${a}>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
int napprox, int *approx, int ndense, int *dense);"
done
done
done
*/

/*
 template __global__ void hmvm_cuda_hybrid2<double,DIV,MUL,ATOMIC>
(T *d_zaut, T *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
int napprox, int *approx, int ndense, int *dense);

for a in 0 1
do
for m in `seq 1 16`
do
for d in 1 2 4 8 16 32
do
echo "template __global__ void hmvm_cuda_hybrid2<double,${d},${m},${a}>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
int napprox, int *approx, int ndense, int *dense);"
done
done
done
*/

#if 1
// ######## ######## ######## ######## ######## ######## ######## ########
/*
  hybrid2
  複数WARP単一GEMV個別行カーネル
  1PMVを担当するのは1TBのまま
  1TBあたりスレッド数を32*mulに増やす
  PMV内の1行を1/div WARPが担当
  <<<napprox+ndense,32*mul>>>
  1 PMV by 1 TB
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - mul：立ち上げるスレッド数(mul*32)、1つのmat-mat-vecまたはmat-vecをmul TBで実行、mul=1,2,3,...,
  - a2t：a2を転置版で計算するか否か
  - a2interchange：a2のループを入れ替えるか否か
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか
*/
template <class T, int div, int atomic>//int div, int mul, int atomic>
__global__ void hmvm_cuda_hybrid2
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense, int mul)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid2 : begin\n");
#endif
  int gid   = blockIdx.x;
  //int tid   = threadIdx.x;
  int bid   = threadIdx.x/(32/div);
  int blen  = mul*div;
  int xid   = threadIdx.x%(32/div);
  int xlen  = (32/div);
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  T tmp = 0.0;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  cg::thread_block_tile<32/div> g = cg::tiled_partition<32/div>(cg::this_thread_block());

  if(gid<napprox){
#ifndef _SKIP_APPROX
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
	  if(xid==0)tmp2[il] = 0.0;
	  tmp = 0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
	  for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
	  if(xid==0)tmp2[il] = tmp;
	}
	head = a2[ip];
	for(il=bid; il<kt; il+=blen){
	  for(it=xid; it<ndl; it+=xlen){
		ill=it+nstrtl-1;
		itl=it+il*ndl;
		myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
	  }
	}
#endif // approx
  }else if(gid-napprox<ndense){
#ifndef _SKIP_DENSE
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
	  if(atomic==0){
		//for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(tmp, offset, warpSize);
		for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
		if(xid==0){
		  myAtomicAdd(&d_zaut[ill], tmp);
		}
	  }else{
		myAtomicAdd(&d_zaut[ill], tmp);
	  }
	}
#endif // dense
  }

#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid2 : end\n");
#endif
}

template <class T>
void hmvm_cuda_hybrid2_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  const int L=10, M=5;
  FILE *F;
  int i, l, lmax;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  cudaError_t ret;
  if(bench==0){lmax=L;}else{lmax=1;}
  for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	switch(atomic){
	case 0:
	  switch(div){
	  case  1: hmvm_cuda_hybrid2<T, 1,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  2: hmvm_cuda_hybrid2<T, 2,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  4: hmvm_cuda_hybrid2<T, 4,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  8: hmvm_cuda_hybrid2<T, 8,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 16: hmvm_cuda_hybrid2<T,16,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 32: hmvm_cuda_hybrid2<T,32,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  }
	  break;
	case 1:
	  switch(div){
	  case  1: hmvm_cuda_hybrid2<T, 1,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  2: hmvm_cuda_hybrid2<T, 2,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  4: hmvm_cuda_hybrid2<T, 4,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  8: hmvm_cuda_hybrid2<T, 8,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 16: hmvm_cuda_hybrid2<T,16,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 32: hmvm_cuda_hybrid2<T,32,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  }
	  break;
	}
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
	printf("TIME %d hmvm_cuda1_hybrid2b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
#endif
}

template
void hmvm_cuda_hybrid2_proxy<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, float *v, float *b, int nd, char *fname, int bench);
template
void hmvm_cuda_hybrid2_proxy<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, double *v, double *b, int nd, char *fname, int bench);
// ######## ######## ######## ######## ######## ######## ######## ########
#endif

#if 1
// ######## ######## ######## ######## ######## ######## ######## ########
/*
  hybrid3
  複数WARP個別GEMVカーネル
  1PMVを1WARPが担当
  1TBあたりスレッド数は32*mul(mul WARP)
  PMV内の1行を1/div WARPが担当
  <<<napprox/mul+ndense/mul, 32*mul>>>
  1 PMV by 1 WARP
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - mul：立ち上げるブロック数スレッド数の係数、1つのmat-mat-vecまたはmat-vecを1WARPで実行、mul=1,2,3,...,
  - a2t：a2を転置版で計算するか否か
  - a2interchange：a2のループを入れ替えるか否か
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか
*/
template <class T, int div, int atomic>
__global__ void hmvm_cuda_hybrid3
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense, int mul)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid3 : begin\n");
#endif
  int gid   = blockIdx.x*mul+threadIdx.x/32;
  //int tid   = threadIdx.x;
  int bid   = ((threadIdx.x/mul)/32)/div;
  int blen  = div;
  int xid   = threadIdx.x%(32/div);
  int xlen  = (32/div);
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  T tmp = 0.0;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  cg::thread_block_tile<32/div> g = cg::tiled_partition<32/div>(cg::this_thread_block());

#if 1
  if(gid<((napprox+mul-1)/mul)){
#ifndef _SKIP_APPROX
	if(gid<napprox){
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
  	    if(xid==0)tmp2[il] = 0.0;
		tmp = 0.0;
		for(it=xid; it<ndt; it+=xlen){
	      itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp += rowmat[head+itl]*d_zu[itt];
        }
		//for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
		for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
		if(xid==0)tmp2[il] = tmp;
      }
	  head = a2[ip];
	  for(il=bid; il<kt; il+=blen){
	    for(it=xid; it<ndl; it+=xlen){
	      ill=it+nstrtl-1;
		  itl=it+il*ndl;
		  myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
        }
      }
    }
#endif // approx
  }else{
#ifndef _SKIP_DENSE
	if(gid-((napprox+mul-1)/mul)<ndense){
	ip = dense[gid-((napprox+mul-1)/mul)];
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
	  if(atomic==0){
	    //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down(tmp, offset);
	    for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
		if(xid==0){
	      atomicAdd(&d_zaut[ill], tmp);
        }
      }else{
  	    atomicAdd(&d_zaut[ill], tmp);
      }
    }
    }
#endif // dense
  }
#endif
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid3 : end\n");
#endif
}

template <class T>
void hmvm_cuda_hybrid3_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  const int L=10, M=5;
  FILE *F;
  int i, l, lmax;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  cudaError_t ret;
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	switch(atomic){
	case 0:
	  switch(div){
	  case  1: hmvm_cuda_hybrid3<T, 1,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  2: hmvm_cuda_hybrid3<T, 2,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  4: hmvm_cuda_hybrid3<T, 4,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  8: hmvm_cuda_hybrid3<T, 8,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 16: hmvm_cuda_hybrid3<T,16,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 32: hmvm_cuda_hybrid3<T,32,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  }
	  break;
	case 1:
	  switch(div){
	  case  1: hmvm_cuda_hybrid3<T, 1,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  2: hmvm_cuda_hybrid3<T, 2,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  4: hmvm_cuda_hybrid3<T, 4,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case  8: hmvm_cuda_hybrid3<T, 8,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 16: hmvm_cuda_hybrid3<T,16,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  case 32: hmvm_cuda_hybrid3<T,32,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, napprox, approx, ndense, dense, mul); break;
	  }
	  break;
	}
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
	printf("TIME %d hmvm_cuda1_hybrid3b%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
#endif
}

template
void hmvm_cuda_hybrid3_proxy<float>
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, float *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, float *v, float *b, int nd, char *fname, int bench);
template
void hmvm_cuda_hybrid3_proxy<double>
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, double *rowmat,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int div, int mul, int atomic, double *v, double *b, int nd, char *fname, int bench);
// ######## ######## ######## ######## ######## ######## ######## ########
#endif

#if 0
// ######## ######## ######## ######## ######## ######## ######## ########
/*
  hybrid4
  複数WARP単一GEMV単一行GEMVカーネル
  <<<napprox+ndense, 32*mul>>>
  1 GEMV by mul TB
  1 line by 1/mul WARP
  バリエーション
  - mul：立ち上げるブロック数スレッド数の係数、1つのmat-mat-vecまたはmat-vecを1WARPで実行、mul=1,2,3,...,
  - a2t：a2を転置版で計算するか否か
  - a2interchange：a2のループを入れ替えるか否か
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか
*/
template <class T>
__global__ void hmvm_cuda_hybrid4
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid4 : begin\n");
#endif
  int gid  = blockIdx.x;
  int tid  = threadIdx.x;
  int tlen = blockDim.x;
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  double tmp;

  ip = dense[gid];
  {
#ifndef _SKIP_APPROX
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
	  for(it=tid; it<ndt; it+=tlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  for (int offset = warpSize/2; offset > 0; offset /= 2)tmp += __shfl_down_sync(tmp, offset, warpSize);
	  if(tid%32==0){
		atomicAdd(&d_zaut[ill], tmp);
	  }
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid4 : end\n");
#endif
}
// ######## ######## ######## ######## ######## ######## ######## ########
#endif
