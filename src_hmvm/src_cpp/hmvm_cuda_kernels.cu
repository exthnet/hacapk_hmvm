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
  バリエーション：
  - a2trans：approxy2を転置版で計算するか否か(0,1)
  - a2interchange：approxy2のループを入れ替えるか否か(0,1)
*/
template <class T, int a2t, int a2i>
__global__ void hmvm_cuda_seq
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_seq : begin\n");
#endif
  int ndl, ndt, nstrtl, nstrtt;
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
	kt = _kt[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
#endif
	head = a1[ip];
	for(il=0; il<kt; il++){
	  tmp2[il] = 0.0;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp2[il] += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp2[il] += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	}
	head = a2[ip];
	if(a2t==0){
	  if(a2i==0){
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
		  }
		}
	  }else{
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  for(il=0; il<kt; il++){
			itl=it+il*ndl;
			myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
		  }
		}
	  }
	}else{
	  if(a2i==0){
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it*kt+il;
			myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
		  }
		}
	  }else{
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  for(il=0; il<kt; il++){
			itl=it*kt+il;
			myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
		  }
		}
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
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
#endif
	head = a1[ip];
	for(il=0; il<ndl; il++){
	  T tmp = (T)0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	  myAtomicAdd(&d_zaut[ill], tmp);
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_seq : end\n");
#endif
}

template <class T, int a2t, int a2i>
void hmvm_cuda_seq_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	hmvm_cuda_seq<T,a2t,a2i><<<1,1,sizeof(T)*ktmax>>>
	  (d_zaut, d_zu, nlf, ktmax,
	   ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense);
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
	printf("TIME %d hmvm_cuda1_seq%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}

void hmvm_cuda_seq_proxy
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, float *rowmat, float *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 float *v, float *b, int nd, char *fname, int bench,
 int a2t, int a2i){
  if(a2t==0 && a2i==0)
	hmvm_cuda_seq_proxy<float,0,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==0 && a2i==1)
	hmvm_cuda_seq_proxy<float,0,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==0)
	hmvm_cuda_seq_proxy<float,1,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==1)
	hmvm_cuda_seq_proxy<float,1,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
}
void hmvm_cuda_seq_proxy
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, double *rowmat, double *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 double *v, double *b, int nd, char *fname, int bench,
 int a2t, int a2i){
  if(a2t==0 && a2i==0)
	hmvm_cuda_seq_proxy<double,0,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==0 && a2i==1)
	hmvm_cuda_seq_proxy<double,0,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==0)
	hmvm_cuda_seq_proxy<double,1,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==1)
	hmvm_cuda_seq_proxy<double,1,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt, _kt,
	   a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
}
#endif
// ######## ######## ######## ######## ######## ######## ######## ########

// ######## ######## ######## ######## ######## ######## ######## ########
/*
  block並列化
  ThreadBlockごとに1つの部分行列積(mat-mat-vecまたはmat-vec)を行う
  ThreadBlock内部は逐次
  - a2trans：approxy2を転置版で計算するか否か(0,1)
  - a2interchange：approxy2のループを入れ替えるか否か(0,1)
 */
template <class T, int a2t, int a2i>
__global__ void hmvm_cuda_block
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD_block : begin\n");
#endif
  int ndl, ndt, nstrtl, nstrtt;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);

  if(blockIdx.x<napprox){
#ifndef _SKIP_APPROX
	// approx
	//for(i=0; i<napprox; i++){
	ip = approx[blockIdx.x];
	ndl = _ndl[ip];
	ndt = _ndt[ip];
	nstrtl = _nstrtl[ip];
	nstrtt = _nstrtt[ip];
	kt = _kt[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
#endif
	head = a1[ip];
	for(il=0; il<kt; il++){
	  tmp2[il] = 0.0;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp2[il] += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp2[il] += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	}
	head = a2[ip];
	if(a2t==0){
	  if(a2i==0){
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
		  }
		}
	  }else{
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  for(il=0; il<kt; il++){
			itl=it+il*ndl;
			myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
		  }
		}
	  }
	}else{
	  if(a2i==0){
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it*kt+il;
			myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
		  }
		}
	  }else{
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  for(il=0; il<kt; il++){
			itl=it*kt+il;
			myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
		  }
		}
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
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
#endif
	head = a1[ip];
	for(il=0; il<ndl; il++){
	  T tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	  myAtomicAdd(&d_zaut[ill], tmp);
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD_block : end\n");
#endif
}

template <class T, int a2t, int a2i>
void hmvm_cuda_block_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	hmvm_cuda_block<T,a2t,a2i><<<napprox+ndense,1,sizeof(T)*ktmax>>>
	  (d_zaut, d_zu, nlf, ktmax,
	   ltmtx, ndt, ndl, nstrtl, nstrtt,
	   kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense);
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
	printf("TIME %d hmvm_cuda1_block%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}

void hmvm_cuda_block_proxy
(float *d_zaut, float *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, float *rowmat, float *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 float *v, float *b, int nd, char *fname, int bench,
 int a2t, int a2i){
  if(a2t==0 && a2i==0)
	hmvm_cuda_block_proxy<float,0,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==0 && a2i==1)
	hmvm_cuda_block_proxy<float,0,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==0)
	hmvm_cuda_block_proxy<float,1,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==1)
	hmvm_cuda_block_proxy<float,1,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
}
void hmvm_cuda_block_proxy
(double *d_zaut, double *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, double *rowmat, double *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 double *v, double *b, int nd, char *fname, int bench,
 int a2t, int a2i){
  if(a2t==0 && a2i==0)
	hmvm_cuda_block_proxy<double,0,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==0 && a2i==1)
	hmvm_cuda_block_proxy<double,0,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==0)
	hmvm_cuda_block_proxy<double,1,0>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
  if(a2t==1 && a2i==1)
	hmvm_cuda_block_proxy<double,1,1>
	  (d_zaut, d_zu, nlf, ktmax,
	   _ltmtx, _ndt, _ndl, _nstrtl, _nstrtt,
	   _kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense,
	   v, b, nd, fname, bench);
}
// ######## ######## ######## ######## ######## ######## ######## ########

// ######## ######## ######## ######## ######## ######## ######## ########
#if 1
/*
  hybrid1
  TBあたりスレッド数は32に固定
  1TBが1つのPMVを担当
  PMV内の1行を1/div WARPが担当
  <<<napprox+ndense,32>>>
  1 PMV by 1 TB(=1WARP)
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - a2t：a2を転置版で計算するか否か(0,1)
  - a2interchange：a2のループを入れ替えるか否か(0,1)
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  6x2x2x2x2=96通り
  divが大きな時にちょっとおかしいかも？
*/
template <class T, int div, int a2t, int a2i, int aatomic, int datomic>
__global__ void hmvm_cuda_hybrid1
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid1 : begin\n");
#endif
  int gid   = blockIdx.x;
  int bid   = threadIdx.x/(32/div);
  int blen  = div;
  int xid   = threadIdx.x%(32/div);
  int xlen  = (32/div);
  int ndl, ndt, nstrtl, nstrtt;
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
	kt = _kt[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
#endif
	head = a1[ip];
	for(il=bid; il<kt; il+=blen){
	  if(xid==0)tmp2[il] = 0.0;
	  tmp = 0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	  for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
	  if(xid==0)tmp2[il] = tmp;
	}
	head = a2[ip];
	if(a2t==0){
	  if(a2i==0){
		for(il=bid; il<kt; il+=blen){
		  for(it=xid; it<ndl; it+=xlen){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
		  }
		}
	  }else{ // a2i==1
		if(aatomic==0){
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			tmp = 0.0;
			for(il=xid; il<kt; il+=xlen){
			  itl=it+il*ndl;
			  tmp += rowmat[head+itl]*tmp2[il];
			}
			for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
			if(xid==0){
			  myAtomicAdd(&d_zaut[ill], tmp);
			}
		  }
		}else{ // aatomic==1
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			for(il=xid; il<kt; il+=xlen){
			  itl=it+il*ndl;
			  myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
			}
		  }
		}
	  }
	}else{ // a2t==1
	  if(a2i==0){
		for(il=bid; il<kt; il+=blen){
		  for(it=xid; it<ndl; it+=xlen){
			ill=it+nstrtl-1;
			itl=it*kt+il;
			myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
		  }
		}
	  }else{ // a2i==1
		if(aatomic==0){
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			tmp = 0.0;
			for(il=xid; il<kt; il+=xlen){
			  itl=it*kt+il;
			  tmp += rowmat_t[head+itl]*tmp2[il];
			}
			for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
			if(xid==0){
			  myAtomicAdd(&d_zaut[ill], tmp);
			}
		  }
		}else{ // aatomic==1
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			for(il=xid; il<kt; il+=xlen){
			  itl=it*kt+il;
			  myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
			}
		  }
		}
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
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
#endif
	head = a1[ip];
	for(il=bid; il<ndl; il+=blen){
	  tmp = 0.0;
	  ill=il+nstrtl-1;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	  if(datomic==0){
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

template <class T, int div, int a2t, int a2i, int aa, int da>
void hmvm_cuda_hybrid1_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms,
 T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	hmvm_cuda_hybrid1<T,div,a2t,a2i,aa,da><<<blocks,threads,shms>>>
	  (d_zaut, d_zu, nlf, ktmax, ltmtx,
	   ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense);
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
	printf("TIME %d hmvm_cuda1_hybrid1%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}
#include "template_hybrid1.hpp"
// ######## ######## ######## ######## ######## ######## ######## ########
#endif

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
  - mul：立ち上げるスレッド数(mul*32)、1つのmat-mat-vecまたはmat-vecをmul TBで実行、mul=1,2,3,...,16
  - a2t：a2を転置版で計算するか否か(0,1)
  - a2interchange：a2のループを入れ替えるか否か(0,1)
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  6x16*2x2x2x2=1536通り
*/
template <class T, int div, int mul, int a2t, int a2i, int aatomic, int datomic>
__global__ void hmvm_cuda_hybrid2
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense)
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
  int ndl, ndt, nstrtl, nstrtt;
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
	kt = _kt[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
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
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
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
	  if(datomic==0){
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

#if 0
template <class T, int div>
__global__ void hmvm_cuda_hybrid2
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int mul, int a2t, int a2i, int aatomic, int datomic)
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
  int ndl, ndt, nstrtl, nstrtt;
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
	kt = _kt[ip];
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
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
#if _DEBUG_LEVEL >= 3
	printf("%d: %d %d %d %d\n", ip, ndl, ndt, nstrtl, nstrtt);
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
	  if(datomic==0){
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
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms,
 T *v, T *b, int nd, char *fname, int bench,
 int div, int mul, int a2t, int a2i, int aa, int da)
{
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	switch(div){
	case 1:
	  hmvm_cuda_hybrid2<T,1><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense , mul, a2t, a2i, aa, da);
	  break;
	case 2:
	  hmvm_cuda_hybrid2<T,2><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense , mul, a2t, a2i, aa, da);
	  break;
	case 4:
	  hmvm_cuda_hybrid2<T,4><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense , mul, a2t, a2i, aa, da);
	  break;
	case 8:
	  hmvm_cuda_hybrid2<T,8><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense , mul, a2t, a2i, aa, da);
	  break;
	case 16:
	  hmvm_cuda_hybrid2<T,16><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense , mul, a2t, a2i, aa, da);
	  break;
	case 32:
	  hmvm_cuda_hybrid2<T,32><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense , mul, a2t, a2i, aa, da);
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
	printf("TIME %d hmvm_cuda1_hybrid2%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}
#else
template <class T, int div, int mul, int a2t, int a2i, int aa, int da>
void hmvm_cuda_hybrid2_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms,
 T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	hmvm_cuda_hybrid2<T,div,mul,a2t,a2i,aa,da><<<blocks,threads,shms>>>
	  (d_zaut, d_zu, nlf, ktmax, ltmtx,
	   ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense);
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
	printf("TIME %d hmvm_cuda1_hybrid2%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}
#include "template_hybrid2.hpp"
#endif
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
template <class T, int div, int mul, int a2t, int a2i, int aatomic, int datomic>
__global__ void hmvm_cuda_hybrid3
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid3 : begin\n");
#endif
  int gid   = blockIdx.x*mul+threadIdx.x/32;
  int bid   = ((threadIdx.x%32)/(32/div));
  int blen  = (32/(32/div));
  int xid   = (threadIdx.x%(32/div));
  int xlen  = (32/div);
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  T tmp = 0.0;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  cg::thread_block_tile<32/div> g = cg::tiled_partition<32/div>(cg::this_thread_block());

#if 0
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
	  if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = 0.0;
	  tmp = 0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
	  for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
	  if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = tmp;
	}
	head = a2[ip];
	for(il=bid; il<kt; il+=blen){
	  for(it=xid; it<ndl; it+=xlen){
		ill=it+nstrtl-1;
		itl=it+il*ndl;
		myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[(threadIdx.x/32)*ktmax+il]);
	  }
	}
  }
#endif

#if 0
  //ip = dense[gid-((napprox+mul-1)/mul)];
  if(gid < ndense){
	ip = dense[gid];
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
#endif

#if 1
  if(gid<((napprox+mul-1)/mul)*mul){
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
  	    if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = 0.0;
		tmp = 0.0;
		for(it=xid; it<ndt; it+=xlen){
	      itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp += rowmat[head+itl]*d_zu[itt];
        }
		//for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
		for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
		if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = tmp;
      }
	  head = a2[ip];
	  for(il=bid; il<kt; il+=blen){
	    for(it=xid; it<ndl; it+=xlen){
	      ill=it+nstrtl-1;
		  itl=it+il*ndl;
		  myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[(threadIdx.x/32)*ktmax+il]);
        }
      }
    }
#endif // approx
  }else{
#ifndef _SKIP_DENSE
	ip = gid-((napprox+mul-1)/mul)*mul;
    if(ip<ndense){
	  ip = dense[ip];
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
		if(datomic==0){
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

#if 1
template <class T, int div>
__global__ void hmvm_cuda_hybrid3
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int mul, int a2t, int a2i, int aatomic, int datomic)
{
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid3 : begin\n");
#endif
  int gid   = blockIdx.x*mul+threadIdx.x/32;
  int bid   = ((threadIdx.x%32)/(32/div));
  int blen  = (32/(32/div));
  int xid   = (threadIdx.x%(32/div));
  int xlen  = (32/div);
  int ndl, ndt, nstrtl, nstrtt, ltmtx;
  int ip, kt, il, it, itt, itl, ill;
  size_t head;
  T tmp = 0.0;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  cg::thread_block_tile<32/div> g = cg::tiled_partition<32/div>(cg::this_thread_block());

#if 0
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
	  if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = 0.0;
	  tmp = 0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		tmp += rowmat[head+itl]*d_zu[itt];
	  }
	  //for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
	  for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
	  if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = tmp;
	}
	head = a2[ip];
	for(il=bid; il<kt; il+=blen){
	  for(it=xid; it<ndl; it+=xlen){
		ill=it+nstrtl-1;
		itl=it+il*ndl;
		myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[(threadIdx.x/32)*ktmax+il]);
	  }
	}
  }
#endif

#if 0
  //ip = dense[gid-((napprox+mul-1)/mul)];
  if(gid < ndense){
	ip = dense[gid];
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
	  if(datomic==0){
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
#endif

#if 1
  if(gid<((napprox+mul-1)/mul)*mul){
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
  	    if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = 0.0;
		tmp = 0.0;
		for(it=xid; it<ndt; it+=xlen){
	      itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp += rowmat[head+itl]*d_zu[itt];
        }
		//for (int offset = warpSize/(2*div); offset > 0; offset /= 2)tmp += __shfl_down_sync(0xffff, tmp, offset, warpSize);
		for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
		if(xid==0)tmp2[(threadIdx.x/32)*ktmax+il] = tmp;
      }
	  head = a2[ip];
	  for(il=bid; il<kt; il+=blen){
	    for(it=xid; it<ndl; it+=xlen){
	      ill=it+nstrtl-1;
		  itl=it+il*ndl;
		  myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[(threadIdx.x/32)*ktmax+il]);
        }
      }
    }
#endif // approx
  }else{
#ifndef _SKIP_DENSE
	ip = gid-((napprox+mul-1)/mul)*mul;
    if(ip<ndense){
	  ip = dense[ip];
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
		if(datomic==0){
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
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms,
 T *v, T *b, int nd, char *fname, int bench,
 int div, int mul, int a2t, int a2i, int aa, int da)
{
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	switch(div){
	case  1: hmvm_cuda_hybrid3<T, 1><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul, a2t, a2i, aa, da);
	  break;
	case  2: hmvm_cuda_hybrid3<T, 2><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul, a2t, a2i, aa, da);
	  break;
	case  4: hmvm_cuda_hybrid3<T, 4><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul, a2t, a2i, aa, da);
	  break;
	case  8: hmvm_cuda_hybrid3<T, 8><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul, a2t, a2i, aa, da);
	  break;
	case 16: hmvm_cuda_hybrid3<T,16><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul, a2t, a2i, aa, da);
	  break;
	case 32: hmvm_cuda_hybrid3<T,32><<<blocks,threads,shms>>>
		(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt,
		 a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul, a2t, a2i, aa, da);
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
	printf("TIME %d hmvm_cuda1_hybrid3%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}
#else
template <class T, int div, int mul, int a2t, int a2i, int aa, int da>
void hmvm_cuda_hybrid3_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms,
 T *v, T *b, int nd, char *fname, int bench)
{
#if 1
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	hmvm_cuda_hybrid3<T,div,mul,a2t,a2i,aa,da><<<blocks,threads,shms>>>
	  (d_zaut, d_zu, nlf, ktmax, ltmtx,
	   ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, rowmat_t,
	   napprox, approx, ndense, dense);
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
	printf("TIME %d hmvm_cuda1_hybrid3%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
#endif
}
#include "template_hybrid3.hpp"
#endif
// ######## ######## ######## ######## ######## ######## ######## ########
#endif

#if 1
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
template <class T, int atomic>
__global__ void hmvm_cuda_hybrid4
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt, int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense, int mul)
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
  T tmp = 0.0;
  extern __shared__ __align__(sizeof(T)) unsigned char my_smem[];
  T *tmp2 = reinterpret_cast<T *>(my_smem);
  cg::thread_block_tile<32> g = cg::tiled_partition<32>(cg::this_thread_block());

#if 1
  if(gid<ndense){
#ifndef _SKIP_APPROX
	ip = dense[gid];
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
	  //for (int offset = warpSize/2; offset > 0; offset /= 2)tmp += __shfl_down_sync(tmp, offset, warpSize);
	  for (int offset = g.size()/2; offset > 0; offset /= 2)tmp += g.shfl_down(tmp, offset);
	  if(tid%32==0){
		atomicAdd(&d_zaut[ill], tmp);
	  }
	}
#endif
  }
#endif

#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_hybrid4 : end\n");
#endif
}

template <class T>
void hmvm_cuda_hybrid4_proxy
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *ltmtx, int *ndt, int *ndl, int *nstrtl, int *nstrtt, int *kt, int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense,
 int blocks, int threads, int shms, int mul, int atomic, T *v, T *b, int nd, char *fname, int bench)
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
	  hmvm_cuda_hybrid4<T,0><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul);
	  break;
	case 1:
	  hmvm_cuda_hybrid4<T,1><<<blocks,threads,shms>>>(d_zaut, d_zu, nlf, ktmax, ltmtx, ndt, ndl, nstrtl, nstrtt, kt, a1, a2, rowmat, rowmat_t, napprox, approx, ndense, dense, mul);
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
	printf("TIME %d hmvm_cuda1_hybrid4%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
#endif
}

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
