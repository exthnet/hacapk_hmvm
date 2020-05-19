// -*- C++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>
#include <cuda_runtime_api.h>
#include <cooperative_groups.h>

#include "hacapk.h"

#define CHECK_DO(act,msg) {ret=act; if(ret!=cudaSuccess){printf("%s failed\n",msg);exit(-1);};}

namespace cg = cooperative_groups;

#if __CUDA_ARCH__ < 600
__device__ static double myAtomicAdd(double* address, double val)
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
__device__ static inline float myAtomicAdd(float* address, float val)
{
  atomicAdd(address, val);
}
#else
#define myAtomicAdd atomicAdd
#endif

// ######## ######## ######## ######## ######## ######## ######## ########
/*
  hybrid1
  TBあたりスレッド数は32に固定
  1TBが1つのPMVを担当
  PMV内の1行を1/div WARPが担当
  <<<napprox+ndense,32>>>
  1 PMV by 1 TB (1 TB = 1 WARP)
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - a2t：a2を転置版で計算するか否か(0,1)
  - a2interchange：a2のループを入れ替えるか否か(0,1)
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  6x2x2x2x2=96通り
*/
template <class T, int div, int a2t, int a2i, int aatomic, int datomic>
__global__ void hmvm_cuda_hybrid1
(T *d_zaut, T *d_zu, int nlf, int ktmax,
 int *_ltmtx, int *_ndt, int *_ndl, int *_nstrtl, int *_nstrtt, int *_kt,
 int *a1, int *a2, T *rowmat, T *rowmat_t,
 int napprox, int *approx, int ndense, int *dense
 )
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
  T tmp1 = (T)0.0;
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
	  if(xid==0)tmp2[il] = (T)0.0;
	  tmp1 = (T)0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp1 += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp1 += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	  for (int offset = g.size()/2; offset > 0; offset /= 2)tmp1 += g.shfl_down(tmp1, offset);
	  if(xid==0)tmp2[il] = tmp1;
	}
	__syncwarp();
	head = a2[ip];
	if(a2t==0){ // a2t==0
	  if(a2i==0){ // a2i==0
		for(il=bid; il<kt; il+=blen){
		  for(it=xid; it<ndl; it+=xlen){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
		  }
		}
	  }else{ // a2i==1
		if(aatomic==0){ // aatomic==0
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			tmp1 = (T)0.0;
			for(il=xid; il<kt; il+=xlen){
			  itl=it+il*ndl;
			  tmp1 += rowmat[head+itl]*tmp2[il];
			}
			for (int offset = g.size()/2; offset > 0; offset /= 2)tmp1 += g.shfl_down(tmp1, offset);
			if(xid==0){
			  myAtomicAdd(&d_zaut[ill], tmp1);
			}
		  }
		}else{ // aatomic==1
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			tmp1 = (T)0.0;
			for(il=xid; il<kt; il+=xlen){
			  itl=it+il*ndl;
			  tmp1 += rowmat[head+itl]*tmp2[il];
			}
			myAtomicAdd(&d_zaut[ill], tmp1);
		  }
		}
	  }
	}else{ // a2t==1
	  if(a2i==0){ // a2i==0
		for(il=bid; il<kt; il+=blen){
		  for(it=xid; it<ndl; it+=xlen){
			ill=it+nstrtl-1;
			itl=it*kt+il;
			myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
		  }
		}
	  }else{ // a2i==1
		if(aatomic==0){ // aatomic==0
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			tmp1 = (T)0.0;
			for(il=xid; il<kt; il+=xlen){
			  itl=it*kt+il;
			  tmp1 += rowmat_t[head+itl]*tmp2[il];
			}
			for (int offset = g.size()/2; offset > 0; offset /= 2)tmp1 += g.shfl_down(tmp1, offset);
			if(xid==0){
			  myAtomicAdd(&d_zaut[ill], tmp1);
			}
		  }
		}else{ // aatomic==1
		  for(it=bid; it<ndl; it+=blen){
			ill=it+nstrtl-1;
			tmp1 = (T)0.0;
			for(il=xid; il<kt; il+=xlen){
			  itl=it*kt+il;
			  tmp1 += rowmat_t[head+itl]*tmp2[il];
			}
			myAtomicAdd(&d_zaut[ill], tmp1);
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
	  ill = il+nstrtl-1;
	  tmp1 = (T)0.0;
	  for(it=xid; it<ndt; it+=xlen){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){ // a2t==0
		  tmp1 += rowmat[head+itl]*d_zu[itt];
		}else{ // a2t==1
		  tmp1 += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	  if(datomic==0){ // datomic==0
		for (int offset = g.size()/2; offset > 0; offset /= 2)tmp1 += g.shfl_down(tmp1, offset);
		if(xid==0){
		  myAtomicAdd(&d_zaut[ill], tmp1);
		}
	  }else{ // datomic==1
		myAtomicAdd(&d_zaut[ill], tmp1);
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
	printf("TIME %d hmvm_cuda_hybrid1%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
}
#include "template_hybrid1.hpp"
// ######## ######## ######## ######## ######## ######## ######## ########

template<class T>
void hmvm_cuda1(matrix2<T> *mat2, T *b, int kernel, int dump_result)
{
  matrix2<T> d_sm;
  int i, nd = mat2->nd, ktmax = mat2->ktmax, nlf = mat2->nlf;
  T *v=NULL;
  T *d_b, *d_v;
  int ip;
  int len;
  cudaError_t ret;
  printf("hmvm_cuda1_%s: begin\n", typeid(T).name()); fflush(stdout);
  v    = new T[nd];
  for(i=0;i<nd;i++){
	v[i] = (T)0.0;
  }
  CHECK_DO(cudaMalloc((void**)&d_b, sizeof(T)*nd),"cudaMalloc d_b");
  CHECK_DO(cudaMalloc((void**)&d_v, sizeof(T)*nd),"cudaMalloc d_v");

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

/*
  hybrid1
  基本スレッド並列化カーネル
  1PMVを1TB(32Thread=1WARP)で計算する
  1行を1/div WARPで計算する
  <<<napprox+ndense,32>>>
  1 GEMV by 1 TB
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - a2t：a2を転置版で計算するか否か(0,1)
  - a2interchange：a2のループを入れ替えるか否か(0,1)
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  6x2x2x2x2=96通り

  何故かdivが8or16の時にa2i=1,aa=1の結果がおかしい
  div8は誤差の範囲かも知れない、div16はよりおかしい
*/
  if(kernel>=1000&&kernel<1096){
	int subkernel = kernel-1000;
	int div, a2t, a2i, aa, da;
	div = subkernel%6;
	a2t = (subkernel/6)%2;
	a2i = ((subkernel/6)/2)%2;
	aa = (((subkernel/6)/2)/2)%2;
	da = ((((subkernel/6)/2)/2)/2)%2;
	div = pow(2,div);
	char name[0xff], fname[0xff];
	snprintf(name,0xff,"hybrid1_div%d_a2t%d_a2i%d_aa%d_da%d_%s", div, a2t, a2i, aa, da, typeid(T).name());
	snprintf(fname,0xff,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	// EXEC
	hmvm_cuda_hybrid1_proxy<T>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   d_sm.napprox+d_sm.ndense, 32, d_sm.ktmax*sizeof(T),
	   v, b, nd, fname, 0,
	   div, a2t, a2i, aa, da);
	// BENCH
	hmvm_cuda_hybrid1_proxy<T>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   d_sm.napprox+d_sm.ndense, 32, d_sm.ktmax*sizeof(T),
	   v, b, nd, fname, 5,
	   div, a2t, a2i, aa, da);
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
  printf("hmvm_cuda1: end\n");
}


// ######## ######## ######## ######## ######## ######## ######## ########
// template関数の実体化のための宣言
// ######## ######## ######## ######## ######## ######## ######## ########
template void hmvm_cuda1<float>(matrix2<float>  *mat2, float *b, int kernel, int dump_result);
template void hmvm_cuda1<double>(matrix2<double> *mat2, double *b, int kernel, int dump_result);
