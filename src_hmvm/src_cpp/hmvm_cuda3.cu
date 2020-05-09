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
  hybrid3
  複数WARP個別GEMVカーネル
  1PMVを1WARPが担当
  1TBあたりスレッド数は32*mul(mul WARP)
  PMV内の1行を1/div WARPが担当
  hybrid1->hybrid2では増やしたWARPを同一PMVの計算に割り当てたが、hybrid3ではそれぞれ個別のPMVを計算する
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

#if 0
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
}
#include "template_hybrid3.hpp"
#endif
// ######## ######## ######## ######## ######## ######## ######## ########

template<class T>
void hmvm_cuda3(matrix2<T> *mat2, T *b, int kernel, int dump_result)
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
  printf("hmvm_cuda3_%s: begin\n", typeid(T).name()); fflush(stdout);
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

#if 1
/*
  hybrid3
  複数WARP個別GEMVカーネル
  1PMVを1WARPが担当
  1TBあたりスレッド数は32*mul(mul WARP)
  PMV内の1行を1/div WARPが担当
  <<<napprox/mul+ndense/mul, 32*mul>>>
  1 GEMV by 1 WARP
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - mul：立ち上げるブロック数スレッド数の係数、1つのmat-mat-vecまたはmat-vecを1WARPで実行、mul=1,2,3,...,
  - a2t：a2を転置版で計算するか否か
  - a2interchange：a2のループを入れ替えるか否か
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか
*/
  if((kernel>=20000)&&(kernel<21536)){
	int subkernel = kernel-20000;
	int div, mul, a2t, a2i, aa, da;
	div = subkernel%6;
	mul = (subkernel/6)%16 + 1;
	a2t = ((subkernel/6)/16)%2;
	a2i = (((subkernel/6)/16)/2)%2;
	aa = ((((subkernel/6)/16)/2)/2)%2;
	da = (((((subkernel/6)/16)/2)/2)/2)%2;
	div = pow(2,div);
	char name[0xff], fname[0xff];
	snprintf(name,0xff,"hybrid3_div%d_mul%d_a2t%d_a2i%d_aa%d_da%d_%s", div, mul, a2t, a2i, aa, da, typeid(T).name());
	snprintf(fname,0xff,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	//printf("DIV = %d, MUL = %d, ATOMIC = %d\n", DIV, MUL, ATOMIC);
	// EXEC
#if 1
	hmvm_cuda_hybrid3_proxy<T>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   (d_sm.napprox+mul-1)/mul+(d_sm.ndense+mul-1)/mul, 32*mul, d_sm.ktmax*sizeof(T)*mul,
	   v, b, nd, fname, 0,
	   div, mul, a2t, a2i, aa, da);
	  //(d_sm.napprox+MUL-1)/MUL,32*MUL,d_sm.ktmax*sizeof(T)*MUL,DIV,MUL,ATOMIC, v, b, nd, fname, 0);
	  //(d_sm.ndense+MUL-1)/MUL,32*MUL,d_sm.ktmax*sizeof(T),DIV,MUL,ATOMIC, v, b, nd, fname, 0);
#endif
	// BENCH
	/*
	if(0){
	  hmvm_cuda_hybrid3_proxy<T>(d_v, d_b, d_sm.nlf, d_sm.ktmax,
								 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
								 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
								 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
								 (d_sm.napprox+MUL-1)/MUL+(d_sm.ndense+MUL-1)/MUL,32*MUL,d_sm.ktmax*sizeof(T),DIV,MUL,ATOMIC, v, b, nd, fname, 1);
	}
	*/
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
  printf("hmvm_cuda3: end\n");
}


// ######## ######## ######## ######## ######## ######## ######## ########
// template関数の実体化のための宣言
// ######## ######## ######## ######## ######## ######## ######## ########
template void hmvm_cuda3<float>(matrix2<float>  *mat2, float *b, int kernel, int dump_result);
template void hmvm_cuda3<double>(matrix2<double> *mat2, double *b, int kernel, int dump_result);