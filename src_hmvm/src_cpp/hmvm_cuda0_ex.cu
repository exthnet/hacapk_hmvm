// -*- C++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>
#include <cuda_runtime_api.h>

#include "hacapk.h"

#define CHECK_DO(act,msg) {ret=act; if(ret!=cudaSuccess){printf("%s failed\n",msg);exit(-1);};}

//namespace cg = cooperative_groups;

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
  T tmp1;
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
	if(a2t==0){
	  for(il=0; il<kt; il++){
		tmp1 = (T)0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += rowmat[head+itl]*d_zu[itt];
		}
		tmp2[il] = tmp1;
	  }
	}else{
	  for(il=0; il<kt; il++){
		tmp1 = (T)0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += rowmat_t[head+itl]*d_zu[itt];
		}
		tmp2[il] = tmp1;
	  }
	}
	head = a2[ip];
	if(a2t==0){
	  if(a2i==0){
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			//myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
			d_zaut[ill] += rowmat[head+itl]*tmp2[il];
		  }
		}
	  }else{ // a2i==1
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(il=0; il<kt; il++){
			itl=it+il*ndl;
			//myAtomicAdd(&d_zaut[ill], rowmat[head+itl]*tmp2[il]);
			tmp1 += rowmat[head+itl]*tmp2[il];
		  }
		  d_zaut[ill] += tmp1;
		}
	  }
	}else{ // a2t==1
	  if(a2i==0){
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it*kt+il;
			//myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
			d_zaut[ill] += rowmat_t[head+itl]*tmp2[il];
		  }
		}
	  }else{ // a2i==1
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(il=0; il<kt; il++){
			itl=it*kt+il;
			//myAtomicAdd(&d_zaut[ill], rowmat_t[head+itl]*tmp2[il]);
			tmp1 += rowmat_t[head+itl]*tmp2[il];
		  }
		  d_zaut[ill] += tmp1;
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
	  //myAtomicAdd(&d_zaut[ill], tmp);
	  d_zaut[ill] += tmp;
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cuda_seq : end\n");
#endif
}

template <class T, int a2t, int a2i>
void hmvm_cuda_seq_proxy
(T *d_zaut, T *d_zu,
 matrix2<T> *h_mat, matrix2<T> *d_mat,
 T *v, T *b, char *fname, int bench)
{
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<h_mat->nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*h_mat->nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*h_mat->nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	hmvm_cuda_seq<T,a2t,a2i><<<1,1,sizeof(T)*h_mat->ktmax>>>
	  (d_zaut, d_zu, d_mat->nlf, d_mat->ktmax,
	   d_mat->ltmtx, d_mat->ndt, d_mat->ndl, d_mat->nstrtl, d_mat->nstrtt,
	   d_mat->kt, d_mat->a1, d_mat->a2, d_mat->rowmat, d_mat->rowmat_t,
	   d_mat->napprox, d_mat->approx, d_mat->ndense, d_mat->dense);
	cudaDeviceSynchronize();
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }
  if(bench==0){
	CHECK_DO(cudaMemcpy(v, d_zaut, sizeof(T)*h_mat->nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v");
	printf("write to %s\n", fname);
	F = fopen(fname, "w");
	for(i=0;i<h_mat->nd;i++)fprintf(F, "%.3E\n", v[i]);
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
	printf("TIME %d hmvm_cuda0_seq%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
}

template <class T>
void hmvm_cuda_seq_proxy
(T *d_zaut, T *d_zu,
 matrix2<T> *h_mat, matrix2<T> *d_mat,
 T *v, T *b, char *fname, int bench,
 int a2t, int a2i)
{
  if(a2t==0 && a2i==0)hmvm_cuda_seq_proxy<T,0,0>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
  if(a2t==0 && a2i==1)hmvm_cuda_seq_proxy<T,0,1>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
  if(a2t==1 && a2i==0)hmvm_cuda_seq_proxy<T,1,0>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
  if(a2t==1 && a2i==1)hmvm_cuda_seq_proxy<T,1,1>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
}

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
  T tmp1;
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
	if(a2t==0){
	  for(il=0; il<kt; il++){
		tmp1 = (T)0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += rowmat[head+itl]*d_zu[itt];
		}
		tmp2[il] = tmp1;
	  }
	}else{
	  for(il=0; il<kt; il++){
		tmp1 = (T)0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += rowmat_t[head+itl]*d_zu[itt];
		}
		tmp2[il] = tmp1;
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
		  tmp1 = (T)0.0;
		  for(il=0; il<kt; il++){
			itl=it+il*ndl;
			tmp1 += rowmat[head+itl]*tmp2[il];
		  }
		  myAtomicAdd(&d_zaut[ill], tmp1);
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
		  tmp1 = (T)0.0;
		  for(il=0; il<kt; il++){
			itl=it*kt+il;
			tmp1 += rowmat_t[head+itl]*tmp2[il];
		  }
		  myAtomicAdd(&d_zaut[ill], tmp1);
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
	  tmp1 = (T)0.0;
	  ill=il+nstrtl-1;
	  for(it=0; it<ndt; it++){
		itt=it+nstrtt-1;
		itl=it+il*ndt;
		if(a2t==0){
		  tmp1 += rowmat[head+itl]*d_zu[itt];
		}else{
		  tmp1 += rowmat_t[head+itl]*d_zu[itt];
		}
	  }
	  myAtomicAdd(&d_zaut[ill], tmp1);
	}
#endif
  }
#if _DEBUG_LEVEL >= 2
  printf("hmvm_cudaD_block : end\n");
#endif
}

template <class T, int a2t, int a2i>
void hmvm_cuda_block_proxy
(T *d_zaut, T *d_zu,
 matrix2<T> *h_mat, matrix2<T> *d_mat,
 T *v, T *b, char *fname, int bench)
{
  int M=5, L=M+bench;
  FILE *F;
  int i, l, lmax;
  double d1, d2, *dtimes, dmin, dmax, davg;
  cudaError_t ret;
  dtimes = new double[L];
  if(bench==0){lmax=1;}else{lmax=L;}
  for(l=0;l<lmax;l++){
	for(i=0;i<h_mat->nd;i++)v[i] = (T)0.0;
	CHECK_DO(cudaMemcpy(d_zaut, v, sizeof(T)*h_mat->nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v");
	CHECK_DO(cudaMemcpy(d_zu, b, sizeof(T)*h_mat->nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b");
	cudaDeviceSynchronize();
	d1 = omp_get_wtime();
	hmvm_cuda_block<T,a2t,a2i><<<h_mat->napprox+h_mat->ndense,1,sizeof(T)*h_mat->ktmax>>>
	  (d_zaut, d_zu, d_mat->nlf, d_mat->ktmax,
	   d_mat->ltmtx, d_mat->ndt, d_mat->ndl, d_mat->nstrtl, d_mat->nstrtt,
	   d_mat->kt, d_mat->a1, d_mat->a2, d_mat->rowmat, d_mat->rowmat_t,
	   d_mat->napprox, d_mat->approx, d_mat->ndense, d_mat->dense);
	cudaDeviceSynchronize();
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }
  if(bench==0){
	CHECK_DO(cudaMemcpy(v, d_zaut, sizeof(T)*h_mat->nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v");
	printf("write to %s\n", fname);
	F = fopen(fname, "w");
	for(i=0;i<h_mat->nd;i++)fprintf(F, "%.3E\n", v[i]);
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
	printf("TIME %d hmvm_cuda0_block%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }
  delete [] dtimes;
}

template <class T>
void hmvm_cuda_block_proxy
(T *d_zaut, T *d_zu,
 matrix2<T> *h_mat, matrix2<T> *d_mat,
 T *v, T *b, char *fname, int bench,
 int a2t, int a2i)
{
  if(a2t==0 && a2i==0)hmvm_cuda_block_proxy<T,0,0>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
  if(a2t==0 && a2i==1)hmvm_cuda_block_proxy<T,0,1>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
  if(a2t==1 && a2i==0)hmvm_cuda_block_proxy<T,1,0>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
  if(a2t==1 && a2i==1)hmvm_cuda_block_proxy<T,1,1>(d_zaut, d_zu, h_mat, d_mat, v, b, fname, bench);
}

// ######## ######## ######## ######## ######## ######## ######## ########
template<class T>
void hmvm_cuda0(matrix2<T> *mat2, T *b, int kernel, int dump_result, int nbench)
{
  matrix2<T> d_sm;
  int i, nd = mat2->nd, ktmax = mat2->ktmax, nlf = mat2->nlf;
  T *v=NULL;
  T *d_b, *d_v;
  int ip;
  int len;
  cudaError_t ret;
  printf("hmvm_cuda0_%s: begin\n", typeid(T).name()); fflush(stdout);
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
	完全逐次
	バリエーション：
	- a2trans：approxy2を転置版で計算するか否か(0,1)
	- a2interchange：approxy2のループを入れ替えるか否か(0,1)
  */
  if(kernel>=0 && kernel<4)
  {
	int subkernel = kernel;
	int a2t, a2i;
	a2t = subkernel%2;
	a2i = (subkernel/2)%2;
	char name[0xff], fname[0xff];
	snprintf(name,0xff,"seq_a2t%d_a2i%d%s", a2t, a2i, typeid(T).name());
	snprintf(fname,0xff,"result_cuda0_%s.txt", name);
	printf("fname = %s\n", fname);
	// EXEC
	if(dump_result)hmvm_cuda_seq_proxy<T>(d_v, d_b, mat2, &d_sm, v, b, fname, 0, a2t, a2i);
	// BENCH
	if(nbench>0)hmvm_cuda_seq_proxy<T>(d_v, d_b, mat2, &d_sm, v, b, fname, nbench, a2t, a2i);
  }

  /*
	block並列化
	ThreadBlockごとに1つの部分行列積(mat-mat-vecまたはmat-vec)を行う
	ThreadBlock内部は逐次
	- a2trans：approxy2を転置版で計算するか否か(0,1)
	- a2interchange：approxy2のループを入れ替えるか否か(0,1)
  */
  if(kernel>=10 && kernel<14)
  {
	int subkernel = kernel-10;
	int a2t, a2i;
	a2t = subkernel%2;
	a2i = (subkernel/2)%2;
	char name[0xff], fname[0xff];
	snprintf(name,0xff,"block_a2t%d_a2i%d%s", a2t, a2i, typeid(T).name());
	snprintf(fname,0xff,"result_cuda0_%s.txt", name);
	printf("fname = %s\n", fname);
	// EXEC
	if(dump_result)hmvm_cuda_block_proxy<T>(d_v, d_b, mat2, &d_sm, v, b, fname, 0, a2t, a2i);
	// BENCH
	if(nbench>0)hmvm_cuda_block_proxy<T>(d_v, d_b, mat2, &d_sm, v, b, fname, nbench, a2t, a2i);
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
  printf("hmvm_cuda0: end\n");
}

// ######## ######## ######## ######## ######## ######## ######## ########
// template関数の実体化のための宣言
// ######## ######## ######## ######## ######## ######## ######## ########
template void hmvm_cuda0<float>(matrix2<float>  *mat2, float *b, int kernel, int dump_result, int nbench);
template void hmvm_cuda0<double>(matrix2<double> *mat2, double *b, int kernel, int dump_result, int nbench);
