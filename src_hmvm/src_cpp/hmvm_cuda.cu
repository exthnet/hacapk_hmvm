// -*- C++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>
#include <cuda_runtime_api.h>

#include "hacapk.h"
#include "hmvm_cuda_kernels.cu"

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
  printf("hmvm_cuda1_%s: begin\n", typeid(T).name()); fflush(stdout);
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

#if 0
#define EXEC(FUNCNAME,BLOCKS,THREADS,S)									\
  for(i=0;i<nd;i++)v[i] = (T)0.0;										\
  CHECK_DO(cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v"); \
  CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b"); \
  cudaDeviceSynchronize();												\
  FUNCNAME<<<BLOCKS,THREADS,S>>>										\
	(d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,			\
	 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,				\
	 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
  cudaDeviceSynchronize();												\
  CHECK_DO(cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v"); \
  printf("write to %s\n", fname);										\
  F = fopen(fname, "w");												\
  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);							\
  fclose(F);

#define BENCH(FUNCNAME,BLOCKS,THREADS,S)

#define EXEC2(FUNCNAME,BLOCKS,THREADS,S,DIV,OPT)							\
  for(i=0;i<nd;i++)v[i] = (T)0.0;										\
  CHECK_DO(cudaMemcpy(d_v, v, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy v to d_v"); \
  CHECK_DO(cudaMemcpy(d_b, b, sizeof(T)*nd, cudaMemcpyHostToDevice),"cudaMemcpy b to d_b"); \
  cudaDeviceSynchronize();												\
  FUNCNAME<T,DIV,OPT><<<BLOCKS,THREADS,S>>>							\
	(d_v, d_b, d_sm.nlf, d_sm.ktmax,									\
	 d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,			\
	 d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,	\
	 d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense);				\
  cudaDeviceSynchronize();												\
  CHECK_DO(cudaMemcpy(v, d_v, sizeof(T)*nd, cudaMemcpyDeviceToHost),"cudaMemcpy d_v to v"); \
  printf("write to %s\n", fname);										\
  F = fopen(fname, "w");												\
  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);							\
  fclose(F);

#define BENCH2(FUNCNAME,BLOCKS,THREADS,S,DIV,OPT)
#endif

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
	char name[64], fname[128];
	snprintf(name,64,"seq_a2t%d_a2i%d%s", a2t, a2i, typeid(T).name());
	snprintf(fname,128,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	// EXEC
	//EXEC((hmvm_cuda_seq<T,1,1>),1,1,d_sm.ktmax*sizeof(T));
	hmvm_cuda_seq_proxy
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   v, b, nd, fname, 0,
	   a2t, a2i);
	if(0)hmvm_cuda_seq_proxy
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   v, b, nd, fname, 5,
	   a2t, a2i);
	printf("TIME %d hmvm_cuda1_seq%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
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
	char name[64], fname[128];
	snprintf(name,64,"block_a2t%d_a2i%d%s", a2t, a2i, typeid(T).name());
	snprintf(fname,128,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	// EXEC
	hmvm_cuda_block_proxy
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   v, b, nd, fname, 0,
	   a2t, a2i);
	if(0)hmvm_cuda_block_proxy
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	   d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   v, b, nd, fname, 5,
	   a2t, a2i);
	printf("TIME %d hmvm_cuda1_block%s min %e max %e avg %e\n", L-M, typeid(T).name(), dmin, dmax, davg);
  }

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
	/*
	hmvm_cuda_hybrid1_proxy
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T),
	   v, b, nd, fname, 0,
	   div, a2t, a2i, aa, da);
	// BENCH
	if(0)hmvm_cuda_hybrid1_proxy
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   d_sm.napprox+d_sm.ndense,32,d_sm.ktmax*sizeof(T),
	   v, b, nd, fname, 5,
	   div, a2t, a2i, aa, da);
	*/
  }

#if 1
/*
  hybrid2
  複数WARP単一GEMV個別行カーネル
  <<<napprox+ndense,32*mul>>>
  1 GEMV by 1 WARP
  1 line by 1/div WARP
  バリエーション
  - div：1行を1/divのWARPで計算する、div=1,2,4,8,16,32
  - mul：立ち上げるスレッド数(mul*32)、1つのmat-mat-vecまたはmat-vecをmul TBで実行、mul=1,2,3,...,16
  - a2t：a2を転置版で計算するか否か(0,1)
  - a2interchange：a2のループを入れ替えるか否か(0,1)
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか(0,1)
  6x16x2x2x2x2=1536通り
  全部大丈夫そう
*/
  if(kernel>=10000&&kernel<11536){
	int subkernel = kernel-10000;
	int div, mul, a2t, a2i, aa, da;
	div = subkernel%6;
	mul = (subkernel/6)%16 + 1;
	a2t = ((subkernel/6)/16)%2;
	a2i = (((subkernel/6)/16)/2)%2;
	aa = ((((subkernel/6)/16)/2)/2)%2;
	da = (((((subkernel/6)/16)/2)/2)/2)%2;
	div = pow(2,div);
	if((32*mul)%div)printf("invalid parameters: 32*%d %% %d\n", mul, div);
	char name[0xff], fname[0xff];
	snprintf(name,0xff,"hybrid2_div%d_mul%d_a2t%d_a2i%d_aa%d_da%d_%s", div, mul, a2t, a2i, aa, da, typeid(T).name());
	snprintf(fname,0xff,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	//printf("DIV = %d, MUL = %d, ATOMIC = %d\n", DIV, MUL, ATOMIC);
	// EXEC
	hmvm_cuda_hybrid2_proxy<T>
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   d_sm.napprox+d_sm.ndense, 32*mul, d_sm.ktmax*sizeof(T),
	   v, b, nd, fname, 0,
	   div, mul, a2t, a2i, aa, da);
	// BENCH
	/*
	if(0)hmvm_cuda_hybrid2_proxy
	  (d_v, d_b, d_sm.nlf, d_sm.ktmax,
	   d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt, d_sm.kt,
	   d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	   d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	   d_sm.napprox+d_sm.ndense,32*mul,d_sm.ktmax*sizeof(T),
	   v, b, nd, fname, 5,
	   div, mul, a2t, a2i, aa, da);
	*/
  }
#endif

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


#if 1
/*
  hybrid4
  複数WARP単一GEMV単一行GEMVカーネル
  <<<napprox+ndense, 32*mul>>>
  1 GEMV by mul TB
  1 line by 1/mul WARP
  バリエーション
  - mul：
  - a2t：a2を転置版で計算するか否か
  - a2interchange：a2のループを入れ替えるか否か
  - aatomic：approxの計算をatomic優先にするかwarp shuffle併用するか
  - datomic：denseの計算をatomic優先にするかwarp shuffle併用するか
*/
  if((kernel>=4000)&&(kernel<4192)){
	int subkernel = kernel-4000;
	int MUL, ATOMIC;
	ATOMIC = subkernel/16;
	MUL = subkernel%16;
	ATOMIC = 1; MUL = 1; // test
	char name[32], fname[64];
	snprintf(name,32,"hybrid4_mul%d_atomic%d_%s", MUL, ATOMIC, typeid(T).name());
	snprintf(fname,64,"result_cuda1_%s.txt", name);
	printf("fname = %s\n", fname);
	printf("MUL = %d, ATOMIC = %d\n", MUL, ATOMIC);
	// EXEC
	{
    hmvm_cuda_hybrid4_proxy<T>(d_v, d_b, d_sm.nlf, d_sm.ktmax,
	d_sm.ltmtx, d_sm.ndt, d_sm.ndl, d_sm.nstrtl, d_sm.nstrtt,
	d_sm.kt, d_sm.a1, d_sm.a2, d_sm.rowmat, d_sm.rowmat_t,
	d_sm.napprox, d_sm.approx, d_sm.ndense, d_sm.dense,
	d_sm.napprox+d_sm.ndense,32*MUL,d_sm.ktmax*sizeof(T)*MUL,MUL,ATOMIC, v, b, nd, fname, 0);
	//(d_sm.napprox+MUL-1)/MUL,32*MUL,d_sm.ktmax*sizeof(T)*MUL,DIV,MUL,ATOMIC, v, b, nd, fname, 0);
	//(d_sm.ndense+MUL-1)/MUL,32*MUL,d_sm.ktmax*sizeof(T),DIV,MUL,ATOMIC, v, b, nd, fname, 0);
	}
	// BENCH
	{
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
