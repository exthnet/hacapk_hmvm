#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>

#include "hacapk.h"

// ######## ######## ######## ########
template<class T>
void hmvm_seq_1(T *v, matrix<T> *mat, T *b)
{
  int i, j;
  T *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1: begin\n");
  printf("hmvm_seq1: nlf=%d\n", mat->nlf);
#endif
  tmp = (T*)malloc(sizeof(T)*mat->ktmax);
  for(i=0;i<mat->nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	ltmtx  = mat->submat[i].ltmtx;
	ndl    = mat->submat[i].ndl;
	ndt    = mat->submat[i].ndt;
	nstrtl = mat->submat[i].nstrtl;
	nstrtt = mat->submat[i].nstrtt;
	kt     = mat->submat[i].kt;
	if(ltmtx==1){
#ifndef _SKIP_APPROX
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat->submat[i].a1[itl] * b[itt];
		}
	  }
	  for(il=0;il<kt;il++){
		for(it=0;it<ndl;it++){
		  ill = it+nstrtl-1;
		  itl=it+il*ndl;
		  v[ill] += mat->submat[i].a2[itl] * tmp[il];
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat->submat[i].a1[itl] * b[itt];
		}
	  }
#endif
	}
  }
  free(tmp);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1: end\n");
#endif
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_1t(T *v, matrix<T> *mat, T *b)
{
  int i, j;
  T *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1_t: begin\n");
  printf("hmvm_seq1_t: nlf=%d\n", mat->nlf);
#endif
  tmp = (T*)malloc(sizeof(T)*mat->ktmax);
  for(i=0;i<mat->nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	ltmtx  = mat->submat[i].ltmtx;
	ndl    = mat->submat[i].ndl;
	ndt    = mat->submat[i].ndt;
	nstrtl = mat->submat[i].nstrtl;
	nstrtt = mat->submat[i].nstrtt;
	kt     = mat->submat[i].kt;
	if(ltmtx==1){
#ifndef _SKIP_APPROX
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat->submat[i].a1[itl] * b[itt];
		}
	  }
	  for(il=0;il<ndl;il++){
		for(it=0;it<kt;it++){
		  ill=il+nstrtl-1;
		  itl=it+il*kt;
		  v[ill] += mat->submat[i].a2t[itl] * tmp[it];
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat->submat[i].a1[itl] * b[itt];
		}
	  }
#endif
	}
  }
  free(tmp);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1_t: end\n");
#endif
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_2(T *v, matrix2<T> *mat, T *b)
{
  int i, j;
  T tmp;
  T *tmp2;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2: begin\n");
  printf("hmvm_seq2: nlf=%d\n", mat->nlf);
#endif
  tmp2 = (T*)malloc(sizeof(T)*mat->ktmax);
  for(i=0;i<mat->nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	int head;
	ltmtx  = mat->ltmtx[i];
	ndl    = mat->ndl[i];
	ndt    = mat->ndt[i];
	nstrtl = mat->nstrtl[i];
	nstrtt = mat->nstrtt[i];
	kt     = mat->kt[i];
	if(ltmtx==1){
#ifndef _SKIP_APPROX
	  for(j=0;j<kt;j++)tmp2[j]=0.0;
	  head = mat->a1[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp2[il] += mat->rowmat[head+itl] * b[itt];
		}
	  }
	  head = mat->a2[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndl;it++){
		  ill = it+nstrtl-1;
		  itl=it+il*ndl;
		  v[ill] += mat->rowmat[head+itl] * tmp2[il];
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  head = mat->a1[i];
	  for(il=0;il<ndl;il++){
		tmp = 0.0;
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp += mat->rowmat[head+itl] * b[itt];
		}
		v[ill] += tmp;
	  }
#endif
	}
  }
  free(tmp2);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2: end\n");
#endif
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_2t(T *v, matrix2<T> *mat, T *b)
{
  int i, j;
  T *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2_t: begin\n");
  printf("hmvm_seq2_t: nlf=%d\n", mat->nlf);
#endif
  tmp = (T*)malloc(sizeof(T)*mat->ktmax);
  for(i=0;i<mat->nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	int head;
	ltmtx  = mat->ltmtx[i];
	ndl    = mat->ndl[i];
	ndt    = mat->ndt[i];
	nstrtl = mat->nstrtl[i];
	nstrtt = mat->nstrtt[i];
	kt     = mat->kt[i];
	if(ltmtx==1){
#ifndef _SKIP_APPROX
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  head = mat->a1[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat->rowmat_t[head+itl] * b[itt];
		}
	  }
	  head = mat->a2[i];
	  for(il=0;il<ndl;il++){
		for(it=0;it<kt;it++){
		  ill=il+nstrtl-1;
		  itl=it+il*kt;
		  v[ill] += mat->rowmat_t[head+itl] * tmp[it];
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  head = mat->a1[i];
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat->rowmat_t[head+itl] * b[itt];
		}
	  }
#endif
	}
  }
  free(tmp);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2_t: end\n");
#endif
}


// ######## ######## ######## ########
template<class T>
void hmvm_seq(matrix<T> *mat, matrix2<T> *mat2, T *b, int dump_result)
{
  int i, nd;
  char fname[0xff];
  FILE *F;
  T *v=NULL;
  printf("hmvm_seq_%s: begin\n", typeid(T).name());
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(T*)malloc(sizeof(T)*nd);

  // hmvm
  if(mat!=NULL){
	printf("hmvm_seq_1\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_1(v, mat, b);
	if(dump_result){
	  snprintf(fname, 0xff, "result_seq_1_%s.txt", typeid(T).name());
	  F = fopen(fname, "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm (loop interchanged)
  if(mat!=NULL){
	printf("hmvm_seq_1t\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_1t(v, mat, b);
	if(dump_result){
	  snprintf(fname, 0xff, "result_seq_1t_%s.txt", typeid(T).name());
	  F = fopen(fname, "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm using rowmat array
  if(mat2!=NULL){
	printf("hmvm_seq_2\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_2(v, mat2, b);
	if(dump_result){
	  snprintf(fname, 0xff, "result_seq_2_%s.txt", typeid(T).name());
	  F = fopen(fname, "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm using rowmat array (loop interchanged)
  if(mat2!=NULL){
	printf("hmvm_seq_2t\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_2t(v, mat2, b);
	if(dump_result){
	  snprintf(fname, 0xff, "result_seq_2t_%s.txt", typeid(T).name());
	  F = fopen(fname, "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  free(v);

  printf("hmvm_seq_%s: end\n", typeid(T).name());
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_bench(matrix<T> *mat, matrix2<T> *mat2, T *b)
{
  const int L=10;
  int i, l, nd;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  T *v=NULL;
  printf("hmvm_seq_%s_bench: begin\n", typeid(T).name());
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(T*)malloc(sizeof(T)*nd);

  // hmvm
  if(mat!=NULL){
	printf("hmvm_seq_1\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_1(v, mat, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  davg += dtimes[i];
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_seq_1 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // hmvm (loop interchanged)
  if(mat!=NULL){
	printf("hmvm_seq_1t\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_1t(v, mat, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  davg += dtimes[i];
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_seq_1t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // hmvm using rowmat array
  if(mat2!=NULL){
	printf("hmvm_seq_2\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_2(v, mat2, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  davg += dtimes[i];
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_seq_2 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // hmvm using rowmat array (loop interchanged)
  if(mat2!=NULL){
	printf("hmvm_seq_2t\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_2t(v, mat2, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  davg += dtimes[i];
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_seq_2t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  free(v);

  printf("hmvm_seq_%s_bench: end\n", typeid(T).name());
}


// ######## ######## ######## ########
template void hmvm_seq<float> (matrix<float>  *mat, matrix2<float>  *mat2, float  *b, int dump_result);
template void hmvm_seq<double>(matrix<double> *mat, matrix2<double> *mat2, double *b, int dump_result);
template void hmvm_seq_bench<float> (matrix<float>  *mat, matrix2<float>  *mat2, float  *b);
template void hmvm_seq_bench<double>(matrix<double> *mat, matrix2<double> *mat2, double *b);
