#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>

#include "hacapk.h"

// ######## ######## ######## ########
template<class T, int a2i>
void hmvm_seq_1(T *v, const matrix<T> *mat, const T *b)
{
  int i, j;
  T tmp1;
  T *tmp2;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1: begin\n");
  printf("hmvm_seq1: nlf=%d\n", mat->nlf);
#endif
  tmp2 = (T*)malloc(sizeof(T)*mat->ktmax);
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
	  for(j=0;j<kt;j++)tmp2[j]=0.0;
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp2[il] += mat->submat[i].a1[itl] * b[itt];
		}
	  }
	  if(a2i==0){
		for(il=0;il<kt;il++){
		  for(it=0;it<ndl;it++){
			ill = it+nstrtl-1;
			itl=it+il*ndl;
			v[ill] += mat->submat[i].a2[itl] * tmp2[il];
		  }
		}
	  }else{
		for(it=0;it<ndl;it++){
		  ill = it+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(il=0;il<kt;il++){
			itl=it+il*ndl;
			tmp1 += mat->submat[i].a2[itl] * tmp2[il];
		  }
		  v[ill] += tmp1;
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		tmp1 = 0.0;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += mat->submat[i].a1[itl] * b[itt];
		}
		v[ill] += tmp1;
	  }
#endif
	}
  }
  free(tmp2);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1: end\n");
#endif
}

// ######## ######## ######## ########
template<class T, int a2i>
void hmvm_seq_1t(T *v, const matrix<T> *mat, const T *b)
{
  int i, j;
  T tmp1;
  T *tmp2;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1_t: begin\n");
  printf("hmvm_seq1_t: nlf=%d\n", mat->nlf);
#endif
  tmp2 = (T*)malloc(sizeof(T)*mat->ktmax);
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
	  for(j=0;j<kt;j++)tmp2[j]=0.0;
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp2[il] += mat->submat[i].a1[itl] * b[itt];
		}
	  }
	  if(a2i==0){
		for(it=0;it<kt;it++){
		  for(il=0;il<ndl;il++){
			ill=il+nstrtl-1;
			itl=it+il*kt;
			v[ill] += mat->submat[i].a2t[itl] * tmp2[it];
		  }
		}
	  }else{
		for(il=0;il<ndl;il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0;it<kt;it++){
			itl=it+il*kt;
			tmp1 += mat->submat[i].a2t[itl] * tmp2[it];
		  }
		  v[ill] += tmp1;
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		tmp1 = (T)0.0;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += mat->submat[i].a1[itl] * b[itt];
		}
		v[ill] += tmp1;
	  }
#endif
	}
  }
  free(tmp2);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1_t: end\n");
#endif
}

// ######## ######## ######## ########
template<class T, int a2i>
void hmvm_seq_2(T *v, const matrix2<T> *mat, const T *b)
{
  int i, j;
  T tmp1;
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
	  if(a2i==0){
		for(il=0;il<kt;il++){
		  for(it=0;it<ndl;it++){
			ill = it+nstrtl-1;
			itl=it+il*ndl;
			v[ill] += mat->rowmat[head+itl] * tmp2[il];
		  }
		}
	  }else{
		for(it=0;it<ndl;it++){
		  ill = it+nstrtl-1;
		  tmp1 = 0.0;
		  for(il=0;il<kt;il++){
			itl=it+il*ndl;
			tmp1 += mat->rowmat[head+itl] * tmp2[il];
		  }
		  v[ill] += tmp1;
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  head = mat->a1[i];
	  for(il=0;il<ndl;il++){
		tmp1 = (T)0.0;
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += mat->rowmat[head+itl] * b[itt];
		}
		v[ill] += tmp1;
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
template<class T, int a2i>
void hmvm_seq_2t(T *v, const matrix2<T> *mat, const T *b)
{
  int i, j;
  T tmp1;
  T *tmp2;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2_t: begin\n");
  printf("hmvm_seq2_t: nlf=%d\n", mat->nlf);
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
		  tmp2[il] += mat->rowmat_t[head+itl] * b[itt];
		}
	  }
	  head = mat->a2[i];
	  if(a2i==0){
		for(it=0;it<kt;it++){
		  for(il=0;il<ndl;il++){
			ill=il+nstrtl-1;
			itl=it+il*kt;
			v[ill] += mat->rowmat_t[head+itl] * tmp2[it];
		  }
		}
	  }else{
		for(il=0;il<ndl;il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0;it<kt;it++){
			itl=it+il*kt;
			tmp1 += mat->rowmat_t[head+itl] * tmp2[it];
		  }
		  v[ill] += tmp1;
		}
	  }
#endif
	}else{
#ifndef _SKIP_DENSE
	  head = mat->a1[i];
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		tmp1 = (T)0.0;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp1 += mat->rowmat_t[head+itl] * b[itt];
		}
		v[ill] += tmp1;
	  }
#endif
	}
  }
  free(tmp2);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2_t: end\n");
#endif
}


// ######## ######## ######## ########
// hmvm for mat1
template<class T>
void hmvm_seq_proxy1(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench,
					 void (*fn)(T *, const matrix<T> *, const T *), const char *subname)
{
  int M=5, L, lmax;
  int i, l, nd;
  double d1, d2, *dtimes, dmin, dmax, davg;
  char fname[0xff];
  FILE *F;
  T *v=NULL;
  printf("hmvm_seq_%s: begin\n", typeid(T).name());
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(T*)malloc(sizeof(T)*nd);
  L = M + nbench;
  dtimes = new double[L];
  if(nbench==0){lmax=1;}else{lmax=L;}

  // hmvm
  {
	printf("hmvm_%s\n", subname);
	for(l=0;l<lmax;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  fn(v, mat, b); // hmvm kernel
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	if(dump_result){
	  snprintf(fname, 0xff, "result_%s_%s.txt", subname, typeid(T).name());
	  F = fopen(fname, "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}else{
	  dmin = 9999.99;	dmax = 0.0;	davg = 0.0;
	  for(i=M;i<L;i++){
		davg += dtimes[i];
		if(dmin>dtimes[i])dmin=dtimes[i];
		if(dmax<dtimes[i])dmax=dtimes[i];
	  }
	  davg /= (L-M);
	  printf("TIME hmvm_%s_%s %d times min %e max %e avg %e\n", subname, typeid(T).name(), L, dmin, dmax, davg);
	}
  }
  free(v);
}

// hmvm for mat2
template<class T>
void hmvm_seq_proxy2(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench,
					 void (*fn)(T *, const matrix2<T> *, const T *), const char *subname)
{
  int M=5, L, lmax;
  int i, l, nd;
  double d1, d2, *dtimes, dmin, dmax, davg;
  char fname[0xff];
  FILE *F;
  T *v=NULL;
  printf("hmvm_seq_%s: begin\n", typeid(T).name());
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(T*)malloc(sizeof(T)*nd);
  L = M + nbench;
  dtimes = new double[L];
  if(nbench==0){lmax=1;}else{lmax=L;}

  // hmvm
  {
	printf("hmvm_%s\n", subname);
	for(l=0;l<lmax;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  fn(v, mat2, b); // hmvm_kernel
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	if(dump_result){
	  snprintf(fname, 0xff, "result_%s_%s.txt", subname, typeid(T).name());
	  F = fopen(fname, "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}else{
	  dmin = 9999.99;	dmax = 0.0;	davg = 0.0;
	  for(i=M;i<L;i++){
		davg += dtimes[i];
		if(dmin>dtimes[i])dmin=dtimes[i];
		if(dmax<dtimes[i])dmax=dtimes[i];
	  }
	  davg /= (L-M);
	  printf("TIME hmvm_%s_%s %d times min %e max %e avg %e\n", subname, typeid(T).name(), L, dmin, dmax, davg);
	}
  }
  free(v);
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_proxy(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench)
{
  // hmvm
  if(mat!=NULL){	hmvm_seq_proxy1(mat, mat2, b, 1, 0, &hmvm_seq_1<T,0>, "seq_1");  }

  // hmvm (loop interchanged)
  if(mat!=NULL){	hmvm_seq_proxy1(mat, mat2, b, 1, 0, &hmvm_seq_1<T,1>, "seq_1i");  }

  // hmvm (trans)
  if(mat!=NULL){	hmvm_seq_proxy1(mat, mat2, b, 1, 0, &hmvm_seq_1t<T,0>, "seq_1t");  }

  // hmvm (trans, loop interchanged)
  if(mat!=NULL){	hmvm_seq_proxy1(mat, mat2, b, 1, 0, &hmvm_seq_1t<T,1>, "seq_1ti");  }

  // hmvm using rowmat array
  if(mat2!=NULL){	hmvm_seq_proxy2(mat, mat2, b, 1, 0, &hmvm_seq_2<T,0>, "seq_2");  }

  // hmvm using rowmat array (loop interchanged)
  if(mat2!=NULL){	hmvm_seq_proxy2(mat, mat2, b, 1, 0, &hmvm_seq_2<T,1>, "seq_2i");  }

  // hmvm using rowmat array (trans)
  if(mat2!=NULL){	hmvm_seq_proxy2(mat, mat2, b, 1, 0, &hmvm_seq_2t<T,0>, "seq_2t");  }

  // hmvm using rowmat array (trans, loop interchanged)
  if(mat2!=NULL){	hmvm_seq_proxy2(mat, mat2, b, 1, 0, &hmvm_seq_2t<T,1>, "seq_2ti");  }
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench)
{
  printf("hmvm_seq_%s: begin\n", typeid(T).name());
  if(dump_result)hmvm_seq_proxy(mat, mat2, b, 1, 0);
  if(nbench>0)hmvm_seq_proxy(mat, mat2, b, 0, nbench);
  printf("hmvm_seq_%s: end\n", typeid(T).name());
}

// ######## ######## ######## ########
template void hmvm_seq<float> (const matrix<float>  *mat, const matrix2<float>  *mat2, const float  *b, int dump_result, int nbench);
template void hmvm_seq<double>(const matrix<double> *mat, const matrix2<double> *mat2, const double *b, int dump_result, int nbench);
