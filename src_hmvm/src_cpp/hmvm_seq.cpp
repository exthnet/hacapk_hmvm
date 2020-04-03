#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "hacapk_cpp.hpp"

// ######## ######## ######## ########
template<class T>
void hmvm_seq_1(T *v, struct matrix<T> mat, T *b)
{
  int i, j;
  T *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1: begin\n");
  printf("hmvm_seq1: nlf=%d\n", mat.nlf);
#endif
  tmp = (T*)malloc(sizeof(T)*mat.ktmax);
  for(i=0;i<mat.nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	ltmtx  = mat.submat[i].ltmtx;
	ndl    = mat.submat[i].ndl;
	ndt    = mat.submat[i].ndt;
	nstrtl = mat.submat[i].nstrtl;
	nstrtt = mat.submat[i].nstrtt;
	kt     = mat.submat[i].kt;
	if(ltmtx==1){
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat.submat[i].a1[itl] * b[itt];
		}
	  }
	  for(il=0;il<kt;il++){
		for(it=0;it<ndl;it++){
		  ill = it+nstrtl-1;
		  itl=it+il*ndl;
		  v[ill] += mat.submat[i].a2[itl] * tmp[il];
		}
	  }
	}else{
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat.submat[i].a1[itl] * b[itt];
		}
	  }
	}
  }
  free(tmp);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1: end\n");
#endif
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_1t(T *v, struct matrix<T> mat, T *b)
{
  int i, j;
  T *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1_t: begin\n");
  printf("hmvm_seq1_t: nlf=%d\n", mat.nlf);
#endif
  tmp = (T*)malloc(sizeof(T)*mat.ktmax);
  for(i=0;i<mat.nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	ltmtx  = mat.submat[i].ltmtx;
	ndl    = mat.submat[i].ndl;
	ndt    = mat.submat[i].ndt;
	nstrtl = mat.submat[i].nstrtl;
	nstrtt = mat.submat[i].nstrtt;
	kt     = mat.submat[i].kt;
	if(ltmtx==1){
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat.submat[i].a1[itl] * b[itt];
		}
	  }
	  for(il=0;il<ndl;il++){
		for(it=0;it<kt;it++){
		  ill=il+nstrtl-1;
		  itl=it+il*kt;
		  v[ill] += mat.submat[i].a2t[itl] * tmp[it];
		}
	  }
	}else{
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat.submat[i].a1[itl] * b[itt];
		}
	  }
	}
  }
  free(tmp);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1_t: end\n");
#endif
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_2(T *v, struct matrix2<T> mat, T *b)
{
  int i, j;
  T *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2: begin\n");
  printf("hmvm_seq2: nlf=%d\n", mat.nlf);
#endif
  tmp = (T*)malloc(sizeof(T)*mat.ktmax);
  for(i=0;i<mat.nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	int head;
	ltmtx  = mat.ltmtx[i];
	ndl    = mat.ndl[i];
	ndt    = mat.ndt[i];
	nstrtl = mat.nstrtl[i];
	nstrtt = mat.nstrtt[i];
	kt     = mat.kt[i];
	if(ltmtx==1){
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  head = mat.a1[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat.rowmat[head+itl] * b[itt];
		}
	  }
	  head = mat.a2[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndl;it++){
		  ill = it+nstrtl-1;
		  itl=it+il*ndl;
		  v[ill] += mat.rowmat[head+itl] * tmp[il];
		}
	  }
	}else{
	  head = mat.a1[i];
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat.rowmat[head+itl] * b[itt];
		}
	  }
	}
  }
  free(tmp);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2: end\n");
#endif
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_2t(T *v, struct matrix2<T> mat, T *b)
{
  int i, j;
  T *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2_t: begin\n");
  printf("hmvm_seq2_t: nlf=%d\n", mat.nlf);
#endif
  tmp = (T*)malloc(sizeof(T)*mat.ktmax);
  for(i=0;i<mat.nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	int il, it, ill, itt, itl;
	int head;
	ltmtx  = mat.ltmtx[i];
	ndl    = mat.ndl[i];
	ndt    = mat.ndt[i];
	nstrtl = mat.nstrtl[i];
	nstrtt = mat.nstrtt[i];
	kt     = mat.kt[i];
	if(ltmtx==1){
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  head = mat.a1[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat.rowmat_t[head+itl] * b[itt];
		}
	  }
	  head = mat.a2[i];
	  for(il=0;il<ndl;il++){
		for(it=0;it<kt;it++){
		  ill=il+nstrtl-1;
		  itl=it+il*kt;
		  v[ill] += mat.rowmat_t[head+itl] * tmp[it];
		}
	  }
	}else{
	  head = mat.a1[i];
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat.rowmat_t[head+itl] * b[itt];
		}
	  }
	}
  }
  free(tmp);
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2_t: end\n");
#endif
}


// ######## ######## ######## ########
template<class T>
void hmvm_seq(struct matrix<T> mat, struct matrix2<T> mat2, T *b)
{
  int i, nd=mat.nd;
  FILE *F;
  T *v=NULL;
  printf("hmvm_seq: begin\n");
  v=(T*)malloc(sizeof(T)*nd);

  // hmvm
  printf("hmvm_seq_1\n");
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_seq_1(v, mat, b);
  F = fopen("seq_1.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  // hmvm (loop interchanged)
  printf("hmvm_seq_1t\n");
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_seq_1t(v, mat, b);
  F = fopen("seq_1t.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  // hmvm using rowmat array
  printf("hmvm_seq_2\n");
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_seq_2(v, mat2, b);
  F = fopen("seq_2.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  // hmvm using rowmat array (loop interchanged)
  printf("hmvm_seq_2t\n");
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_seq_2t(v, mat2, b);
  F = fopen("seq_2t.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  free(v);

  printf("hmvm_seq: end\n");
}

// ######## ######## ######## ########
template<class T>
void hmvm_seq_bench(struct matrix<T> mat, struct matrix2<T> mat2, T *b)
{
  const int L=10;
  int i, l, nd=mat.nd;
  FILE *F;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  T *v=NULL;
  printf("hmvm_seq_bench: begin\n");
  v=(T*)malloc(sizeof(T)*nd);

  // hmvm
  {
	printf("hmvm_seq_1\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_1(v, mat, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	F = fopen("seq_1.txt", "w");
	for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
	fclose(F);
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
  {
	printf("hmvm_seq_1t\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_1t(v, mat, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	F = fopen("seq_1t.txt", "w");
	for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
	fclose(F);
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
  {
	printf("hmvm_seq_2\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_2(v, mat2, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	F = fopen("seq_2.txt", "w");
	for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
	fclose(F);
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
  {
	printf("hmvm_seq_2t\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_seq_2t(v, mat2, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	F = fopen("seq_2t.txt", "w");
	for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
	fclose(F);
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

  printf("hmvm_seq_bench: end\n");
}
