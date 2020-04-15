#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "hacapk_c.h"

// ######## ######## ######## ########
void hmvm_seq_1(double *v, matrix *mat, double *b)
{
  int i, j;
  double *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1: begin\n");
  printf("hmvm_seq1: nlf=%d\n", mat->nlf);
#endif
  tmp = (double*)malloc(sizeof(double)*mat->ktmax);
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
	}else{
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat->submat[i].a1[itl] * b[itt];
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
void hmvm_seq_1t(double *v, matrix *mat, double *b)
{
  int i, j;
  double *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq1_t: begin\n");
  printf("hmvm_seq1_t: nlf=%d\n", mat->nlf);
#endif
  tmp = (double*)malloc(sizeof(double)*mat->ktmax);
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
	}else{
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat->submat[i].a1[itl] * b[itt];
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
void hmvm_seq_2(double *v, matrix2 *mat, double *b)
{
  int i, j;
  double *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2: begin\n");
  printf("hmvm_seq2: nlf=%d\n", mat->nlf);
#endif
  tmp = (double*)malloc(sizeof(double)*mat->ktmax);
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
	  for(j=0;j<kt;j++)tmp[j]=0.0;
	  head = mat->a1[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  tmp[il] += mat->rowmat[head+itl] * b[itt];
		}
	  }
	  head = mat->a2[i];
	  for(il=0;il<kt;il++){
		for(it=0;it<ndl;it++){
		  ill = it+nstrtl-1;
		  itl=it+il*ndl;
		  v[ill] += mat->rowmat[head+itl] * tmp[il];
		}
	  }
	}else{
	  head = mat->a1[i];
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat->rowmat[head+itl] * b[itt];
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
void hmvm_seq_2t(double *v, matrix2 *mat, double *b)
{
  int i, j;
  double *tmp;
#if _DEBUG_LEVEL >= 1
  printf("hmvm_seq2_t: begin\n");
  printf("hmvm_seq2_t: nlf=%d\n", mat->nlf);
#endif
  tmp = (double*)malloc(sizeof(double)*mat->ktmax);
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
	}else{
	  head = mat->a1[i];
	  for(il=0;il<ndl;il++){
		ill=il+nstrtl-1;
		for(it=0;it<ndt;it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  v[ill] += mat->rowmat_t[head+itl] * b[itt];
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
void hmvm_seq(matrix *mat, matrix2 *mat2, double *b, int dump_result)
{
  int i, nd;
  FILE *F;
  double *v=NULL;
  printf("hmvm_seq: begin\n"); fflush(stdout);
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(double*)malloc(sizeof(double)*nd);

  // hmvm
  if(mat!=NULL){
	printf("hmvm_seq_1\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_1(v, mat, b);
	if(dump_result){
	  F = fopen("result_seq_1_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm (loop interchanged)
  if(mat!=NULL){
	printf("hmvm_seq_1t\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_1t(v, mat, b);
	if(dump_result){
	  F = fopen("result_seq_1t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm using rowmat array
  if(mat2!=NULL){
	printf("hmvm_seq_2\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_2(v, mat2, b);
	if(dump_result){
	  F = fopen("result_seq_2_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm using rowmat array (loop interchanged)
  if(mat2!=NULL){
	printf("hmvm_seq_2t\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_seq_2t(v, mat2, b);
	if(dump_result){
	  F = fopen("result_seq_2t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  free(v);

  printf("hmvm_seq: end\n"); fflush(stdout);
}

// ######## ######## ######## ########
void hmvm_seq_bench(matrix *mat, matrix2 *mat2, double *b)
{
  const int L=10;
  int i, l, nd;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  double *v=NULL;
  printf("hmvm_seq_bench: begin\n"); fflush(stdout);
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(double*)malloc(sizeof(double)*nd);

  // hmvm
  if(mat!=NULL){
	printf("hmvm_seq_1\n"); fflush(stdout);
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
	printf("TIME %d hmvm_seq_1 min %e max %e avg %e\n", L, dmin, dmax, davg); fflush(stdout);
  }

  // hmvm (loop interchanged)
  if(mat!=NULL){
	printf("hmvm_seq_1t\n"); fflush(stdout);
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
	printf("TIME %d hmvm_seq_1t min %e max %e avg %e\n", L, dmin, dmax, davg); fflush(stdout);
  }

  // hmvm using rowmat array
  if(mat2!=NULL){
	printf("hmvm_seq_2\n"); fflush(stdout);
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
	printf("TIME %d hmvm_seq_2 min %e max %e avg %e\n", L, dmin, dmax, davg); fflush(stdout);
  }

  // hmvm using rowmat array (loop interchanged)
  if(mat2!=NULL){
	printf("hmvm_seq_2t\n"); fflush(stdout);
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
	printf("TIME %d hmvm_seq_2t min %e max %e avg %e\n", L, dmin, dmax, davg); fflush(stdout);
  }

  free(v);

  printf("hmvm_seq_bench: end\n"); fflush(stdout);
}
