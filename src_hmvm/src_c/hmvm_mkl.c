#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mkl.h>

#include "hacapk.h"

// ######## ######## ######## ########
// mkl blas, parallel (multi-threaded blas)
// ######## ######## ######## ########
void  hmvm_blas_p_calc_1
(double *zau, matrix *mat, double *zu, double *zbut)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;
  int zero = 0;
  int one = 1;
  double dzero = 0.0;
  double done = 1.0;

  nlf=mat->nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat->submat[ip].ndl;
    ndt    = mat->submat[ip].ndt;
    nstrtl = mat->submat[ip].nstrtl;
    nstrtt = mat->submat[ip].nstrtt;

    if(mat->submat[ip].ltmtx==1){
      kt=mat->submat[ip].kt;
	  for(il=0;il<kt;il++)zbut[il]=0.0;
	  dgemv_("t", &ndt, &kt, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
	  dgemv_("n", &ndl, &kt, &done, mat->submat[ip].a2, &ndl, zbut, &one, &done, &zau[nstrtl-1], &one);
    } else if(mat->submat[ip].ltmtx==2){
	  dgemv_("t", &ndt, &ndl, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, &zau[nstrtl-1], &one);
    }
  }
}

void  hmvm_blas_p_calc_1t
(double *zau, matrix *mat, double *zu, double *zbut)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;
  int zero = 0;
  int one = 1;
  double dzero = 0.0;
  double done = 1.0;

  nlf=mat->nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat->submat[ip].ndl;
    ndt    = mat->submat[ip].ndt;
    nstrtl = mat->submat[ip].nstrtl;
    nstrtt = mat->submat[ip].nstrtt;

    if(mat->submat[ip].ltmtx==1){
      kt=mat->submat[ip].kt;
	  for(il=0;il<kt;il++)zbut[il]=0.0;
	  dgemv_("t", &ndt, &kt, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
	  dgemv_("t", &kt, &ndl, &done, mat->submat[ip].a2t, &kt, zbut, &one, &done, &zau[nstrtl-1], &one);
    } else if(mat->submat[ip].ltmtx==2){
	  dgemv_("t", &ndt, &ndl, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, &zau[nstrtl-1], &one);
    }
  }
}

void  hmvm_blas_p_calc_2
(double *zau, matrix2 *mat2, double *zu, double *zbut)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill,head;
  int zero = 0;
  int one = 1;
  double dzero = 0.0;
  double done = 1.0;

  nlf=mat2.nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat2.ndl[ip];
    ndt    = mat2.ndt[ip];
    nstrtl = mat2.nstrtl[ip];
    nstrtt = mat2.nstrtt[ip];

    if(mat2.ltmtx[ip]==1){
      kt=mat2.kt[ip];
	  for(il=0;il<kt;il++)zbut[il]=0.0;
	  head = mat2.a1[ip];
	  dgemv_("t", &ndt, &kt, &done, &mat2.rowmat[head], &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
	  head = mat2.a2[ip];
	  dgemv_("n", &ndl, &kt, &done, &mat2.rowmat[head], &ndl, zbut, &one, &done, &zau[nstrtl-1], &one);
    } else if(mat2.ltmtx[ip]==2){
	  head = mat2.a1[ip];
	  dgemv_("t", &ndt, &ndl, &done, &mat2.rowmat[head], &ndt, &zu[nstrtt-1], &one, &done, &zau[nstrtl-1], &one);
    }
  }
}

void  hmvm_blas_p_calc_2t
(double *zau, matrix2 *mat2, double *zu, double *zbut)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill,head;
  int zero = 0;
  int one = 1;
  double dzero = 0.0;
  double done = 1.0;

  nlf=mat2.nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat2.ndl[ip];
    ndt    = mat2.ndt[ip];
    nstrtl = mat2.nstrtl[ip];
    nstrtt = mat2.nstrtt[ip];

    if(mat2.ltmtx[ip]==1){
      kt=mat2.kt[ip];
	  for(il=0;il<kt;il++)zbut[il]=0.0;
	  head = mat2.a1[ip];
	  dgemv_("t", &ndt, &kt, &done, &mat2.rowmat_t[head], &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
	  head = mat2.a2[ip];
	  dgemv_("t", &kt, &ndl, &done, &mat2.rowmat_t[head], &kt, zbut, &one, &done, &zau[nstrtl-1], &one);
    } else if(mat2.ltmtx[ip]==2){
	  head = mat2.a1[ip];
	  dgemv_("t", &ndt, &ndl, &done, &mat2.rowmat_t[head], &ndt, &zu[nstrtt-1], &one, &done, &zau[nstrtl-1], &one);
    }
  }
}

// mkl blas interface
void hmvm_blas_p(matrix *mat, matrix2 *mat2, double *b, int dump_result)
{
  FILE *F;
  int i, l, nd, tmpkt;
  double *v=NULL, *tmp=NULL;
  printf("hmvm_blas_p: begin\n"); fflush(stdout);
  if(mat!=NULL){nd=mat->nd;tmpkt=mat->ktmax;}else{nd=mat2->nd;tmpkt=mat2->ktmax;}
  v=(double*)malloc(sizeof(double)*nd);
  tmp=(double*)malloc(sizeof(double)*tmpkt);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  // blas_p_1
  if(mat!=NULL){
	printf("blas_p_1\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_p_calc_1(v, mat, b, tmp);
	if(dump_result){
	  F = fopen("blas_p_1_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // blas_p_1t
  if(mat!=NULL){
	printf("blas_p_1t\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_p_calc_1t(v, mat, b, tmp);
	if(dump_result){
	  F = fopen("blas_p_1t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // blas_p_2
  if(mat2!=NULL){
	printf("blas_p_2\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_p_calc_2(v, mat2, b, tmp);
	if(dump_result){
	  F = fopen("blas_p_2_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // blas_p_2t
  if(mat2!=NULL){
	printf("blas_p_2t\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_p_calc_2t(v, mat2, b, tmp);
	if(dump_result){
	  F = fopen("blas_p_2t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  free(v); free(tmp);
  printf("hmvm_blas_p: end\n"); fflush(stdout);
}

// mkl blas benchmark interface
void hmvm_blas_p_bench(matrix *mat, matrix2 *mat2, double *b)
{
  const int L=10;
  FILE *F;
  int i, l, nd, tmpkt;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  double *v=NULL, *tmp=NULL;
  printf("hmvm_blas_p_bench: begin\n"); fflush(stdout);
  if(mat!=NULL){nd=mat->nd;tmpkt=mat->ktmax;}else{nd=mat2->nd;tmpkt=mat2->ktmax;}
  v=(double*)malloc(sizeof(double)*nd);
  tmp=(double*)malloc(sizeof(double)*tmpkt);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  // blas_p_1
  {
	printf("blas_p_1\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_p_calc_1(v, mat, b, tmp);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_p_1 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // blas_p_1t
  {
	printf("blas_p_1t\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_p_calc_1t(v, mat, b, tmp);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_p_1t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // blas_p_2
  {
	printf("blas_p_2\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_p_calc_2(v, mat2, b, tmp);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_p_2 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // blas_p_2t
  {
	printf("blas_p_2t\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_p_calc_2t(v, mat2, b, tmp);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_p_2t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  free(v); free(tmp);
  printf("hmvm_blas_p_bench: end\n"); fflush(stdout);
}


// ######## ######## ######## ########
// mkl blas sequential, called from threads
// ######## ######## ######## ########
void  hmvm_blas_s_calc_1
(double *zau, matrix *mat, double *zu)
{
  mkl_set_num_threads(1);
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int zero = 0;
	int one = 1;
	double dzero = 0.0;
	double done = 1.0;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (double*)malloc(sizeof(double)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (double*)malloc(sizeof(double)*mat->ktmax);
	ls = nd;
	le = 1;
#pragma omp for
	for(ip=0; ip<nlf; ip++){
	  ndl   =mat->submat[ip].ndl;
	  ndt   =mat->submat[ip].ndt;
	  nstrtl=mat->submat[ip].nstrtl;
	  nstrtt=mat->submat[ip].nstrtt;
	  if(nstrtl<ls)ls=nstrtl;
	  if(nstrtl+ndl-1>le)le=nstrtl+ndl-1;
	  if(mat->submat[ip].ltmtx==1){
		kt=mat->submat[ip].kt;
		for(il=0;il<kt;il++)zbut[il]=0.0;
		dgemv_("t", &ndt, &kt, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
		dgemv_("n", &ndl, &kt, &done, mat->submat[ip].a2, &ndl, zbut, &one, &done, &zaut[nstrtl-1], &one);
	  } else if(mat->submat[ip].ltmtx==2){
		dgemv_("t", &ndt, &ndl, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, &zaut[nstrtl-1], &one);
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}

void  hmvm_blas_s_calc_1t
(double *zau, matrix *mat, double *zu)
{
  mkl_set_num_threads(1);
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int zero = 0;
	int one = 1;
	double dzero = 0.0;
	double done = 1.0;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (double*)malloc(sizeof(double)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (double*)malloc(sizeof(double)*mat->ktmax);
	ls = nd;
	le = 1;
#pragma omp for
	for(ip=0; ip<nlf; ip++){
	  ndl   =mat->submat[ip].ndl;
	  ndt   =mat->submat[ip].ndt;
	  nstrtl=mat->submat[ip].nstrtl;
	  nstrtt=mat->submat[ip].nstrtt;
	  if(nstrtl<ls)ls=nstrtl;
	  if(nstrtl+ndl-1>le)le=nstrtl+ndl-1;
	  if(mat->submat[ip].ltmtx==1){
		kt=mat->submat[ip].kt;
		for(il=0;il<kt;il++)zbut[il]=0.0;
		dgemv_("t", &ndt, &kt, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
		dgemv_("t", &kt, &ndl, &done, mat->submat[ip].a2t, &kt, zbut, &one, &done, &zaut[nstrtl-1], &one);
	  } else if(mat->submat[ip].ltmtx==2){
		dgemv_("t", &ndt, &ndl, &done, mat->submat[ip].a1, &ndt, &zu[nstrtt-1], &one, &done, &zaut[nstrtl-1], &one);
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}

void  hmvm_blas_s_calc_2
(double *zau, matrix2 *mat2, double *zu)
{
  mkl_set_num_threads(1);
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill, head;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat2.nd;
	int nlf = mat2.nlf;
	int zero = 0;
	int one = 1;
	double dzero = 0.0;
	double done = 1.0;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (double*)malloc(sizeof(double)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (double*)malloc(sizeof(double)*mat2.ktmax);
	ls = nd;
	le = 1;
#pragma omp for
	for(ip=0; ip<nlf; ip++){
	  ndl   =mat2.ndl[ip];
	  ndt   =mat2.ndt[ip];
	  nstrtl=mat2.nstrtl[ip];
	  nstrtt=mat2.nstrtt[ip];
	  if(nstrtl<ls)ls=nstrtl;
	  if(nstrtl+ndl-1>le)le=nstrtl+ndl-1;
	  if(mat2.ltmtx[ip]==1){
		kt=mat2.kt[ip];
		for(il=0;il<kt;il++)zbut[il]=0.0;
		head = mat2.a1[ip];
		dgemv_("t", &ndt, &kt, &done, &mat2.rowmat[head], &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
		head = mat2.a2[ip];
		dgemv_("n", &ndl, &kt, &done, &mat2.rowmat[head], &ndl, zbut, &one, &done, &zaut[nstrtl-1], &one);
	  } else if(mat2.ltmtx[ip]==2){
		head = mat2.a1[ip];
		dgemv_("t", &ndt, &ndl, &done, &mat2.rowmat[head], &ndt, &zu[nstrtt-1], &one, &done, &zaut[nstrtl-1], &one);
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}

void  hmvm_blas_s_calc_2t
(double *zau, matrix2 *mat2, double *zu)
{
  mkl_set_num_threads(1);
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill, head;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat2.nd;
	int nlf = mat2.nlf;
	int zero = 0;
	int one = 1;
	double dzero = 0.0;
	double done = 1.0;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (double*)malloc(sizeof(double)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (double*)malloc(sizeof(double)*mat2.ktmax);
	ls = nd;
	le = 1;
#pragma omp for
	for(ip=0; ip<nlf; ip++){
	  ndl   =mat2.ndl[ip];
	  ndt   =mat2.ndt[ip];
	  nstrtl=mat2.nstrtl[ip];
	  nstrtt=mat2.nstrtt[ip];
	  if(nstrtl<ls)ls=nstrtl;
	  if(nstrtl+ndl-1>le)le=nstrtl+ndl-1;
	  if(mat2.ltmtx[ip]==1){
		kt=mat2.kt[ip];
		for(il=0;il<kt;il++)zbut[il]=0.0;
		head = mat2.a1[ip];
		dgemv_("t", &ndt, &kt, &done, &mat2.rowmat_t[head], &ndt, &zu[nstrtt-1], &one, &done, zbut, &one);
		head = mat2.a2[ip];
		dgemv_("t", &kt, &ndl, &done, &mat2.rowmat_t[head], &kt, zbut, &one, &done, &zaut[nstrtl-1], &one);
	  } else if(mat2.ltmtx[ip]==2){
		head = mat2.a1[ip];
		dgemv_("t", &ndt, &ndl, &done, &mat2.rowmat_t[head], &ndt, &zu[nstrtt-1], &one, &done, &zaut[nstrtl-1], &one);
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}

// mkl blas interface
void hmvm_blas_s(matrix *mat, matrix2 *mat2, double *b, int dump_result)
{
  const int L=15;
  FILE *F;
  int i, l, nd = mat->nd;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  double *v=NULL;
  printf("hmvm_blas_s: begin\n"); fflush(stdout);
  if(mat!=NULL)nd=mat->nd; else nd=mat2->nd
  v=(double*)malloc(sizeof(double)*nd);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  // blas_s_1
  if(mat!=NULL){
	printf("blas_s_1\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_s_calc_1(v, mat, b);
	if(dump_result){
	  F = fopen("blas_s_1_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }
  // blas_s_1t
  if(mat!=NULL){
	printf("blas_s_1t\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_s_calc_1t(v, mat, b);
	if(dump_result){
	  F = fopen("blas_s_1t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }
  // blas_s_2
  if(mat2!=NULL){
	printf("blas_s_2\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_s_calc_2(v, mat2, b);
	if(dump_result){
	  F = fopen("blas_s_2_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }
  // blas_s_2t
  if(mat2!=NULL){
	printf("blas_s_2t\n"); fflush(stdout);
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_blas_s_calc_2t(v, mat2, b);
	if(dump_result){
	  F = fopen("blas_s_2t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  free(v);
  printf("hmvm_blas_s: end\n"); fflush(stdout);
}

// mkl blas benchmark interface
void hmvm_blas_s_bench(matrix *mat, matrix2 *mat2, double *b)
{
  const int L=10;
  int i, l, nd;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  double *v=NULL;
  printf("hmvm_blas_s_bench: begin\n"); fflush(stdout);
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(double*)malloc(sizeof(double)*nd);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  // blas_s_1
  if(mat!=NULL){
	printf("blas_s_1\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_s_calc_1(v, mat, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_s_1 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }
  // blas_s_1t
  if(mat!=NULL){
	printf("blas_s_1t\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_s_calc_1t(v, mat, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_s_1t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }
  // blas_s_2
  if(mat2!=NULL){
	printf("blas_s_2\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_s_calc_2(v, mat2, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_s_2 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }
  // blas_s_2t
  if(mat2!=NULL){
	printf("blas_s_2t\n"); fflush(stdout);
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_blas_s_calc_2t(v, mat2, b);
	  d2 = omp_get_wtime();
	  dtimes[l] = d2-d1;
	}
	dmin = 9999.99;
	dmax = 0.0;
	davg = 0.0;
	for(i=5;i<L;i++){
	  if(dmin>dtimes[i])dmin=dtimes[i];
	  if(dmax<dtimes[i])dmax=dtimes[i];
	  davg += dtimes[i];
	}
	davg /= (L-5);
	printf("TIME %d hmvm_blas_s_2t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  free(v);
  printf("hmvm_blas_s_bench: end\n"); fflush(stdout);
}

