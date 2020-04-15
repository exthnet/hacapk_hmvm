#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "hacapk.h"

// ######## ######## ######## ########
void  hmvm_omp_1
(double *zau, matrix *mat, double *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;

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
		//for(il=0;il<kt;il++)zbut[il]=0.0;
		for(il=0; il<kt; il++){
		  zbut[il]=0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zbut[il] += mat->submat[ip].a1[itl]*zu[itt];
		  }
		}
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			zaut[ill] += mat->submat[ip].a2[itl]*zbut[il];
		  }
		}
	  } else if(mat->submat[ip].ltmtx==2){
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat->submat[ip].a1[itl]*zu[itt];
		  }
		}
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}

// ######## ######## ######## ########
void  hmvm_omp_1t
(double *zau, matrix *mat, double *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;

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
		//for(il=0;il<kt;il++)zbut[il]=0.0;
		for(il=0; il<kt; il++){
		  zbut[il]=0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zbut[il] += mat->submat[ip].a1[itl]*zu[itt];
		  }
		}
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<kt; it++){
			itl=it+il*kt;
			zaut[ill] += mat->submat[ip].a2t[itl]*zbut[it];
		  }
		}
	  } else if(mat->submat[ip].ltmtx==2){
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat->submat[ip].a1[itl]*zu[itt];
		  }
		}
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}

// ######## ######## ######## ########
void  hmvm_omp_2
(double *zau, matrix2 *mat, double *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int head;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (double*)malloc(sizeof(double)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (double*)malloc(sizeof(double)*mat->ktmax);
	ls = nd;
	le = 1;
#pragma omp for
	for(ip=0; ip<nlf; ip++){
	  ndl   =mat->ndl[ip];
	  ndt   =mat->ndt[ip];
	  nstrtl=mat->nstrtl[ip];
	  nstrtt=mat->nstrtt[ip];
	  if(nstrtl<ls)ls=nstrtl;
	  if(nstrtl+ndl-1>le)le=nstrtl+ndl-1;

	  if(mat->ltmtx[ip]==1){
		kt=mat->kt[ip];
		//for(il=0;il<kt;il++)zbut[il]=0.0;
		head=mat->a1[ip];
		for(il=0; il<kt; il++){
		  zbut[il]=0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zbut[il] += mat->rowmat[head+itl]*zu[itt];
		  }
		}
		head=mat->a2[ip];
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			zaut[ill] += mat->rowmat[head+itl]*zbut[il];
		  }
		}
	  } else if(mat->ltmtx[ip]==2){
		head=mat->a1[ip];
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat->rowmat[head+itl]*zu[itt];
		  }
		}
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}

// ######## ######## ######## ########
void  hmvm_omp_2t
(double *zau, matrix2 *mat, double *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int head;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (double*)malloc(sizeof(double)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (double*)malloc(sizeof(double)*mat->ktmax);
	ls = nd;
	le = 1;
#pragma omp for
	for(ip=0; ip<nlf; ip++){
	  ndl   =mat->ndl[ip];
	  ndt   =mat->ndt[ip];
	  nstrtl=mat->nstrtl[ip];
	  nstrtt=mat->nstrtt[ip];
	  if(nstrtl<ls)ls=nstrtl;
	  if(nstrtl+ndl-1>le)le=nstrtl+ndl-1;

	  if(mat->ltmtx[ip]==1){
		kt=mat->kt[ip];
		//for(il=0;il<kt;il++)zbut[il]=0.0;
		head=mat->a1[ip];
		for(il=0; il<kt; il++){
		  zbut[il]=0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zbut[il] += mat->rowmat_t[head+itl]*zu[itt];
		  }
		}
		head=mat->a2[ip];
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<kt; it++){
			itl=it+il*kt;
			zaut[ill] += mat->rowmat_t[head+itl]*zbut[it];
		  }
		}
	  } else if(mat->ltmtx[ip]==2){
		head=mat->a1[ip];
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat->rowmat_t[head+itl]*zu[itt];
		  }
		}
	  }
	}
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
    free(zaut); free(zbut);
  }
}


// ######## ######## ######## ########
void hmvm_omp(matrix *mat, matrix2 *mat2, double *b, int dump_result)
{
  int i, nd;
  FILE *F;
  double *v=NULL;
  printf("hmvm_omp: begin\n");
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(double*)malloc(sizeof(double)*nd);

  // hmvm
  if(mat!=NULL){
	printf("hmvm_omp_1\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_omp_1(v, mat, b);
	if(dump_result){
	  F = fopen("result_omp_1_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm (loop tranposed)
  if(mat!=NULL){
	printf("hmvm_omp_1t\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_omp_1t(v, mat, b);
	if(dump_result){
	  F = fopen("result_omp_1t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm using rowmat array
  if(mat2!=NULL){
	printf("hmvm_omp_2\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_omp_2(v, mat2, b);
	if(dump_result){
	  F = fopen("result_omp_2_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  // hmvm using rowmat array (loop tranposed)
  if(mat2!=NULL){
	printf("hmvm_omp_2t\n");
	for(i=0;i<nd;i++)v[i] = 0.0;
	hmvm_omp_2t(v, mat2, b);
	if(dump_result){
	  F = fopen("result_omp_2t_d.txt", "w");
	  for(i=0;i<nd;i++)fprintf(F, "%.3E\n", v[i]);
	  fclose(F);
	}
  }

  free(v);

  printf("hmvm_omp: end\n");
}

// ######## ######## ######## ########
void hmvm_omp_bench(matrix *mat, matrix2 *mat2, double *b)
{
  const int L=10;
  int i, l, nd;
  double d1, d2, dtimes[L], dmin, dmax, davg;
  double *v=NULL;
  printf("hmvm_omp_bench: begin\n");
  if(mat!=NULL)nd=mat->nd;else nd=mat2->nd;
  v=(double*)malloc(sizeof(double)*nd);

  // hmvm
  if(mat!=NULL){
	printf("hmvm_omp_1\n");
	for(l=0;l<L;l++){
	  for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_omp_1(v, mat, b);
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
	printf("TIME %d hmvm_omp_1 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // hmvm (loop transposed)
  if(mat!=NULL){
	printf("hmvm_omp_1t\n");
	for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	  hmvm_omp_1t(v, mat, b);
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
	printf("TIME %d hmvm_omp_1t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // hmvm using rowmat array
  if(mat!=NULL){
	printf("hmvm_omp_2\n");
	for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	hmvm_omp_2(v, mat2, b);
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
	printf("TIME %d hmvm_omp_2 min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  // hmvm using rowmat array (loop tranposed)
  if(mat!=NULL){
	printf("hmvm_omp_2t\n");
	for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = 0.0;
	  d1 = omp_get_wtime();
	hmvm_omp_2t(v, mat2, b);
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
	printf("TIME %d hmvm_omp_2t min %e max %e avg %e\n", L, dmin, dmax, davg);
  }

  free(v);

  printf("hmvm_omp_bench: end\n");
}
