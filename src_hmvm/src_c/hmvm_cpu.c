#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "hacapk_c.h"

// ######## ######## ######## ########
// sequential
// zau : result(nd)
// mat : matrix
// zu  : vector(nd)
// zbu : tmp(nd)
// ######## ######## ######## ########
void  hmvm_cpu_calc_1
(double *zau, matrix mat, double *zu, double *zbu)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

  nlf=mat.nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat.submat[ip].ndl;
    ndt    = mat.submat[ip].ndt;
    nstrtl = mat.submat[ip].nstrtl;
    nstrtt = mat.submat[ip].nstrtt;

    if(mat.submat[ip].ltmtx==1){
      kt=mat.submat[ip].kt;
      for(il=0; il<kt; il++){
		zbu[il]=0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zbu[il] += mat.submat[ip].a1[itl]*zu[itt];
		}
      }
      for(il=0; il<kt; il++){
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  itl=it+il*ndl;
		  zau[ill] += mat.submat[ip].a2[itl]*zbu[il];
		}
      }
    } else if(mat.submat[ip].ltmtx==2){
      for(il=0; il<ndl; il++){
		ill=il+nstrtl-1;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zau[ill] += mat.submat[ip].a1[itl]*zu[itt];
		}
      }
    }
  }
}

// ######## ######## ######## ########
void  hmvm_cpu_calc_1t
(double *zau, matrix mat, double *zu, double *zbu)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

  nlf=mat.nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat.submat[ip].ndl;
    ndt    = mat.submat[ip].ndt;
    nstrtl = mat.submat[ip].nstrtl;
    nstrtt = mat.submat[ip].nstrtt;

    if(mat.submat[ip].ltmtx==1){
      kt=mat.submat[ip].kt;
      for(il=0; il<kt; il++){
		zbu[il]=0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zbu[il] += mat.submat[ip].a1[itl]*zu[itt];
		}
      }
      for(il=0; il<ndl; il++){
		ill=il+nstrtl-1;
		for(it=0; it<kt; it++){
		  itl=it+il*kt;
		  zau[ill] += mat.submat[ip].a2t[itl]*zbu[it];
		}
      }
    } else if(mat.submat[ip].ltmtx==2){
      for(il=0; il<ndl; il++){
		ill=il+nstrtl-1;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zau[ill] += mat.submat[ip].a1[itl]*zu[itt];
		}
      }
    }
  }
}

// ######## ######## ######## ########
void  hmvm_cpu_calc_2
(double *zau, matrix2 mat, double *zu, double *zbu)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;
  int head;

  nlf=mat.nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat.ndl[ip];
    ndt    = mat.ndt[ip];
    nstrtl = mat.nstrtl[ip];
    nstrtt = mat.nstrtt[ip];

    if(mat.ltmtx[ip]==1){
      kt=mat.kt[ip];
	  head = mat.a1[ip];
      for(il=0; il<kt; il++){
		zbu[il]=0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zbu[il] += mat.rowmat[head+itl]*zu[itt];
		}
      }
	  head = mat.a2[ip];
      for(il=0; il<kt; il++){
		for(it=0; it<ndl; it++){
		  ill=it+nstrtl-1;
		  itl=it+il*ndl;
		  zau[ill] += mat.rowmat[head+itl]*zbu[il];
		}
      }
    } else if(mat.ltmtx[ip]==2){
	  head = mat.a1[ip];
      for(il=0; il<ndl; il++){
		ill=il+nstrtl-1;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zau[ill] += mat.rowmat[head+itl]*zu[itt];
		}
      }
    }
  }
}

// ######## ######## ######## ########
void  hmvm_cpu_calc_2t
(double *zau, matrix2 mat, double *zu, double *zbu)
{
  int ip,il,it;
  int nlf,ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;
  int head;

  nlf=mat.nlf;

  for(ip=0; ip<nlf; ip++){
    ndl    = mat.ndl[ip];
    ndt    = mat.ndt[ip];
    nstrtl = mat.nstrtl[ip];
    nstrtt = mat.nstrtt[ip];

    if(mat.ltmtx[ip]==1){
      kt=mat.kt[ip];
	  head=mat.a1[ip];
      for(il=0; il<kt; il++){
		zbu[il]=0.0;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zbu[il] += mat.rowmat_t[head+itl]*zu[itt];
		}
      }
	  head=mat.a2[ip];
      for(il=0; il<ndl; il++){
		ill=il+nstrtl-1;
		for(it=0; it<kt; it++){
		  itl=it+il*kt;
		  zau[ill] += mat.rowmat_t[head+itl]*zbu[it];
		}
      }
    } else if(mat.ltmtx[ip]==2){
	  head=mat.a1[ip];
      for(il=0; il<ndl; il++){
		ill=il+nstrtl-1;
		for(it=0; it<ndt; it++){
		  itt=it+nstrtt-1;
		  itl=it+il*ndt;
		  zau[ill] += mat.rowmat_t[head+itl]*zu[itt];
		}
      }
    }
  }
}

// ######## ######## ######## ######## ######## ######## ######## ########

void hmvm_cpu(matrix mat, matrix2 mat2, double *b)
{
  const int L=15;
  FILE *F;
  int i, l, nd = mat.nd;
  double d1, d2, dtimes[L], dmin, dmax, davg1, davg2;
  double *v=NULL, *tmp=NULL;
  printf("hmvm_cpu: begin\n");
  v=(double*)malloc(sizeof(double)*mat.nd);
  tmp=(double*)malloc(sizeof(double)*mat.nd);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  // cpu 1
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_cpu_calc_1(v, mat, b, tmp);
  F = fopen("hmvm_cpu_1.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = 0.0;
	d1 = omp_get_wtime();
	hmvm_cpu_calc_1(v, mat, b, tmp);
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }

  dmin = 9999.99;
  dmax = 0.0;
  davg1 = 0.0;
  davg2 = 0.0;
  for(i=0;i<L;i++){
	if(dmin>dtimes[i])dmin=dtimes[i];
	if(dmax<dtimes[i])dmax=dtimes[i];
	davg1 += dtimes[i];
  }
  for(i=2;i<L-2;i++){
	davg2 += dtimes[i];
  }
  davg1 /= L;
  davg2 /= (L-4);
  printf("TIME %d hmvm_cpu_1 min %e max %e avg1 %e avg2 %e\n", L, dmin, dmax, davg1, davg2);

  // cpu 1t
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_cpu_calc_1t(v, mat, b, tmp);
  F = fopen("hmvm_cpu_1t.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = 0.0;
	d1 = omp_get_wtime();
	hmvm_cpu_calc_1t(v, mat, b, tmp);
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }

  dmin = 9999.99;
  dmax = 0.0;
  davg1 = 0.0;
  davg2 = 0.0;
  for(i=0;i<L;i++){
	if(dmin>dtimes[i])dmin=dtimes[i];
	if(dmax<dtimes[i])dmax=dtimes[i];
	davg1 += dtimes[i];
  }
  for(i=2;i<L-2;i++){
	davg2 += dtimes[i];
  }
  davg1 /= L;
  davg2 /= (L-4);
  printf("TIME %d hmvm_cpu_1t min %e max %e avg1 %e avg2 %e\n", L, dmin, dmax, davg1, davg2);

  // cpu 2
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_cpu_calc_2(v, mat2, b, tmp);
  F = fopen("hmvm_cpu_2.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = 0.0;
	d1 = omp_get_wtime();
	hmvm_cpu_calc_2(v, mat2, b, tmp);
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }

  dmin = 9999.99;
  dmax = 0.0;
  davg1 = 0.0;
  davg2 = 0.0;
  for(i=0;i<L;i++){
	if(dmin>dtimes[i])dmin=dtimes[i];
	if(dmax<dtimes[i])dmax=dtimes[i];
	davg1 += dtimes[i];
  }
  for(i=2;i<L-2;i++){
	davg2 += dtimes[i];
  }
  davg1 /= L;
  davg2 /= (L-4);
  printf("TIME %d hmvm_cpu_2 min %e max %e avg1 %e avg2 %e\n", L, dmin, dmax, davg1, davg2);

  // cpu2t
  for(i=0;i<nd;i++)v[i] = 0.0;
  hmvm_cpu_calc_2t(v, mat2, b, tmp);
  F = fopen("hmvm_cpu_2t.txt", "w");
  for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]);
  fclose(F);

  for(l=0;l<L;l++){
	for(i=0;i<nd;i++)v[i] = 0.0;
	d1 = omp_get_wtime();
	hmvm_cpu_calc_2t(v, mat2, b, tmp);
	d2 = omp_get_wtime();
	dtimes[l] = d2-d1;
  }

  dmin = 9999.99;
  dmax = 0.0;
  davg1 = 0.0;
  davg2 = 0.0;
  for(i=0;i<L;i++){
	if(dmin>dtimes[i])dmin=dtimes[i];
	if(dmax<dtimes[i])dmax=dtimes[i];
	davg1 += dtimes[i];
  }
  for(i=2;i<L-2;i++){
	davg2 += dtimes[i];
  }
  davg1 /= L;
  davg2 /= (L-4);
  printf("TIME %d hmvm_cpu_2t min %e max %e avg1 %e avg2 %e\n", L, dmin, dmax, davg1, davg2);

  free(v); free(tmp);
  printf("hmvm_cpu: end\n");
}


// ######## ######## ######## ########
// OpenMP
// ######## ######## ######## ########


// ######## ######## ######## ########
void  hmvm_omp_calc_3
(double *zau, matrix mat, double *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	double *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat.nd;
	int nlf = mat.nlf;
	int b1, b2, b3;
	int e1, e2, e3;

	b1 = mat.b1;	b2 = mat.b2;	b3 = mat.b3;
	e1 = mat.e1;	e2 = mat.e2;	e3 = mat.e3;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (double*)malloc(sizeof(double)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (double*)malloc(sizeof(double)*mat.ktmax);

	for(ip=b1; ip<e1; ip++){
	  ndl   =mat.submat[ip].ndl;
	  ndt   =mat.submat[ip].ndt;
	  nstrtl=mat.submat[ip].nstrtl;
	  nstrtt=mat.submat[ip].nstrtt;

	  if(mat.submat[ip].ltmtx==1){
		kt=mat.submat[ip].kt;
		//for(il=0;il<kt;il++)zbut[il]=0.0;
		for(il=0; il<kt; il++){
		  zbut[il]=0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zbut[il] += mat.submat[ip].a1[itl]*zu[itt];
		  }
		}
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			zaut[ill] += mat.submat[ip].a2[itl]*zbut[il];
		  }
		}
	  } else if(mat.submat[ip].ltmtx==2){
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat.submat[ip].a1[itl]*zu[itt];
		  }
#pragma omp atomic
		  zau[il] += zaut[il];
		}
	  }
	}

	for(ip=b2; ip<e2; ip++){
	  ndl   =mat.submat[ip].ndl;
	  ndt   =mat.submat[ip].ndt;
	  nstrtl=mat.submat[ip].nstrtl;
	  nstrtt=mat.submat[ip].nstrtt;

	  if(mat.submat[ip].ltmtx==1){
		kt=mat.submat[ip].kt;
		//for(il=0;il<kt;il++)zbut[il]=0.0;
		for(il=0; il<kt; il++){
		  zbut[il]=0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zbut[il] += mat.submat[ip].a1[itl]*zu[itt];
		  }
		}
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			zaut[ill] += mat.submat[ip].a2[itl]*zbut[il];
		  }
		}
	  } else if(mat.submat[ip].ltmtx==2){
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat.submat[ip].a1[itl]*zu[itt];
		  }
		  //#pragma omp atomic
		  zau[il] += zaut[il];
		}
	  }
	}

	for(ip=b3; ip<e3; ip++){
	  ndl   =mat.submat[ip].ndl;
	  ndt   =mat.submat[ip].ndt;
	  nstrtl=mat.submat[ip].nstrtl;
	  nstrtt=mat.submat[ip].nstrtt;

	  if(mat.submat[ip].ltmtx==1){
		kt=mat.submat[ip].kt;
		//for(il=0;il<kt;il++)zbut[il]=0.0;
		for(il=0; il<kt; il++){
		  zbut[il]=0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zbut[il] += mat.submat[ip].a1[itl]*zu[itt];
		  }
		}
		for(il=0; il<kt; il++){
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			itl=it+il*ndl;
			zaut[ill] += mat.submat[ip].a2[itl]*zbut[il];
		  }
		}
	  } else if(mat.submat[ip].ltmtx==2){
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat.submat[ip].a1[itl]*zu[itt];
		  }
#pragma omp atomic
		  zau[il] += zaut[il];
		}
	  }
	}

    free(zaut); free(zbut);
  }
}
// ######## ######## ######## ########

void hmvm_omp(matrix mat, matrix2 mat2, double *b, int kernel)
{
  const int L=5, M=5;
  FILE *F;
  int i, l, nd = mat.nd;
  double d1, d2, dtimes[L+M], dmin, dmax, davg1, davg2;
  double *v=NULL, *tmp=NULL;
  printf("hmvm_omp: begin\n"); fflush(stdout);
  v=(double*)malloc(sizeof(double)*mat.nd);
  tmp=(double*)malloc(sizeof(double)*mat.nd);
  for(i=0;i<nd;i++){
	b[i] = sin((double)(i+1));
  }

  char name[8];
  char fname[32];

#define BENCH(FUNCNAME,MATRIX)						\
	for(i=0;i<nd;i++)v[i] = 0.0; \
    FUNCNAME(v, MATRIX, b);				\
	F = fopen(fname, "w"); \
	for(i=0;i<nd;i++)fprintf(F, "%E\n", v[i]); \
	fclose(F); \
	for(l=0;l<M+L;l++){ \
	  for(i=0;i<nd;i++)v[i] = 0.0; \
	  d1 = omp_get_wtime(); \
	  FUNCNAME(v, MATRIX, b);				\
	  d2 = omp_get_wtime(); \
	  dtimes[l] = d2-d1; \
	} \
 \
	dmin = 9999.99; \
	dmax = 0.0; \
	davg1 = 0.0; \
	davg2 = 0.0; \
	for(i=0;i<M+L;i++)davg1 += dtimes[i];		\
	for(i=M;i<M+L;i++){ \
	  if(dmin>dtimes[i])dmin=dtimes[i]; \
	  if(dmax<dtimes[i])dmax=dtimes[i]; \
	  davg2 += dtimes[i]; \
	} \
	davg1 /= (M+L);								\
	davg2 /= L; \
	printf("TIME omp %06x %d %s min %e max %e avg1 %e avg2 %e |", kernel, M+L, name, dmin, dmax, davg1, davg2); \
	for(i=0;i<M+L;i++)printf(" %e", dtimes[i]);							\
	printf("\n");

  switch(kernel){
  case 0:
  // omp 1
  {
	snprintf(fname,32,"hmvm_omp_1.txt");
	//BENCH(hmvm_omp_calc_1, mat);
	printf("TIME kernel omp1 %d hmvm_omp_1 min %e max %e avg1 %e avg2 %e\n", M+L, dmin, dmax, davg1, davg2);
  }

  // omp 1t
  {
	snprintf(fname,32,"hmvm_omp_1t.txt");
	//BENCH(hmvm_omp_calc_1t, mat);
	printf("TIME kernel omp1t %d hmvm_omp_1t min %e max %e avg1 %e avg2 %e\n", M+L, dmin, dmax, davg1, davg2);
  }

  // omp 2
  {
	snprintf(fname,32,"hmvm_omp_2.txt");
	//BENCH(hmvm_omp_calc_2, mat2);
	printf("TIME kernel omp2 %d hmvm_omp_2 min %e max %e avg1 %e avg2 %e\n", M+L, dmin, dmax, davg1, davg2);
  }

  // omp 2t
  {
	snprintf(fname,32,"hmvm_omp_2t.txt");
	//BENCH(hmvm_omp_calc_2t, mat2);
	printf("TIME kernel omp2t %d hmvm_omp_2t min %e max %e avg1 %e avg2 %e\n", M+L, dmin, dmax, davg1, davg2);
  }
	break;

  case 0x0100:
	snprintf(name,8,"omp1");
	snprintf(fname,32,"hmvm_omp_1.txt");
	BENCH(hmvm_omp_calc_1, mat);
	//printf("TIME kernel %06x %d hmvm_omp_1 min %e max %e avg1 %e avg2 %e\n", kernel, M+L, dmin, dmax, davg1, davg2);
	break;
  case 0x0101:
	snprintf(name,8,"omp1a");
	snprintf(fname,32,"hmvm_omp_1_approx.txt");
	//BENCH(hmvm_omp_calc_1_approx, mat);
	//printf("TIME kernel %06x %d hmvm_omp_1_approx min %e max %e avg1 %e avg2 %e\n", kernel, M+L, dmin, dmax, davg1, davg2);
	break;
  case 0x0102:
	snprintf(name,8,"omp1d");
	snprintf(fname,32,"hmvm_omp_1_dense.txt");
	//BENCH(hmvm_omp_calc_1_dense, mat);
	//printf("TIME kernel %06x %d hmvm_omp_1_dense min %e max %e avg1 %e avg2 %e\n", kernel, M+L, dmin, dmax, davg1, davg2);
	break;

  case 0x0200:
	snprintf(name,8,"omp2");
	snprintf(fname,32,"hmvm_omp_2.txt");
	BENCH(hmvm_omp_calc_2, mat2);
	//printf("TIME kernel %06x %d hmvm_omp_2 min %e max %e avg1 %e avg2 %e\n", kernel, M+L, dmin, dmax, davg1, davg2);
	break;
  case 0x0201:
	snprintf(name,8,"omp2a");
	snprintf(fname,32,"hmvm_omp_2_approx.txt");
	BENCH(hmvm_omp_calc_2_approx, mat2);
	//printf("TIME kernel %06x %d hmvm_omp_2_approx min %e max %e avg1 %e avg2 %e\n", kernel, M+L, dmin, dmax, davg1, davg2);
	break;
  case 0x0202:
	snprintf(name,8,"omp2d");
	snprintf(fname,32,"hmvm_omp_2_dense.txt");
	BENCH(hmvm_omp_calc_2_dense, mat2);
	//printf("TIME kernel %06x %d hmvm_omp_2_dense min %e max %e avg1 %e avg2 %e\n", kernel, M+L, dmin, dmax, davg1, davg2);
	break;
  }


  free(v); free(tmp);
  printf("hmvm_omp: end\n"); fflush(stdout);
}


