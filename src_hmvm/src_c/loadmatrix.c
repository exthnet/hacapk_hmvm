#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hacapk_c.h"


int loadmatrix(const char *fname, matrix *mat, matrix2 *mat2)
{
  int i;
  FILE *F;
  int irecord;
  int nd, nlf, ktmax, len, offset;
  printf("loadmatrix: begin\n"); fflush(stdout);
  F = fopen(fname, "r");
  if(F==NULL){
	printf("fopen %s failed\n", fname);
	return -1;
  }
  fread(&irecord, sizeof(int), 1, F);
  fread(&nd, sizeof(int), 1, F);
  fread(&nlf, sizeof(int), 1, F);
  fread(&ktmax, sizeof(int), 1, F);
  fread(&irecord, sizeof(int), 1, F);
  printf("nd =  %d\n", nd);
  printf("nlf = %d\n", nlf);
  printf("ktmax = %d\n", ktmax);
  mat->nd = nd;
  mat->nlf = nlf;
  mat->ktmax = ktmax;

  mat->submat = (submatrix*)malloc(sizeof(submatrix)*mat->nlf);
  if(mat->submat==NULL){
	printf("malloc mat->submat failed\n");
	return -1;
  }

  len = 0;

  printf("loadmatrix: load rowdata\n"); fflush(stdout);
  for(i=0;i<mat->nlf;i++){
	int ltmtx, ndl, ndt, nstrtl, nstrtt, kt;
	fread(&irecord, sizeof(int), 1, F);
	fread(&ltmtx, sizeof(int), 1, F);
	fread(&ndl, sizeof(int), 1, F);
	fread(&ndt, sizeof(int), 1, F);
	fread(&nstrtl, sizeof(int), 1, F);
	fread(&nstrtt, sizeof(int), 1, F);
	fread(&kt, sizeof(int), 1, F);
	fread(&irecord, sizeof(int), 1, F);
	//printf("%d: %d %d %d %d %d %d\n", i,ltmtx, ndl, ndt, nstrtl, nstrtt, kt);

	mat->submat[i].ltmtx = ltmtx;
	mat->submat[i].ndl = ndl;
	mat->submat[i].ndt = ndt;
	mat->submat[i].nstrtl = nstrtl;
	mat->submat[i].nstrtt = nstrtt;
	mat->submat[i].kt = kt;

	if(ltmtx==1){
	  mat->submat[i].a1 = (double*)malloc(sizeof(double)*kt*ndt);
	  mat->submat[i].a2 = (double*)malloc(sizeof(double)*kt*ndl);
	  mat->submat[i].a2t = (double*)malloc(sizeof(double)*kt*ndl);
	  fread(&irecord, sizeof(int), 1, F);
	  fread(mat->submat[i].a1, sizeof(double), kt*ndt, F);
	  fread(&irecord, sizeof(int), 1, F);
	  fread(&irecord, sizeof(int), 1, F);
	  fread(mat->submat[i].a2, sizeof(double), kt*ndl, F);
	  fread(&irecord, sizeof(int), 1, F);
	  len += kt*ndt + kt*ndl;
	  int x, y;
	  for(y=0;y<kt;y++){
		for(x=0;x<ndl;x++){
		  mat->submat[i].a2t[y+x*kt] = mat->submat[i].a2[y*ndl+x];
		}
	  }
	}else{
	  mat->submat[i].a1 = (double*)malloc(sizeof(double)*ndl*ndt);
	  mat->submat[i].a2 = NULL;
	  mat->submat[i].a2t = NULL;
	  fread(&irecord, sizeof(int), 1, F);
	  fread(mat->submat[i].a1, sizeof(double), ndl*ndt, F);
	  fread(&irecord, sizeof(int), 1, F);
	  len += ndl*ndt;
	}
  }

  printf("loadmatrix: make matrix2\n"); fflush(stdout);
  mat2->nd = nd;
  mat2->nlf = nlf;
  mat2->ktmax = ktmax;
  mat2->ltmtx = (int*)malloc(sizeof(int)*nlf);
  mat2->kt = (int*)malloc(sizeof(int)*nlf);
  mat2->ndl = (int*)malloc(sizeof(int)*nlf);
  mat2->ndt = (int*)malloc(sizeof(int)*nlf);
  mat2->nstrtl = (int*)malloc(sizeof(int)*nlf);
  mat2->nstrtt = (int*)malloc(sizeof(int)*nlf);
  mat2->a1 = (int*)malloc(sizeof(int)*nlf);
  mat2->a2 = (int*)malloc(sizeof(int)*nlf);
  mat2->rowmat = (double*)malloc(sizeof(double)*len);
  mat2->rowmat_t = (double*)malloc(sizeof(double)*len);
  offset = 0;
  for(i=0;i<nlf;i++){
	mat2->ltmtx[i] = mat->submat[i].ltmtx;
	mat2->kt[i] = mat->submat[i].kt;
	mat2->ndl[i] = mat->submat[i].ndl;
	mat2->ndt[i] = mat->submat[i].ndt;
	mat2->nstrtl[i] = mat->submat[i].nstrtl;
	mat2->nstrtt[i] = mat->submat[i].nstrtt;
	int ltmtx = mat->submat[i].ltmtx;
	if(ltmtx==1){
	  mat2->a1[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndt);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndt);
	  offset += mat->submat[i].kt*mat->submat[i].ndt;
	  mat2->a2[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a2, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndl);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a2t, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndl);
	  offset += mat->submat[i].kt*mat->submat[i].ndl;
	}else{
	  mat2->a1[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].ndl*mat->submat[i].ndt);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].ndl*mat->submat[i].ndt);
	  offset += mat->submat[i].ndl*mat->submat[i].ndt;
	  mat2->a2[i] = 0;
	}
  }
  mat2->len = offset;

  printf("loadmatrix: end\n"); fflush(stdout);
  return 0;
}

double frand(double d)
{
  return d*(rand()%1000)/1000.0;
}

// create ndt*ndl matrix * n*n
int dummymatrix(matrix *mat, matrix2 *mat2, int ndl, int ndt, int n)
{
  int i, x, y, z;
  int offset;
  int len;

  printf("dummymatrix: begin\n");
  mat->nd = n*ndl;
  mat->nlf = n*n;
  mat->ktmax = 0;
  mat->submat = (submatrix*)malloc(sizeof(submatrix)*mat->nlf);
  if(mat->submat==NULL){
	printf("malloc mat->submat failed\n");
	return -1;
  }

  printf("nd =  %d\n", mat->nd);
  printf("nlf = %d\n", mat->nlf);
  printf("ktmax = %d\n", mat->ktmax);

  len = 0;

  srand(0);
  for(y=0;y<n;y++){
	for(x=0;x<n;x++){
	  mat->submat[y*n+x].ltmtx = 2;
	  mat->submat[y*n+x].kt = 0;
	  mat->submat[y*n+x].ndl = ndl;
	  mat->submat[y*n+x].ndt = ndt;
	  mat->submat[y*n+x].nstrtl = ndl*x + 1;
	  mat->submat[y*n+x].nstrtt = ndt*y + 1;
	  mat->submat[y*n+x].a1 = (double*)malloc(sizeof(double)*ndl*ndt);
	  for(z=0;z<ndl*ndt;z++)mat->submat[y*n+x].a1[z] = (double)frand(100.0);
	  mat->submat[y*n+x].a2 = NULL;
	  mat->submat[y*n+x].a2t = NULL;
	  len += ndl*ndt;
	}
  }

  printf("dummymatrix: make matrix2\n");
  mat2->nd = mat->nd;
  mat2->nlf = mat->nlf;
  mat2->ktmax = mat->ktmax;
  mat2->ltmtx = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->kt = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->ndl = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->ndt = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->nstrtl = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->nstrtt = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->a1 = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->a2 = (int*)malloc(sizeof(int)*mat->nlf);
  mat2->rowmat = (double*)malloc(sizeof(double)*len);
  mat2->rowmat_t = (double*)malloc(sizeof(double)*len);
  offset = 0;
  for(i=0;i<mat->nlf;i++){
	mat2->ltmtx[i] = mat->submat[i].ltmtx;
	mat2->kt[i] = mat->submat[i].kt;
	mat2->ndl[i] = mat->submat[i].ndl;
	mat2->ndt[i] = mat->submat[i].ndt;
	mat2->nstrtl[i] = mat->submat[i].nstrtl;
	mat2->nstrtt[i] = mat->submat[i].nstrtt;
	int ltmtx = mat->submat[i].ltmtx;
	if(ltmtx==1){
	  mat2->a1[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndt);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndt);
	  offset += mat->submat[i].kt*mat->submat[i].ndt;
	  mat2->a2[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a2, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndl);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a2t, sizeof(double)*mat->submat[i].kt*mat->submat[i].ndl);
	  offset += mat->submat[i].kt*mat->submat[i].ndl;
	}else{
	  mat2->a1[i] = offset;
	  memcpy(&mat2->rowmat[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].ndl*mat->submat[i].ndt);
	  memcpy(&mat2->rowmat_t[offset], mat->submat[i].a1, sizeof(double)*mat->submat[i].ndl*mat->submat[i].ndt);
	  offset += mat->submat[i].ndl*mat->submat[i].ndt;
	  mat2->a2[i] = 0;
	}
  }
  mat2->len = offset;

  printf("dummymatrix: end\n");
  return 0;
}
