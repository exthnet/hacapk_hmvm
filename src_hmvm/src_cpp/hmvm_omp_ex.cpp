#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <omp.h>

#include "hacapk.h"

// ######## ######## ######## ########
template<class T, int a2i>
void hmvm_omp_1(T *zau, const matrix<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (T*)malloc(sizeof(T)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(il=0; il<kt; il++){
			for(it=0; it<ndl; it++){
			  ill=it+nstrtl-1;
			  itl=it+il*ndl;
			  zaut[ill] += mat->submat[ip].a2[itl]*zbut[il];
			}
		  }
		}else{
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			tmp1 = (T)0.0;
			for(il=0; il<kt; il++){
			  itl=it+il*ndl;
			  tmp1 += mat->submat[ip].a2[itl]*zbut[il];
			}
			zaut[ill] += tmp1;
		  }
		}
#endif
	  } else if(mat->submat[ip].ltmtx==2){
#ifndef _SKIP_DENSE
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			tmp1 += mat->submat[ip].a1[itl]*zu[itt];
		  }
		  zaut[ill] += tmp1;
		}
#endif
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
template<class T, int a2i>
void hmvm_omp_1_atomic(T *zau, const matrix<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (T*)malloc(sizeof(T)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(il=0; il<kt; il++){
			for(it=0; it<ndl; it++){
			  ill=it+nstrtl-1;
			  itl=it+il*ndl;
			  //zaut[ill] += mat->submat[ip].a2[itl]*zbut[il];
#pragma omp atomic
			  zau[ill] += mat->submat[ip].a2[itl]*zbut[il];
			}
		  }
		}else{
		  for(it=0; it<ndl; it++){
			for(il=0; il<kt; il++){
			  ill=it+nstrtl-1;
			  itl=it+il*ndl;
			  //zaut[ill] += mat->submat[ip].a2[itl]*zbut[il];
#pragma omp atomic
			  zau[ill] += mat->submat[ip].a2[itl]*zbut[il];
			}
		  }
		}
#endif
	  } else if(mat->submat[ip].ltmtx==2){
#ifndef _SKIP_DENSE
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			//zaut[ill] += mat->submat[ip].a1[itl]*zu[itt];
			tmp1 +=  mat->submat[ip].a1[itl]*zu[itt];
		  }
#pragma omp atomic
		  zau[ill] += tmp1;
		}
#endif
	  }
	}
	/*
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
	*/
    free(zaut); free(zbut);
  }
}

// ######## ######## ######## ########
template<class T, int a2i>
void hmvm_omp_1t(T *zau, const matrix<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (T*)malloc(sizeof(T)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(it=0; it<kt; it++){
			for(il=0; il<ndl; il++){
			  ill=il+nstrtl-1;
			  itl=it+il*kt;
			  zaut[ill] += mat->submat[ip].a2t[itl]*zbut[it];
			}
		  }
		}else{
		  for(il=0; il<ndl; il++){
			ill=il+nstrtl-1;
			tmp1 = (T)0.0;
			for(it=0; it<kt; it++){
			  itl=it+il*kt;
			  tmp1 += mat->submat[ip].a2t[itl]*zbut[it];
			}
			zaut[ill] += tmp1;
		  }
		}
#endif
	  } else if(mat->submat[ip].ltmtx==2){
#ifndef _SKIP_DENSE
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			zaut[ill] += mat->submat[ip].a1[itl]*zu[itt];
		  }
		}
#endif
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
template<class T, int a2i>
void hmvm_omp_1t_atomic(T *zau, const matrix<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (T*)malloc(sizeof(T)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(it=0; it<kt; it++){
			for(il=0; il<ndl; il++){
			  ill=il+nstrtl-1;
			  itl=it+il*kt;
			  //zaut[ill] += mat->submat[ip].a2t[itl]*zbut[it];
#pragma omp atomic
			  zau[ill] += mat->submat[ip].a2t[itl]*zbut[it];
			}
		  }
		}else{
		  for(il=0; il<ndl; il++){
			ill=il+nstrtl-1;
			tmp1 = (T)0.0;
			for(it=0; it<kt; it++){
			  itl=it+il*kt;
			  //zaut[ill] += mat->submat[ip].a2t[itl]*zbut[it];
			  tmp1 += mat->submat[ip].a2t[itl]*zbut[it];
			}
#pragma omp atomic
			zau[ill] += tmp1;
		  }
		}
#endif
	  } else if(mat->submat[ip].ltmtx==2){
#ifndef _SKIP_DENSE
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			//zaut[ill] += mat->submat[ip].a1[itl]*zu[itt];
			tmp1 += mat->submat[ip].a1[itl]*zu[itt];
		  }
#pragma omp atomic
		  zau[ill] += tmp1;
		}
#endif
	  }
	}
	/*
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
	*/
    free(zaut); free(zbut);
  }
}

// ######## ######## ######## ########
template<class T, int a2i>
void hmvm_omp_2(T *zau, const matrix2<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int head;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (T*)malloc(sizeof(T)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(il=0; il<kt; il++){
			for(it=0; it<ndl; it++){
			  ill=it+nstrtl-1;
			  itl=it+il*ndl;
			  zaut[ill] += mat->rowmat[head+itl]*zbut[il];
			}
		  }
		}else{
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			tmp1 = (T)0.0;
			for(il=0; il<kt; il++){
			  itl=it+il*ndl;
			  tmp1 += mat->rowmat[head+itl]*zbut[il];
			}
			zaut[ill] += tmp1;
		  }
		}
#endif
	  } else if(mat->ltmtx[ip]==2){
#ifndef _SKIP_DENSE
		head=mat->a1[ip];
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			tmp1 += mat->rowmat[head+itl]*zu[itt];
		  }
		  zaut[ill] += tmp1;
		}
#endif
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
template<class T, int a2i>
void hmvm_omp_2_atomic(T *zau, const matrix2<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int head;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	//zaut = (T*)malloc(sizeof(T)*nd);
	//for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(il=0; il<kt; il++){
			for(it=0; it<ndl; it++){
			  ill=it+nstrtl-1;
			  itl=it+il*ndl;
			  //zaut[ill] += mat->rowmat[head+itl]*zbut[il];
#pragma omp atomic
			  zau[ill] += mat->rowmat[head+itl]*zbut[il];
			}
		  }
		}else{
		  for(it=0; it<ndl; it++){
			ill=it+nstrtl-1;
			tmp1 = (T)0.0;
			for(il=0; il<kt; il++){
			  itl=it+il*ndl;
			  //zaut[ill] += mat->rowmat[head+itl]*zbut[il];
			  tmp1 += mat->rowmat[head+itl]*zbut[il];
			}
#pragma omp atomic
			zau[ill] += tmp1;
		  }
		}
#endif
	  } else if(mat->ltmtx[ip]==2){
#ifndef _SKIP_DENSE
		head=mat->a1[ip];
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			//zaut[ill] += mat->rowmat[head+itl]*zu[itt];
			tmp1 += mat->rowmat[head+itl]*zu[itt];
		  }
#pragma omp atomic
		  zau[ill] += tmp1;
		  //printf("atomicAdd %d %e\n", ill, tmp);
		}
#endif
	  }
	}
	/*
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
	*/
    //free(zaut);
	free(zbut);
  }
}

// ######## ######## ######## ########
template<class T, int a2i>
void hmvm_omp_2t(T *zau, const matrix2<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int head;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (T*)malloc(sizeof(T)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(it=0; it<kt; it++){
			for(il=0; il<ndl; il++){
			  ill=il+nstrtl-1;
			  itl=it+il*kt;
			  zaut[ill] += mat->rowmat_t[head+itl]*zbut[it];
			}
		  }
		}else{
		  for(il=0; il<ndl; il++){
			ill=il+nstrtl-1;
			tmp1 = (T)0.0;
			for(it=0; it<kt; it++){
			  itl=it+il*kt;
			  tmp1 += mat->rowmat_t[head+itl]*zbut[it];
			}
			zaut[ill] += tmp1;
		  }
		}
#endif
	  }else if(mat->ltmtx[ip]==2){
#ifndef _SKIP_DENSE
		head=mat->a1[ip];
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			tmp1 += mat->rowmat_t[head+itl]*zu[itt];
		  }
		  zaut[ill] += tmp1;
		}
#endif
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
template<class T, int a2i>
void hmvm_omp_2t_atomic(T *zau, const matrix2<T> *mat, const T *zu)
{
#pragma omp parallel
  {
	int ip,il,it;
	int ndl,ndt,nstrtl,nstrtt,kt,itl,itt,ill;

	T tmp1;
	T *zaut, *zbut;
	int ls, le;
	int i;
	int nd = mat->nd;
	int nlf = mat->nlf;
	int head;

#pragma omp for
	for(i=0;i<nd;i++)zau[i]=0.0;

	zaut = (T*)malloc(sizeof(T)*nd);
	for(il=0;il<nd;il++)zaut[il]=0.0;
	zbut = (T*)malloc(sizeof(T)*mat->ktmax);
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
#ifndef _SKIP_APPROX
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
		if(a2i==0){
		  for(it=0; it<kt; it++){
			for(il=0; il<ndl; il++){
			  ill=il+nstrtl-1;
			  itl=it+il*kt;
			  //zaut[ill] += mat->rowmat_t[head+itl]*zbut[it];
#pragma omp atomic
			  zau[ill] += mat->rowmat_t[head+itl]*zbut[it];
			}
		  }
		}else{
		  for(il=0; il<ndl; il++){
			ill=il+nstrtl-1;
			tmp1 = (T)0.0;
			for(it=0; it<kt; it++){
			  itl=it+il*kt;
			  //zaut[ill] += mat->rowmat_t[head+itl]*zbut[it];
			  tmp1 += mat->rowmat_t[head+itl]*zbut[it];
			}
#pragma omp atomic
			zau[ill] += tmp1;
		  }
		}
#endif
	  }else if(mat->ltmtx[ip]==2){
#ifndef _SKIP_DENSE
		head=mat->a1[ip];
		for(il=0; il<ndl; il++){
		  ill=il+nstrtl-1;
		  tmp1 = (T)0.0;
		  for(it=0; it<ndt; it++){
			itt=it+nstrtt-1;
			itl=it+il*ndt;
			//zaut[ill] += mat->rowmat_t[head+itl]*zu[itt];
			tmp1 += mat->rowmat_t[head+itl]*zu[itt];
		  }
#pragma omp atomic
			zau[ill] += tmp1;
		}
#endif
	  }
	}
	/*
    for(il=ls-1;il<=le-1;il++){
#pragma omp atomic
      zau[il] += zaut[il];
    }
	*/
    free(zaut); free(zbut);
  }
}


// ######## ######## ######## ########
// hmvm for mat1
template<class T>
void hmvm_omp_proxy1(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench,
					 void (*fn)(T *, const matrix<T> *, const T *), const char *subname)
{
  int M=5, L, lmax;
  int i, l, nd;
  double d1, d2, *dtimes, dmin, dmax, davg;
  char fname[0xff];
  FILE *F;
  T *v=NULL;
  printf("hmvm_omp_%s: begin\n", typeid(T).name());
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
void hmvm_omp_proxy2(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench,
					 void (*fn)(T *, const matrix2<T> *, const T *), const char *subname)
{
  int M=5, L, lmax;
  int i, l, nd;
  double d1, d2, *dtimes, dmin, dmax, davg;
  char fname[0xff];
  FILE *F;
  T *v=NULL;
  printf("hmvm_omp_%s: begin\n", typeid(T).name());
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
	  fn(v, mat2, b); // hmvm kernel
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
void hmvm_omp_proxy(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench)
{
  // mat1
  // hmvm
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1<T,0>, "omp_1");  }

  // hmvm (loop interchanged)
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1<T,1>, "omp_1i");  }

  // hmvm (trans)
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1t<T,0>, "omp_1t");  }

  // hmvm (trans, loop interchanged)
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1t<T,1>, "omp_1ti");  }

  // mat1 + atomic
  // hmvm
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1_atomic<T,0>, "omp_1_atomic");  }

  // hmvm (loop interchanged)
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1_atomic<T,1>, "omp_1i_atomic");  }

  // hmvm (trans)
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1t_atomic<T,0>, "omp_1t_atomic");  }

  // hmvm (trans, loop interchanged)
  if(mat!=NULL){	hmvm_omp_proxy1(mat, mat2, b, 1, 0, &hmvm_omp_1t_atomic<T,1>, "omp_1ti_atomic");  }

  // mat2
  // hmvm using rowmat array
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2<T,0>, "omp_2");  }

  // hmvm using rowmat array (loop interchanged)
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2<T,1>, "omp_2i");  }

  // hmvm using rowmat array (trans)
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2t<T,0>, "omp_2t");  }

  // hmvm using rowmat array (trans, loop interchanged)
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2t<T,1>, "omp_2ti");  }

  // mat2 + atomic
  // hmvm using rowmat array
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2_atomic<T,0>, "omp_2_atomic");  }

  // hmvm using rowmat array (loop interchanged)
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2_atomic<T,1>, "omp_2i_atomic");  }

  // hmvm using rowmat array (trans)
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2t_atomic<T,0>, "omp_2t_atomic");  }

  // hmvm using rowmat array (trans, loop interchanged)
  if(mat2!=NULL){	hmvm_omp_proxy2(mat, mat2, b, 1, 0, &hmvm_omp_2t_atomic<T,1>, "omp_2ti_atomic");  }
}

// ######## ######## ######## ########
template<class T>
void hmvm_omp(const matrix<T> *mat, const matrix2<T> *mat2, const T *b, int dump_result, int nbench)
{
  printf("hmvm_omp_%s: begin\n", typeid(T).name());
  if(dump_result)hmvm_omp_proxy(mat, mat2, b, 1, 0);
  if(nbench>0)hmvm_omp_proxy(mat, mat2, b, 0, nbench);
  printf("hmvm_omp_%s: end\n", typeid(T).name());
}


// ######## ######## ######## ########
template void hmvm_omp<float> (const matrix<float>  *mat, const matrix2<float>  *mat2, const float  *b, int dump_result, int nbench);
template void hmvm_omp<double>(const matrix<double> *mat, const matrix2<double> *mat2, const double *b, int dump_result, int nbench);
