#ifndef _TEMPLATE_HYBRID1_H
#define _TEMPLATE_HYBRID1_H

#if 0
/*
template <class T>
void hmvm_cuda_hybrid1_proxy
(T *d_zaut, T *d_zu, matrix2<T> *h_mat, matrix2<T> *d_mat,
 int blocks, int threads, int shms, T *v, T *b, char *fname, int bench,
 int div, int a2t, int a2i, int aatomic, int datomic)
{

for div in 1 2 4 8 16 32
do
for a2t in 0 1
do
for a2i in 0 1
do
for aatomic in 0 1
do
for datomic in 0 1
do
echo "
if(div==${div} && a2t==${a2t} && a2i==${a2i} && aatomic==${aatomic} && datomic==${datomic})
hmvm_cuda_hybrid1_proxy<T,${div},${a2t},${a2i},${aatomic},${datomic}>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, fname, bench);
"
done
done
done
done
done

}
*/
#endif

#if 1
template <class T>
void hmvm_cuda_hybrid1_proxy
(T *d_zaut, T *d_zu, matrix2<T> *h_mat, matrix2<T> *d_mat,
 int blocks, int threads, int shms, T *v, T *b, char *name, char *fname, int bench,
 int div, int a2t, int a2i, int aatomic, int datomic)
{
if(div==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,1,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,1,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,2,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,2,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,4,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,4,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,8,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,8,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,16,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,16,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid1_proxy<T,32,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid1_proxy<T,32,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat,
blocks, threads, shms, v, b, name, fname, bench);
}
#endif

#endif
