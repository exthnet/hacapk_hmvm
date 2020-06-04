#ifndef _TEMPLATE_HYBRID3_H
#define _TEMPLATE_HYBRID3_H

#if 0
/*
template <class T>
void hmvm_cuda_hybrid3_proxy
(T *d_zaut, T *d_zu, matrix2<T> *h_mat, matrix2<T> *d_mat,
 int blocks, int threads, int shms, T *v, T *b, char *fname, int bench,
 int div, int mul, int a2t, int a2i, int aatomic, int datomic)
{

for div in 1 2 4 8 16 32
do
for mul in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
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
if(div==${div} && mul==${mul} && a2t==${a2t} && a2i==${a2i} && aatomic==${aatomic} && datomic==${datomic})
hmvm_cuda_hybrid3_proxy<T,${div},${mul},${a2t},${a2i},${aatomic},${datomic}>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, fname, bench);
"
done
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
void hmvm_cuda_hybrid3_proxy
(T *d_zaut, T *d_zu, matrix2<T> *h_mat, matrix2<T> *d_mat,
 int blocks, int threads, int shms, T *v, T *b, char *name, char *fname, int bench,
 int div, int mul, int a2t, int a2i, int aatomic, int datomic)
{
if(div==1 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,1,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,1,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,2,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,2,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,3,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,3,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,4,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,4,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,5,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,5,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,6,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,6,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,7,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,7,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,8,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,8,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,9,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,9,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,10,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,10,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,11,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,11,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,12,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,12,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,13,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,13,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,14,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,14,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,15,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,15,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,1,16,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==1 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,1,16,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,1,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,1,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,2,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,2,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,3,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,3,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,4,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,4,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,5,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,5,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,6,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,6,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,7,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,7,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,8,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,8,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,9,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,9,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,10,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,10,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,11,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,11,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,12,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,12,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,13,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,13,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,14,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,14,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,15,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,15,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,2,16,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==2 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,2,16,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,1,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,1,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,2,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,2,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,3,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,3,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,4,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,4,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,5,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,5,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,6,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,6,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,7,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,7,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,8,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,8,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,9,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,9,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,10,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,10,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,11,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,11,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,12,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,12,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,13,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,13,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,14,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,14,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,15,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,15,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,4,16,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==4 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,4,16,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,1,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,1,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,2,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,2,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,3,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,3,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,4,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,4,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,5,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,5,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,6,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,6,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,7,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,7,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,8,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,8,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,9,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,9,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,10,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,10,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,11,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,11,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,12,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,12,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,13,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,13,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,14,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,14,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,15,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,15,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,8,16,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==8 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,8,16,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,1,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,1,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,2,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,2,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,3,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,3,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,4,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,4,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,5,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,5,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,6,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,6,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,7,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,7,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,8,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,8,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,9,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,9,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,10,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,10,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,11,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,11,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,12,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,12,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,13,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,13,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,14,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,14,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,15,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,15,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,16,16,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==16 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,16,16,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,1,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==1 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,1,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,2,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==2 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,2,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,3,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==3 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,3,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,4,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==4 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,4,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,5,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==5 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,5,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,6,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==6 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,6,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,7,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==7 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,7,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,8,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==8 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,8,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,9,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==9 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,9,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,10,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==10 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,10,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,11,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==11 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,11,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,12,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==12 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,12,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,13,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==13 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,13,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,14,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==14 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,14,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,15,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==15 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,15,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,0,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,0,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,0,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,0,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,0,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,0,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,0,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==0 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,0,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,1,0,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==0 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,1,0,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,1,0,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==0 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,1,0,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,1,1,0,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==1 && aatomic==0 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,1,1,0,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==0)
hmvm_cuda_hybrid3_proxy<T,32,16,1,1,1,0>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);


if(div==32 && mul==16 && a2t==1 && a2i==1 && aatomic==1 && datomic==1)
hmvm_cuda_hybrid3_proxy<T,32,16,1,1,1,1>
(d_zaut, d_zu, h_mat, d_mat, blocks, threads, shms, v, b, name, fname, bench);

}
#endif

#endif
