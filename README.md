# hacapk_hmvm

## src_hdump

### 概要
HACApK 1.0.0を改造し、H行列をバイナリファイルに出力するようにしたもの。
m_HACApK_hdump.f90以外のコードの多くはHACApK 1.0.0のコードを流用している。

### 使い方
makeするとbem-bb-SCM.outファイルが生成される。
MPI + C 環境であればどこでも使える、と期待される。

./bem-bb-SCM.out inputmatrix
として実行すると、hmatrix_1.binとhmatrix_2.binいうバイナリファイルが出力される。
出力されたバイナリファイルはsrc_hmvmの入力ファイルとして使える。
hmatrix_1.binはFortranの構造体の形状を素直に出力したもの。
hmatrix_2.binはCUDA版実装で扱いやすい形式で出力したもの。
(loadmatrix関数で1から2へ変換できる。)

元のHACApKの都合でMPIが必要だが、このバイナリ出力プログラムが通信を行うわけではないため、
実行時にはmpirun等を噛ませる必要はない。
逆に、MPIマルチプロセス実行することを全く想定していない点には注意。
src_hmvm側も1プロセスで全計算を行う想定で実装されている。

### ソースコード概説 (HACApK 1.0.0からの差分)
- m_ppohBEM_bembb2hacapk.f90
末尾でhacapk_gensolv関数を呼び出していたところを、
hacapk_hdumpを呼び出すように変更してある。
実体はm_HACApK_hdump.f90に書かれている。
- m_HACApK_hdump.f90
HACApKのgensolvとsolveの構造を真似して作成した行列出力処理関数群。

## src_hmvm

### 概要
src_hdumpで出力したバイナリファイルを用いてH行列ベクトル積を行うテストコード群。

### 構成と使い方
#### src_c ディレクトリ
シンプルなC言語コード、CPU(OpenMP含む)向け、倍精度計算のみ。

Makefile.gnuを使ってmake(make -f Makefile.gnu)するとhmvm1とhmvm2が生成される。
src_hdumpのbem-bb-SCM.outで出力したhmatrix_1.binを使って
./hmvm1 ./hmatrix_1.bin
として実行すると計算が行われる。
引数をもう1つ追加し
./hmvm1 ./hmatrix_1.bin 1
とすると、各実装による計算結果がそれぞれファイルに出力される。
ファイル名はresult_omp_1_d.txtなど。

各計算関数には1, 1t, 2, 2tのバリエーションが存在する。
1と2の違いはhmatrix_1.binとhmatrix_2.binの違いに対応しており、
1はarray of structure（matrix構造体とsubmatrix構造体により構成）、
2はstructure of array（matrix2構造体により構成）となっている。
tは近似行列ベクトル積を計算する際に、転置された行列データを用い、ループの交換も行う。
（src_cppでは転置とループ交換それぞれの効果を個別に測定できる。）

- Makefile
  - メッセージ表示とcleanのみ。実際のmakeはMakefile.gnuとMakefile.intelを使って行う。
- Makefile.gnu
  - GNU 環境向けのMakefile。
- Makefile.gnu
  - GNU 環境向けのMakefile。
- hacapk.h
  - H行列の定義。
- loadmatrix.cpp
  - hdumpで出力したバイナリファイルを読み込む機能。
- hmvm1.cpp
  - main関数。
- hmvm_seq.cpp
  - 逐次HMVMコード。
- hmvm_omp.cpp
  - OpenMP版HMVMコード。
- hmvm_mkl.c, hmvm_mkl_cblas.c
  - MKL版HMVMコード。

#### src_cpp ディレクトリ
C++コード、CPU(OpenMP)とGPU(CUDA)、単精度と倍精度それぞれで計算。
論文でやったパラメタチューニング的なものは適用していない。（コードが長くて複雑になるのでどうしようかなぁ。）

Makefile.gnuを使ってmake(make -f Makefile.gnu)するとhmvm_gpuまたはhmvm_cpuが生成される。
src_hdumpのbem-bb-SCM.outで出力したhmatrix_0.binを使って
./hmvm_gpu ./hmatrix_0.bin
として実行すると計算が行われる。
引数をもう1つ追加し
./hmvm_gpu ./hmatrix_0.bin 1
とすると、各実装による計算結果がそれぞれファイルに出力される。

- Makefile.gnu
  - GNU + CUDA 環境向けのMakefile。
    先頭のGPU=1を有効化するとGPU実装が有効化される。
- hacapk.h
  - H行列の定義。
- loadmatrix.cpp
  - hdumpで出力したバイナリファイルを読み込む機能。
- hmvm.cpp
  - main関数。
- hmvm_seq.cpp, hmvm_seq.h
  - 逐次HMVMコード。
- hmvm_omp.cpp, hmvm_omp.h
  - OpenMP版HMVMコード。
- hmvm_cuda.cu, hmvm_cuda.h
  - CUDA版HMVMコード。
- hmvm_cuda_kernels.cu, hmvm_cuda_kernels.h
  - CUDAカーネル。


