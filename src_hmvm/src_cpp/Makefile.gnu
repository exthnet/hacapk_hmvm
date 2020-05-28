# -*- makefile -*-

# debug
#OPTS=-D_SKIP_APPROX
#OPTS=-D_SKIP_DENSE

# comment out if you don't want to use GPU
ifndef GPU
	GPU=0
endif

ifeq ($(GPU),1)
CPP  = g++ -std=c++11
CCFLAGS  = -O3 -fopenmp -D_USE_CUDA
NVCC = nvcc
INCS = -I/uvdata/usr/center/a49979a/opt/magma-2.5.3-patched/include
# for Volta
#NVCCFLAGS = -g -lineinfo -O3 -Xcompiler "-O3 -fopenmp" \
#--generate-code arch=compute_70,code=sm_70
NVCCFLAGS = -O3 -Xcompiler "-O3 -fopenmp" \
--generate-code arch=compute_70,code=sm_70
LIBS      = -L/uvdata/usr/center/a49979a/opt/magma-2.5.3-patched/lib -lmagma
# for Volta (Pascal mode)
#NVCCFLAGS = -g -lineinfo -O3 -Xcompiler "-O3 -fopenmp" \
#--generate-code arch=compute_60,code=sm_70
# LIBS      = -L/uvdata/usr/center/a49979a/opt/magma-2.5.3-patched/lib -lmagma
# for many GPUs
#NVCCFLAGS = -g -lineinfo -O3 -Xcompiler "-O3 -fopenmp" \
#--generate-code arch=compute_53,code=sm_53 \
#--generate-code arch=compute_60,code=sm_60 \
#--generate-code arch=compute_70,code=sm_70
# 53=JetsonTX1, 60=Pascal, 70=Volta
# # test
# NVCCFLAGS = -g -lineinfo -O3 -Xcompiler "-O3 -fopenmp -rdynamic" \
# --generate-code arch=compute_70,code=sm_70
# TARGET    = hmvm_gpu1 hmvm_gpu1_ex hmvm_gpu2
# LIBS      = -L/uvdata/usr/center/a49979a/opt/magma-2.5.3-patched/lib -lmagma
else
CPP       = g++ -std=c++11
CCFLAGS   = -O3 -fopenmp
TARGET    = hmvm_cpu1 hmvm_cpu1_ex hmvm_cpu2
endif

######################
# Object files (default)
OBJS1  = loadmatrix.o hmvm1.o hmvm_omp.o hmvm_seq.o
OBJS2  = loadmatrix.o hmvm2.o hmvm_omp.o hmvm_seq.o
# C++ (GPU)
CUOBJS	= hmvm_cuda0.o hmvm_cuda1.o    hmvm_cuda2.o    hmvm_cuda3.o
CUOBJS += hmvm_magma.o hmvm_magma_batched1.o hmvm_magma_batched2.o

# Object files (EX)
OBJS1_EX    = loadmatrix.o hmvm1.o hmvm_omp_ex.o hmvm_seq_ex.o
# C++ (GPU) template extension, requires many long make time
CUOBJS_EX	= hmvm_cuda0_ex.o hmvm_cuda1_ex.o hmvm_cuda2_ex.o hmvm_cuda3_ex.o
CUOBJS_EX  += hmvm_magma_ex.o hmvm_magma_batched1.o hmvm_magma_batched2.o

# common header file(s)
HEADERS  = hacapk.h

######################
# Compile cmmands
.SUFFIXES:
.SUFFIXES: .o .cpp
.SUFFIXES: .o .cu

#$(TARGET): $(OBJS)
#	$(LINK) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

.cpp.o: *.cpp $(HOBJS)
	$(CPP) -c $(CCFLAGS)    $(INCS) $(OPTS) $<
.cu.o: *.cu
	$(NVCC) -c $(NVCCFLAGS) $(INCS) $(OPTS) $<

all: $(TARGET)

hmvm1.o: hacapk.h hmvm_seq.h hmvm_omp.h hmvm_cuda.h hmvm_magma.h
hmvm2.o: hacapk.h hmvm_seq.h hmvm_omp.h hmvm_cuda.h hmvm_magma.h

hmvm_cpu1: $(OBJS1) $(HEADERS)
	$(CPP) -o $@ $(OBJS1) -lm $(CCFLAGS)
hmvm_cpu1_ex: $(OBJS1_EX) $(HEADERS)
	$(CPP) -o $@ $(OBJS1_EX) -lm $(CCFLAGS)
hmvm_cpu2: $(OBJS2) $(HEADERS)
	$(CPP) -o $@ $(OBJS2) -lm $(CCFLAGS)

hmvm_gpu1: $(OBJS1) $(HEADERS) $(CUOBJS)
	$(NVCC) -o $@ $(OBJS1) $(CUOBJS) -lm $(NVCCFLAGS) $(LIBS)
hmvm_gpu1_ex: $(OBJS1_EX) $(HEADERS) $(CUOBJS_EX)
	$(NVCC) -o $@ $(OBJS1_EX) $(CUOBJS_EX) -lm $(NVCCFLAGS) $(LIBS)
hmvm_gpu2: $(OBJS2) $(HEADERS) $(CUOBJS)
	$(NVCC) -o $@ $(OBJS2) $(CUOBJS) -lm $(NVCCFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod *~
distclean:
	rm -f *.o *.mod $(TARGET) *~

rmod:
	rm -f m_*.o *.mod *~

