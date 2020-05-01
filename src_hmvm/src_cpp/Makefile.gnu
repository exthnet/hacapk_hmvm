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
#NVCC = /usr/local/cuda/bin/nvcc
NVCC = nvcc
# for Volta
NVCCFLAGS = -g -lineinfo -O3 -Xcompiler "-O3 -fopenmp" \
--generate-code arch=compute_70,code=sm_70
# for Volta (Pascal mode)
#NVCCFLAGS = -g -lineinfo -O0 -Xcompiler "-O3 -fopenmp" \
#--generate-code arch=compute_60,code=sm_70
# for many GPUs
#NVCCFLAGS = -O -Xcompiler "-O3 -fopenmp" \
#--generate-code arch=compute_53,code=sm_53 \
#--generate-code arch=compute_60,code=sm_60 \
#--generate-code arch=compute_70,code=sm_70
# 53=JetsonTX1, 60=Pascal, 70=Volta
TARGET = hmvm_gpu1 hmvm_gpu2
else
CPP  = g++ -std=c++11
CCFLAGS  = -O3 -fopenmp
TARGET = hmvm_cpu1 hmvm_cpu2
endif

######################
# Object files
OBJS1  = loadmatrix.o hmvm1.o hmvm_omp.o hmvm_seq.o
OBJS2  = loadmatrix.o hmvm2.o hmvm_omp.o hmvm_seq.o
CUOBJS = hmvm_cuda.o #hmvm_cuda_kernels.o
HOBJS = hacapk.h

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

hmvm_cpu1: $(OBJS1) $(HOBJS)
	$(CPP) -o hmvm_cpu1 $(OBJS1) -lm $(CCFLAGS)
hmvm_cpu2: $(OBJS2) $(HOBJS)
	$(CPP) -o hmvm_cpu2 $(OBJS2) -lm $(CCFLAGS)
hmvm_gpu1: $(OBJS1) $(HOBJS) $(CUOBJS)
	$(NVCC) -o hmvm_gpu1 $(OBJS1) $(CUOBJS) -lm $(NVCCFLAGS)
hmvm_gpu2: $(OBJS2) $(HOBJS) $(CUOBJS)
	$(NVCC) -o hmvm_gpu2 $(OBJS2) $(CUOBJS) -lm $(NVCCFLAGS)

clean:
	rm -f *.o *.mod *~
distclean:
	rm -f *.o *.mod $(TARGET) *~

rmod:
	rm -f m_*.o *.mod *~

