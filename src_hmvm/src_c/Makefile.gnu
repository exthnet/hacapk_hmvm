# -*- makefile -*-

CC   = gcc
CCFLAGS  += -Wall -O3 -fopenmp


######################
# Object files
OBJS1  = loadmatrix.o hmvm1.o hmvm_seq.o hmvm_omp.o
OBJS2  = loadmatrix.o hmvm2.o hmvm_seq.o hmvm_omp.o
HEADERS = hacapk.h

######################
# Compile cmmands
.SUFFIXES: .o .c

$(TARGET): $(OBJS)
	$(LINK) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

.c.o: *.c $(HEADERS)
	$(CC)  -c $(CCFLAGS)  $(INCS) $<
.f90.o: *.f90
	$(F90) -c $(F90FLAGS) $(INCS) $<

all: hmvm_cpu1 hmvm_cpu2

hmvm_cpu1: $(OBJS1)
	$(CC) -o hmvm_cpu1 $(OBJS1) -lm $(CCFLAGS)
hmvm_cpu2: $(OBJS2)
	$(CC) -o hmvm_cpu2 $(OBJS2) -lm $(CCFLAGS)

clean:
	rm -f *.o *.mod *~
distclean:
	rm -f *.o *.mod $(TARGET) *~

rmod:
	rm -f m_*.o *.mod *~

