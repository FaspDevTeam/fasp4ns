########################################################################
# Fast Auxiliary Space Preconditioners (FASP)
# 
# This is the Makefile for (Navier-)Stokes test problems! 
# Do NOT put your personal definitions in this general file. Use make.inc instead.
#
########################################################################

########################################################################
# Compilers and dependences
########################################################################
AR=ar ruc
CC=gcc-mp-5
CPP=g++-mp-5
FC=gfortran-mp-5

CSRCDIR=./src
FSRCDIR=./src
INCLUDE=-I ../faspsolver/base/include -I ./include
FASPLIB=../faspsolver/lib/libfasp.a
TESTLIB=./lib/libfasp4ns.a

UMFPACKDIR = /opt/local

BLASLIB = -framework Accelerate

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-g -pg -O3 #-Wall #-fopenmp

COPTS=$(BOPT)
CDEFS=-DWITH_BLAS=1 -DWITH_UMFPACK=1 -DWITH_SuperLU=0
CINCLUDES=$(INCLUDE) $(UMFPACKINCLUDE)
CFLAGS=$(CDEFS) $(COPTS) $(CINCLUDES)

FOPTS=$(BOPT)
FDEFS=$(CDEFS)
FINCLUDES=$(CINCLUDES)
FFLAGS=$(FDEFS) $(FOPTS) $(FINCLUDES)

########################################################################
# Set libraries
########################################################################
LIBS=$(TESTLIB) $(FASPLIB) $(BLASLIB) $(UMFPACKLIB) 

########################################################################
# Load user-defined parameters
########################################################################
include ./make.inc

########################################################################
# Link options
########################################################################
LINKOPTS=$(BOPT)
CLFLAGS=-lstdc++ $(LINKOPTS) $(LIBS)
FLFLAGS=-lm $(LINKOPTS) $(LIBS)

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .cc .cpp .for .f .f77 .f90 .f95
#
FSRC := $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.for))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f77))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f90))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f95))
CSRC := $(foreach dir,$(CSRCDIR),$(wildcard $(CSRCDIR)/*.c))
CSRC += $(foreach dir,$(EXTRDIR),$(wildcard $(EXTRDIR)/*.c))
#
OBJSF := $(patsubst %.for,%.o,$(FSRC))
OBJSF += $(patsubst %.f,%.o,$(FSRC))
OBJSF += $(patsubst %.f77,%.o,$(FSRC))
OBJSF += $(patsubst %.f90,%.o,$(FSRC))
OBJSF += $(patsubst %.f95,%.o,$(FSRC))
OBJSC := $(patsubst %.c,%.o,$(CSRC))
#
.for.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f77.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f90.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f95.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) 
	$(AR) $(TESTLIB) $@
#
.cc.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(TESTLIB) $@
#
.cpp.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(TESTLIB) $@
#
########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(TESTLIB) ns

########################################################################
# Link
########################################################################

all: $(ALLPROG) ns

Default: 
	ns

headers: 
	/bin/cat $(CSRCDIR)/*.c \
	| awk -v name="fasp4ns_functs.h" -f ./bin/mkheaders.awk > ./include/fasp4ns_functs.h

$(TESTLIB): $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)

lib: $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)

########################################################################
# Some test problems
########################################################################

ns: 
	$(CC) $(CFLAGS) -c main/ns.c -o main/ns.o
	$(FC) $(LOPT) main/ns.o $(FLFLAGS) -o ns.ex
	$(FC) $(CFLAGS) -c main/ns.f90 -o main/nsf.o
	$(FC) -o nsf.ex main/nsf.o $(FLFLAGS)

########################################################################
# Clean up
########################################################################

.PHONY : clean distclean help

clean:
	rm -f $(CSRCDIR)/*.o
	rm -f $(FSRCDIR)/*.o
	rm -f main/*.o

distclean:
	make clean
	rm -f lib/*.a
	rm -f *~ *.ex *.out
	rm -f $(CSRCDIR)/*~
	rm -f $(FSRCDIR)/*~

help:
	@echo "======================================================"
	@echo " Fast Auxiliary Space Preconditioners (FASP)"
	@echo "======================================================"
	@echo " "
	@echo " make            : build all exe files "
	@echo " make headers    : build the header file automatically"
	@echo " make clean      : clean all obj files "
	@echo " make allclean   : clean all obj, exe, bak, out files "
	@echo " make help       : show this screen "
	@echo " "
