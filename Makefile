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
CC=gcc
CPP=g++
FC=gfortran

CSRCDIR=./src
FSRCDIR=./src
INCLUDE=-I ../faspsolver/base/include -I ./include
FASPLIB=../faspsolver/lib/libfasp.a
TESTLIB=./lib/libfasp4ns.a
BLASLIB=-framework Accelerate

########################################################################      
# Directory to UMFPACK
########################################################################        
UMFPACKDIR = /opt/local/
# or to a user lib location like
#    /Users/zhangcs/Packages/

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-O3 #-g -Wall #-fopenmp

COPTS=$(BOPT)
CDEFS=-DWITH_BLAS=1 -DWITH_UMFPACK=0 -DWITH_SuperLU=0
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
ifneq ($(wildcard ./make.inc),)
	include ./make.inc
endif

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
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f90.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F90 object $@'
	@$(AR) $(TESTLIB) $@
#
.f95.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F95 object $@'
	@$(AR) $(TESTLIB) $@
#
.c.o:
	@$(CC) -c $< -o $@ $(CFLAGS)
	@echo 'Building C object $@'
	@$(AR) $(TESTLIB) $@
#
.cpp.o:
	@$(CPP) -c $< -o $@ $(CFLAGS)
	@echo 'Building CPP object $@'
	@$(AR) $(TESTLIB) $@
#
.cc.o:
	@$(CPP) -c $< -o $@ $(CFLAGS)
	@echo 'Building CPP object $@'
	@$(AR) $(TESTLIB) $@
#
########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(TESTLIB) NavierStokes

########################################################################
# Link
########################################################################

all: $(ALLPROG) NavierStokes NavierStokesF StokesBrinkman2 StokesBrinkman3

Default: 
	NavierStokes

headers: 
	cat $(CSRCDIR)/*.c \
	| awk -v name="fasp4ns_functs.h" -f ./util/mkheaders.awk > ./include/fasp4ns_functs.h

$(TESTLIB): $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)

lib: $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)

########################################################################
# Some test problems
########################################################################

NavierStokes:
	@$(CC) $(CFLAGS) -c main/NavierStokes.c -o main/NavierStokes.o
	@$(FC) $(LOPT) -o ns.ex main/NavierStokes.o $(FLFLAGS)
	@echo 'Building executable $@'

NavierStokesF:
	@$(FC) $(CFLAGS) -c main/NavierStokes.f90 -o main/NavierStokesF.o
	@$(FC) -o nsf.ex main/NavierStokesF.o $(FLFLAGS)
	@echo 'Building executable $@'

StokesBrinkman2:
	@$(FC) $(CFLAGS) -c main/SB2DInterfaceBLC.f90 -o main/SB2DInterfaceBLC.o
	@$(FC) $(CFLAGS) -c main/StokesBrinkman2D.f90 -o main/StokesBrinkman2D.o
	@$(FC) -o sb2.ex main/SB2DInterfaceBLC.o main/StokesBrinkman2D.o $(FLFLAGS)
	@echo 'Building executable $@'

StokesBrinkman3:
	@$(FC) $(CFLAGS) -c main/SB3DInterfaceBLC.f90 -o main/SB3DInterfaceBLC.o
	@$(FC) $(CFLAGS) -c main/StokesBrinkman3D.f90 -o main/StokesBrinkman3D.o
	@$(FC) -o sb3.ex main/SB3DInterfaceBLC.o main/StokesBrinkman3D.o $(FLFLAGS)
	@echo 'Building executable $@'

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
	rm -rf *.ex.dSYM

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
