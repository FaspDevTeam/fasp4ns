########################################################################
# Compilers and dependences
########################################################################
AR=ar ruc


##########################
# Use GNU compilers
##########################
#CC=gcc
#CPP=g++
#FC=gfortran

##########################
# Use Intel compilers
##########################
CC=icx
CPP=icx
FC=ifort

INCLUDE=-I ../../faspsolver/base/include -I ../include
FASPLIB=../../faspsolver/lib/libfasp.a
TESTLIB=../lib/libfasp4ns.a
FSRCDIR=.

# For Mac OSX users, uncomment the following line:
#BLASLIB=-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.1.sdk -framework Accelerate

##########################
# Use UMFPACK solver
##########################
#SOLVERDIR_UMFPACK = /Users/zhangcs/Packages
#SOLVERLIB_UMFPACK = -L$(SOLVERDIR_UMFPACK)/lib -lsuitesparseconfig -lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd

##########################
# Use MKL solver
##########################
#SOLVERDIR_MKL = /public1/soft/oneAPI/2022.1/mkl/latest
#SOLVERLIB_MKL = /public1/soft/oneAPI/2022.1/mkl/latest/lib/intel64/libmkl_rt.so /public1/soft/oneAPI/2022.1/compiler/latest/linux/compiler/lib/intel64/libiomp5.so /public1/soft/oneAPI/2022.1/compiler/latest/linux/compiler/lib/intel64/libimf.so -lm -Wl,-rpath,::::::::::::::
#SOLVERINCLUDE_MKL = -I$(SOLVERDIR_MKL)/include

SOLVERDIR_MKL = /opt/intel/oneapi/mkl/2021.4.0
SOLVERLIB_MKL = /opt/intel/oneapi/mkl/2021.4.0/lib/intel64/libmkl_rt.so /opt/intel/oneapi/compiler/2021.4.0/linux/compiler/lib/intel64/libiomp5.so /opt/intel/oneapi/compiler/2021.4.0/linux/compiler/lib/intel64/libimf.so -lm -Wl,-rpath,::::::::::::::
SOLVERINCLUDE_MKL = -I$(SOLVERDIR_MKL)/include

##########################
# Use STRUMPACK solver
##########################
SOLVERDIR_STRUM = /data/shared/Solver/STRUMPACK/install/
SOLVERLIB_STRUM = $(SOLVERDIR_STRUM)/lib/libstrumpack.a /data/shared/Solver/metis-5.1.0/install/lib/libmetis.so -lstdc++
SOLVERINCLUDE_STRUM = -I$(SOLVERDIR_STRUM)/include

SOLVERLIB = $(SOLVERLIB_UMFPACK) $(SOLVERLIB_MKL) $(SOLVERLIB_STRUM)
SOLVERINCLUDE = $(SOLVERINCLUDE_UMFPACK) $(SOLVERINCLUDE_MKL) $(SOLVERINCLUDE_STRUM)

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-O3 -fpscomp filesfromcmd -qmkl=sequential 
