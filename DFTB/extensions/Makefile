# If the parallelized version of the code exists with SEGFAULT
# you need to increase the stack with
#   export OMP_STACKSIZE=1g
#   ulimit -s unlimited
#
all: thomson.so tddftb.so mulliken.so slako.so grad.so cosmo.so

## Choose if you wish to compile the extension modules with
## 1) GNU's gfortran
SYS=GNU
## or 2) with Intel's compilers and the Math-Kernel library (MKL)
#SYS=INTEL

ifeq ($(SYS),INTEL)
  # use Intel's implementation of OpenMP
  # MKLROOT needs to be set to the installation directory of the MKL library. On the fuxcs and wuxcs this is done
  # by executing the following commands
  #
  #   module load compiler/intel/icc
  #   source $(dirname $(which icc))/compilervars.sh intel64
  #   source $(dirname $(which icc))/../mkl/bin/mklvars.sh intel64

  OPTIMIZATION= --opt="-O3"
  F2PY_OPTIONS= --compiler=intelem --fcompiler=intelem --f90flags="-openmp-report1 -openmp -mkl=parallel -I$(MKLROOT)/include" -DPARALLEL_OMP -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -liomp5 -lmkl_intel_thread -lmkl_mc3 -lmkl_def -lpthread -lm
else
  # GNU's implementation of OpenMP
  # WARNING: If I link with libgomp I get a segmentation
  # fault on the cluster.
  OpenMP=-lgomp
  # Basic linear algebra routines
  BLAS=-lblas
  #
  OPTIMIZATION= --opt="-O3"
  F2PY_OPTIONS= --fcompiler=gfortran --f90flags="-Wall -fopenmp $(OpenMP) -fexternal-blas $(BLAS)" $(OpenMP) $(BLAS) $(OPTIMIZATION)
endif

thomson.so: thomson.f90
	f2py -c thomson.f90 -m thomson $(F2PY_OPTIONS)

tddftb.so: tddftb.f90
	f2py -c tddftb.f90 -m tddftb $(F2PY_OPTIONS)

mulliken.so: mulliken.f90
	f2py -c mulliken.f90 -m mulliken $(F2PY_OPTIONS)

slako.so: slako.f90 splines.f90 splder.f
	f2py -c slako.f90 -m slako $(F2PY_OPTIONS)

grad.so: grad.f90
	f2py -c grad.f90 -m grad $(F2PY_OPTIONS)

cosmo.so: cosmo.f90
	f2py -c cosmo.f90 -m cosmo $(F2PY_OPTIONS)

clean:
	rm thomson.so tddftb.so mulliken.so slako.so grad.so cosmo.so

find_bug.exe:
	gfortran -g -O0 -Wall find_bug.f90 grad.f90 $(BLAS) -fopenmp $(OpenMP) -o find_bug.exe

