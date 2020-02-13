# directory where PSPFFT is installed, you probably need to change
# this path
PSPFFT_DIR ?= /net/storage/humeniuka/software/pspfft-0.1/1.0
# directory where FFTW3 library is installed, this probably also needs to be changes
FFTW3_LIB_DIR ?= /usr/lib/x86_64-linux-gnu/
#


LIBS= -L $(PSPFFT_DIR)/lib     -L $(FFTW3_LIB_DIR)
INCS= -I $(PSPFFT_DIR)/include

F2PY_OPTIONS= --fcompiler=gfortran --f90flags="-Wall"

all: ../poisson_iterative.so ../poisson_pspfft.x ../io_pspfft.so

../poisson_iterative.so: poisson_iterative.f90
	f2py -c poisson_iterative.f90 -m poisson_iterative $(F2PY_OPTIONS)
	mv poisson_iterative.so ../

../poisson_pspfft.x: poisson_pspfft.f90
	mpif90 $(LIBS) $(INCS) poisson_pspfft.f90 -o poisson_pspfft.x -l pspfft -l fftw3
	mv poisson_pspfft.x ../

../io_pspfft.so: io_pspfft.f90
	f2py -c io_pspfft.f90 -m io_pspfft $(F2PY_OPTIONS)
	mv io_pspfft.so ../

clean:
	rm -f ../poisson_iterative.so ../poisson_pspfft.x ../io.so
