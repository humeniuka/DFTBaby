#
# compile Fortran extensions into python modules
#

all: ../Wigner3j.so ../coul90.so ../mod_bessel.so ../sphharm.so ../cms.so ../ms.so ../photo.so ../numerov.so

# Wigner 3j symbols
../Wigner3j.so: Wigner3j.f95
	# create signature file
	python -m numpy.f2py Wigner3j.f95 -m Wigner3j -h Wigner3j.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c Wigner3j.pyf Wigner3j.f95
	# move module 
	mv Wigner3j.so ../

# Barnett's implementation of Coulomb and spherical Bessel functions
../coul90.so: COUL90.FOR coul90_wrapper.f90
	# create signature file
	python -m numpy.f2py coul90_wrapper.f90 COUL90.FOR -m coul90 -h coul90.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c coul90.pyf coul90_wrapper.f90 COUL90.FOR
	# move module
	mv coul90.so ../

# modified spherical Bessel functions of 1st and 2nd kind
../mod_bessel.so: mod_bessel.f90
	# create signature file
	python -m numpy.f2py mod_bessel.f90 -m mod_bessel -h mod_bessel.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c mod_bessel.pyf mod_bessel.f90 --f90flags="-fbounds-check"
	# move module
	mv mod_bessel.so ../

# spherical harmonics 
../sphharm.so: sphharm.f90
	# create signature file
	python -m numpy.f2py sphharm.f90 -m sphharm -h sphharm.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c sphharm.pyf sphharm.f90
	# move module
	mv sphharm.so ../

# continuum multiple scattering method
../cms.so: cms.f90 sphharm.f90 coul90_wrapper.f90 COUL90.FOR Wigner3j.f95 cart2sph.f90 splines.f90 splder.f
	# create signature file
	python -m numpy.f2py cms.f90 -m cms -h cms.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c cms.pyf cms.f90 sphharm.f90 coul90_wrapper.f90 COUL90.FOR Wigner3j.f95 cart2sph.f90 splines.f90 --f90flags="-fbounds-check" -llapack
	# move module
	mv cms.so ../

# multiple scattering method (for bound orbitals)
../ms.so: ms.f90 sphharm.f90 mod_bessel.f90 coul90_wrapper.f90 COUL90.FOR Wigner3j.f95 cart2sph.f90 splines.f90 splder.f
	# create signature file
	python -m numpy.f2py ms.f90 -m ms -h ms.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c ms.pyf ms.f90 sphharm.f90 mod_bessel.f90 coul90_wrapper.f90 COUL90.FOR Wigner3j.f95 cart2sph.f90 splines.f90 --f90flags="-fbounds-check" -llapack
	# move module
	mv ms.so ../


# continuum multiple scattering method
../photo.so: photo.f90 Wigner3j.f95 sphharm.f90
	# create signature file
	python -m numpy.f2py photo.f90 -m photo -h photo.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c photo.pyf photo.f90 Wigner3j.f95 sphharm.f90 --f90flags="-fbounds-check"
	# move module
	mv photo.so ../

# radial Schroedinger equation is solved using Numerov's scheme
../numerov.so: numerov.f90
	# create signature file
	python -m numpy.f2py numerov.f90 -m numerov -h numerov.pyf --overwrite-signature
	# compile module
	python -m numpy.f2py -c numerov.pyf numerov.f90 --f90flags="-fbounds-check"
	# move module
	mv numerov.so ../



clean:
	rm -f ../Wigner3j.so ../coul90.so ../mod_bessel.so ../sphharm.so ../cms.so ../ms.so ../photo.so ../numerov.so Wigner3j.pyf coul90.pyf sphharm.pyf cms.pyf photo.pyf numerov.pyf

