#!/usr/bin/env python
# run 
#     python setup.py install --user
# This will install python modules into
#  ~/.local/lib/python2.6/site-packages/
# and scripts into
#  ~/.local/bin/
#

# To gain as much parallelism as possible you should compile 
# numpy and scipy with Intel's Math Kernel Library (MKL)
#
import setuptools
from numpy.distutils.core import Extension, setup

compiler_kwds = {"extra_f90_compile_args": ["-fopenmp", "-fexternal-blas"],
                 "extra_link_args": ["-lgomp", "-lblas", "-llapack"]}

##
# If you wish to compile the extension modules with Intel's compilers as well you should replace 'compiler_kwds' above
# with the following definition:
#
"""
compiler_kwds = {"extra_f90_compile_args": ["-fopenmp -lgomp",
                                            "-fexternal-blas",
                                            "-openmp-report1",
                                            "-openmp",
                                            "-mkl=parallel",
                                            "-I$(MKLROOT)/include",
                                            "-DPARALLEL_OMP"],
                 "extra_link_args": ["-L$(MKLROOT)/lib/intel64",
                                     "-lgomp",
                                     "-lmkl_intel_ilp64",
                                     "-lmkl_core",
                                     "-liomp5",
                                     "-lmkl_intel_thread",
                                     "-lmkl_mc3",
                                     "-lmkl_def",
                                     "-lpthread",
                                     "-lm"]}
"""
##

# extensions written in Fortran

# this extension solves the Thomson problem for creating caps
# of carbon nanotubes
thomson_ext = Extension(name="DFTB.extensions.thomson",
                       sources=["DFTB/extensions/thomson.f90"],
                       **compiler_kwds)

# this extension calculates transition densities faster
tddftb_ext = Extension(name="DFTB.extensions.tddftb",
                       sources=["DFTB/extensions/tddftb.f90"],
                       **compiler_kwds)

# this extension calculates Mulliken monopoles and dipoles faster
mulliken_ext = Extension(name="DFTB.extensions.mulliken",
                       sources=["DFTB/extensions/mulliken.f90"],
                       **compiler_kwds)

# implementation of Slater-Koster rules in Fortran 
slako_ext = Extension(name="DFTB.extensions.slako",
                       sources=["DFTB/extensions/slako.f90"],
                       **compiler_kwds)


# some functions needed for analytic gradients
grad_ext = Extension(name="DFTB.extensions.grad",
                       sources=["DFTB/extensions/grad.f90"],
                       **compiler_kwds)

# force field for periodic QM/MM calculations
ff_ext = Extension(name="DFTB.ForceField.ff",
                   sources=["DFTB/ForceField/src/ff_pythonmodule.c",
                            "DFTB/ForceField/src/ff.c",
                            "DFTB/ForceField/src/linked_list.c",
                            "DFTB/ForceField/src/input.c"],
                   **compiler_kwds)

# iterative Poisson solver for electrostatic potential fitting
poisson_ext = Extension(name="DFTB.Poisson.poisson_iterative",
                        sources=["DFTB/Poisson/src/poisson_iterative.f90"],
                        **compiler_kwds)

# some functions needed for implicit solvent model (COSMO)
cosmo_ext = Extension(name="DFTB.extensions.cosmo",
                       sources=["DFTB/extensions/cosmo.f90"],
                       **compiler_kwds)                   

packages = ['DFTB',
            'DFTB.RepulsivePotential', 'DFTB.RepulsivePotential.reppot_tables',
            'DFTB.SlaterKoster', 'DFTB.SlaterKoster.slako_tables', 'DFTB.SlaterKoster.confined_pseudo_atoms', 'DFTB.SlaterKoster.free_pseudo_atoms',
            'DFTB.Dynamics', 'DFTB.Dynamics.Analyse', 'DFTB.Dynamics.Analyse.Viewer',
            'DFTB.MetaDynamics', 
            'DFTB.Modeling',
            'DFTB.Analyse', 'DFTB.Analyse.blender', 'DFTB.Analyse.vmd', 'DFTB.Analyse.mayavi',
            'DFTB.Formats',
            'DFTB.Formats.Gaussian2py',
            'DFTB.extensions',
            'DFTB.ForceField',
            'DFTB.Poisson',
            'DFTB.MolecularIntegrals',
            'DFTB.MultipleScattering',
            'DFTB.Optimize']

ext_modules = [thomson_ext, tddftb_ext, mulliken_ext, slako_ext, grad_ext, ff_ext, poisson_ext, cosmo_ext]

# current version
version = '0.1.0'
# add a time stamp to the version string and save it in the file VERSION
# so that we know when the distribution was built, even if the version number has not changed
from sys import argv
if "dist" in argv or "sdist" in argv:
    import datetime
    timestamp=datetime.datetime.today().strftime('%Y-%m-%d')
    # write 
    fh=open("VERSION","w")
    fh.write("%s.date-%s\n" % (version, timestamp))
    fh.close()
    
setup(name="DFTBaby",
      version=version,
      description='tight binding density functional theory for excited states',
      author='Alexander Humeniuk',
      author_email='alexander.humeniuk@gmail.com',
      packages=packages,
      ext_modules = ext_modules,
      #package_data={'DFTB': ['hotbit_parameters/*', 'dftbaby.cfg']},
      package_data={'DFTB':                         ['dftbaby.cfg'],
                    'DFTB.Dynamics.Analyse.Viewer': ['*.ui']},
      scripts=['DFTB/DFTB2.py', 'DFTB/LR_TDDFTB.py',
               'DFTB/optimize.py', 'DFTB/optimize_neb.py', 'DFTB/optimize_meci.py', 'DFTB/scan.py',
               'DFTB/Modeling/split_xyz.py', 'DFTB/Modeling/show_molcoords.py', 'DFTB/Modeling/modify_internal.py', 
               'DFTB/Modeling/cut_sphere.py', 'DFTB/Modeling/convert_xyz.py', 
               'DFTB/Modeling/nanotube_builder.py',
               'DFTB/Modeling/filter_fragments.py',
               'DFTB/Dynamics/check_initial_distribution.py',
               'DFTB/Dynamics/random_velocities.py',
               'DFTB/Analyse/absorption_spectrum.py',
               'DFTB/Analyse/ph_correlation.py',
               'DFTB/Dynamics/SurfaceHopping.py',
               'DFTB/Dynamics/Analyse/populations.py',
               'DFTB/Dynamics/Analyse/Viewer/trajViewer2.py',
               'DFTB/MetaDynamics/reconstruct.py',
               'DFTB/Dynamics/dftbaby.sh',
               'DFTB/Optimize/GeometryOptimization.py',
               'DFTB/ForceField/assign_atom_types.py', 'DFTB/ForceField/copy_atom_types.py', 'DFTB/ForceField/ff_optimize.py',
               'DFTB/ForceField/fit_electric_charges.py', 'DFTB/ForceField/fit_magnetic_dipoles.py',
               'DFTB/Poisson/electrostatic_potential.py', 'DFTB/Poisson/vector_potential.py',
               'DFTB/MultipleScattering/electroSka.py']
      )
