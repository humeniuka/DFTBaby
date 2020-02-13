from distutils.core import setup, Extension

compiler_kwds = {"extra_link_args": ["-shared", "-L/local_scratch/humeniuka/software/qc/libxc/lib", "-lxc"],
                 "extra_compile_args": ["-I/usr/include/python2.7/numpy", "-I/local_scratch/humeniuka/software/qc/libxc/include"]}

libxc_extension = Extension(name="_pylibxc",
                            sources=["pylibxc.cpp"],
                            **compiler_kwds)

setup(name = "_pylibxc",
      version = "0.0.1",
      description = "python wrapper for library of exchange-correlation functionals",
      ext_modules = [libxc_extension])
