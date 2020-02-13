"""
solve Poisson equations on a 3-dimensional grid using the PSPFFT code.

The grid with the source data is written to a binary file. Then the external
program ``poisson_pspfft.x`` is called, which writes the solution to a binary file,
which is then read back in.
"""

from DFTB.Poisson import io_pspfft
import numpy as np
import os

def poisson3d(xvec,yvec,zvec, source, dummy, **kwds):
    """
    solve poisson equation
      __2
      \/  potential(x,y,z) = source(x,y,z)

    potential is assumed to decay to 0 for r->oo. 

    Parameters
    ----------
    xvec,yvec,zvec: coordinates of grid points on each axis
    source: source distribution, source[i,j,k]=source(xvec[i],yvec[j],zvec[k])

    Returns
    -------
    potential: solution, potential[i,j,k] = potential(xvec[i],yvec[j],zvec[k])
    """
    # The dummy argument is there, so that the interface is the same as for `poisson_iterative.poisson3d`
    
    # write grid to binary file
    io_pspfft.write_grid(xvec,yvec,zvec, source, 'source.dat')
    # solve Poisson equation
    # poisson_pspfft.x reads the source distribution from 'source.dat'
    # and writes the potential to 'potential.dat'
    os.system('poisson_pspfft.x')
    potential = io_pspfft.read_grid('potential.dat', xvec,yvec,zvec)
    # remove temporary files
    os.remove('source.dat')
    os.remove('source.dat.header')
    os.remove('potential.dat')
    os.remove('potential.dat.header')

    return potential

