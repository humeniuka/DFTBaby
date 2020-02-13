"""
test iterative Poisson solvers
"""
import numpy as np
from DFTB.Poisson.poisson_iterative import poisson3d

def test_poisson3d():
    ## 1. no force
    nxpts = 100
    nypts = 100
    nzpts = 100
    xvec = np.linspace(-2.0, 2.0, nxpts)
    yvec = np.linspace(-2.0, 2.0, nypts)
    zvec = np.linspace(-2.0, 2.0, nzpts)
    x,y,z = np.meshgrid(xvec,yvec,zvec,indexing='ij')
    f = 0*x
    u_exact = 2*x+y-3*z

    # boundary conditions...
    u0 = 0*x
    # ...on z-planes
    u0[ :, :, 0] = u_exact[ :, :, 0]
    u0[ :, :,-1] = u_exact[ :, :,-1]
    # ...on y-planes
    u0[ :, 0, :] = u_exact[ :, 0, :]
    u0[ :,-1, :] = u_exact[ :,-1, :]
    # ...on x-planes
    u0[ 0, :, :] = u_exact[ 0, :, :]
    u0[-1, :, :] = u_exact[-1, :, :]    

    # solve
    u = poisson3d(xvec,yvec,zvec, f, u0, maxiter=500000)
    err = np.sum(abs(u-u_exact))/float(x.size)
    print "err = %s" % err
    
    
    
if __name__ == "__main__":
    test_poisson3d()
