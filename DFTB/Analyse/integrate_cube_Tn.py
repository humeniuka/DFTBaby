#!/usr/bin/env python
"""
integrate a cube-file over the positive quadrante spanned by 2 axes so that a curve along the remaining axis results:

  f(y) = integrate_(0<x<inf, 0<z<inf) rho(x,y,z)

The transition densities of the longer porphyrin tapes are difficult to visualize because one needs a magnifying class to see the density variations. Integrating the transition densities along the short axes might give a better idea of the extension of the excitation along the long axis.
"""

from DFTB.Analyse import Cube
import numpy as np
import numpy.linalg as la
import numpy.fft as fft

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print "Usage: %s <cube file> <axes to integrate over> <output>\n" % sys.argv[0]
        print "  integrates the cube file over the positive quadrant spanned by two axes"
        print "  The planes can be specified by combinations of 'x','y' and 'z', for instance"
        print "    'xy' -> rho(z) = integrate_(0<x<inf, 0<y<inf) f(x,y,z)"
        print "  The 1D function rho(z) is written as a table with rows 'zi rho(zi)' to <output>."
        exit(-1)
    cube_file = sys.argv[1]
    axes_str = sys.argv[2]
    out_file = sys.argv[3]
    
    # Which axes should be integrated out?
    assert len(axes_str) == 2
    ax2index = {"x": 0, "y": 1, "z": 2}
    axes_integrate = [ax2index[s] for s in axes_str]
    # Which axis remains?
    axis_remains = None
    for i in range(0,3):
        if not (i in axes_integrate):
            axis_remains = i
            break
    print "axes %s are integrated out" % str(axes_integrate)
    print "remaining axis: %s" % axis_remains

    # load the cube file
    atomlist, origin, axes, data = Cube.readCube(cube_file)
    axes = [np.array(ax) for ax in axes]

    x0 = np.array(origin)
    nx,ny,nz = data.shape
    ns = [nx,ny,nz]

    # remaining coordinate is called u
    du = la.norm(axes[axis_remains])
    print "du = %s" % du
    nu = ns[axis_remains]
    u = np.linspace(origin[axis_remains], nu*du, nu)
    # shift axis so that the center is located at u=0
    u = u-u[len(u)/2]
    # rho(u), 1D curve after integrating
    rho_u = np.zeros(ns[axis_remains])
    print u.shape
    print rho_u.shape

    for i in xrange(nx):
        for j in xrange(ny):
            for k in xrange(nz):
                index_i = [i,j,k]
                xi = x0 + i*axes[0] + j*axes[1] + k*axes[2]
                for ai in axes_integrate:
                    if xi[ai] < 0.0:
                        break
                else:
                    # xi is in the positive quadrant
                    rho_u[index_i[axis_remains]] += data[i,j,k]

    data_1d = np.vstack((u,rho_u)).transpose()
    fh = open(out_file, "w")
    print>>fh, "# u  in bohr              rho(u)"
    np.savetxt(fh, data_1d)
    fh.close()
    print "Curve saved to %s" % out_file

    # Now the curve is Fourier transformed
    rho_k = fft.rfft(rho_u)
    k = 2.0*np.pi * fft.rfftfreq(len(rho_u),du)
    data_kspace_1d = np.vstack((k.real,rho_k.real,rho_k.imag)).transpose()
    
    out_file_kspace = out_file+"_kspace"
    fh = open(out_file_kspace, "w")
    print>>fh, "# k  in 1/bohr              rho(k)"
    np.savetxt(fh, data_kspace_1d)
    fh.close()
    print "Fourier transform of curve saved to %s" % out_file_kspace
    
    # Analyse rho(k) to find an approximate wavevector for the state
    rho_k2 = abs(rho_k)**2
    # find the first maximum
    for ik in range(1,len(k)):
        if (rho_k2[ik] > rho_k2[ik-1]) and (rho_k2[ik] > rho_k2[ik+1]):
            break
    kmax = k[ik]
    rho_kmax = rho_k[ik]
    print "k max = %s" % kmax
    print "rho(k kmax) = %s" % rho_kmax
    print "wave length lambda = 2*pi/kmax = %s bohr" % (2.0*np.pi/kmax)

    from matplotlib import pyplot as plt
    plt.plot(u,rho_u)
    density_wave = 1.0/float(len(k)) * rho_kmax*(np.exp(1.0j*kmax*u) - np.exp(-1.0j*kmax*u))
    plt.plot(u, density_wave.real)
    plt.plot(u, density_wave.imag)
    plt.show()


    # 
    
