#!/usr/bin/env python
"""
compare libxc's implementation GGA_X_LB with my own for 
exponentially decaying densities of the type

    rho(r) = N exp(-a * r)

and plot the radial dependence of v_xc[rho]. This reveals numerical
problems for large distances r and narrow densities (large alphas).

van Leeuwen & Baerends proposed a correction to the xc-potential which
recovers the correct asymptotic behaviour (Ref. [1]). The functional
has the form
                 1/3 
   v  [rho] = rho    f(x)
    xc
           __
          |\/ rho|
where x = ---------  is a dimensionless quantity and
          rho^{4/3}
                        x^2
    f(x) = -b ---------------------
               1 + 3 b x arcsinh(x)

This function is designed such that for asymptotically decaying
densities, rho = rho0 * exp(-a*r), the exchange-correlation potential
turns into a Coulomb potential at large distances

   v  [rho] ----->  -1/r       for x ---> oo
    xc

Reference
---------
[1]  R. van Leeuwen and E. J. Baerends,
     "Exchange-correlation potential with correct asymptotic behavior",
     Phys. Rev. A 49, 2421 
"""
from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_integration, multicenter_gradient2, atomlist2arrays
from DFTB.MolecularIntegrals import settings
from DFTB.SlaterKoster import XCFunctionals

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # hydrogen atom
    atomlist = [(1, (0.0, 0.0, 0.0))]

    # choose resolution of multicenter grids for bound orbitals
    settings.radial_grid_factor = 120      # controls size of radial grid  
    settings.lebedev_order = 25          # controls size of angular grid

    
    plt.ylim((-10.0, 2.0))
    plt.title(r"exchange potential for $\rho(r) = N \exp(-\alpha r)$")
    plt.xlabel(r"distance r / bohr")
    plt.ylabel(r"potential / Hartree")

    # define radial grid
    r = np.linspace(0.1, 100.0, 1000)
    # correct asymptotic form of v_xc potential
    plt.plot(r, -1/r, "-", color="black", lw=2, label=r"$-1/r$")

    # As the density becomes narrower (alpha increases), the computation of
    # the LB exchange correlation functional runs into numerical problems,
    # because both (grad rho) and rho are extremely small so that
    #  x = |grad rho|/rho^{4/3} cannot be calculated reliably. Increasing the
    # resolution of the grid does not solve the problem, because with finite
    # precision there will always be a distance at which the the numerical v_xc
    # blows up.
    alphas = [0.1, 0.5, 1.0]
    for alpha in alphas:
        print "density rho(r) = N*exp(-alpha*r) decays with exponential  alpha= %e" % alpha
        # define density function 
        def rho(x,y,z):
            r = np.sqrt(x*x+y*y+z*z)
            n = np.exp(-alpha*r)
            # normalization constant such that
            #   /
            # N | exp(-alpha*r) dV = 1
            #   /
            nrm = alpha**3/(8.0*np.pi)
            return nrm * n

        print "(grad rho)^2 ..."
        # GGA xc-functionals need the gradient squared of the density
        #          __     2
        # sigma = (\/ rho)
        atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
        sigma = multicenter_gradient2(rho,
                                      atomic_coordinates, atomic_numbers,
                                      radial_grid_factor=settings.radial_grid_factor,
                                      lebedev_order=settings.lebedev_order)

        s = sigma(r,0*r,0*r)
        n = rho(r,0*r,0*r)

        # dimensional quantity  x = |grad rho|/rho^{4/3}
        x = np.sqrt(s) / n**(4.0/3.0)

        # van Leeuven and Baerends correction to exchange-correction
        # potential with asymptotic -1/r behaviour. The functional
        # contains one parameter  beta
        beta = 0.05
        vx_LB_correction = - beta * n**(1.0/3.0) * x**2 / (1.0 + 3*beta*x*np.arcsinh(x))

        # LDA exchange
        xc_LDA = XCFunctionals.libXCFunctional("lda_x", "lda_c_xalpha")
        vx_LDA = xc_LDA.func_x.vxc(n.flatten(), s.flatten())
        
        vx_LB = vx_LDA + vx_LB_correction
        
        # libxc's implementation
        xc_LB = XCFunctionals.libXCFunctional("gga_x_lb", "lda_c_xalpha")
        # only exchange part
        vx_LB_libxc = xc_LB.func_x.vxc(n.flatten(), s.flatten())

        # LB approximation
        l, = plt.plot(r, vx_LB,
                      label=r"$v_{x}^{LB}$ $\alpha=%e$" % alpha)
        plt.plot(     r, vx_LB_libxc,
                      ls="-.", label=r"$v_{x}^{LB}$ (libxc) $\alpha=%e$" % alpha, color=l.get_color())

        #
        # Here I tried to define a functional that depends on the
        # center of charge
        #         /                / /
        #   <r> = | rho(r) r d^3r /  | rho(r) d^3r
        #         /              /   /
        # and the width of the charge distribution
        #
        #   dr = sqrt(<(r-<r>)^2>)
        #
        # The ratio chi = dr/<r> is dimensionless.
        #
        # The exchange functional
        #
        #    v_x[rho] = - f(xi) / |r-<r>[rho]|
        # 
        # has the asymptotic behaviour    v_x ---> -1/|r - <r>[rho]|
        # provided that f(chi) ---> 1 for chi ---> 0.
        # The functional also fulfills the scaling relation
        #
        #    v_x(l^3 rho(l*r); r) = l v_x(rho; l*r)
        #
        # since
        #
        #   <r>[l^3 rho(l*r)] = 1/l <r>[rho]
        #
        # and 1/|r-<r>[l^3 rho(l*r)]| = 1/|r - 1/l <r>[rho]|
        #                             = l/[l*r - <r>[rho]|
        #
        # The switching function f(chi) = exp(-omega * chi) satisfies
        #   f(0) = 1   (large radii)   and 
        #   f(inf) = 0  (small radii)
        #
        
        # expectation values <r^2>
        def integrand(x,y,z):
            r2 = x**2+y**2+z**2
            return rho(x,y,z) * r2

        dr2 = multicenter_integration(integrand,
                                atomic_coordinates, atomic_numbers,
                                radial_grid_factor=settings.radial_grid_factor,
                                lebedev_order=settings.lebedev_order)

        dr = np.sqrt(dr2)
        print "sqrt(<r^2>)= %e" % dr
        
        omega = 100.0
        chi = dr/r
        fswitch = np.exp(-omega* chi)
        plt.plot(r, vx_LDA - fswitch / r, ls="--", color=l.get_color(), label=r"$v_{x}^{LDA}$ - $\exp(-\omega (\frac{\Delta r}{r})) \frac{1}{r}$")
        
    plt.legend()
    plt.show()

    
