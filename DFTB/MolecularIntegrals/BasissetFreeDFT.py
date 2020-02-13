#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
basis-set-free local density functional calculations as described
in Ref. [3] and [4]

The electron-electron interaction is treated differently depending
on the number of electrons `n`:

  n = 1: In one-electron systems there is no electron-electron
         interaction, the single-particle Hamiltonian reads  
           H = T + Vnuc
  n = 2: For closed-shell two-electron systems Hartree-Fock theory is 
         used, since the exchange integral removes the self-interaction:
           H = T + Vnuc + 1/2 Vcoul[rho]
  n > 2: For many-electron systems we use DFT with a local 
         exchange-correlation functional:
           H = T + Vnuc + Vcoul[rho] + Vxc[rho]

References
----------
[3] A.Becke, R.Dickson, "Numerical solution of Schroedinger's equation in polyatomic molecules",
    J.Chem.Phys. 92, 3610 (1990)
[4] R.Dickson, A.Becke, "Basis-set-free local density-functional calculations of geometries of polyatomic molecules",
    J.Chem.Phys. 99, 3898 (1993)
"""
import numpy as np
import scipy.linalg as sla
from scipy import optimize

from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_poisson, multicenter_laplacian, multicenter_integration, multicenter_inhomogeneous_schroedinger, multicenter_continuum_schroedinger, multicenter_spherical_remainder, multicenter_operation, multicenter_gradient2, multicenter_residual, spherical_average_func, radial_component_func, atomlist2arrays, print_grid_summary
from DFTB.MolecularIntegrals.Ints1e import integral, overlap, kinetic, nuclear, nuclear_repulsion
from DFTB.MolecularIntegrals import settings
from DFTB.SlaterKoster import XCFunctionals
from DFTB.utils import annotated_matrix
from DFTB import AtomicData
# for continuum orbitals
from DFTB.MolecularIntegrals.SphericalCoords import cartesian2spherical
from DFTB.MolecularIntegrals.LebedevQuadrature import spherical_harmonics_it

def density_func(orbitals):
    """
    create function for evaluating the restricted (rho_a = rho_b) electron density
    for a closed shell singlet ground state

    NOTE:
        The charge of the electron is negative, but the electron density
        is positive everywhere and integrates to the number of electrons.

    Parameters
    ==========
    orbitals      :  doubly occupied orbitals, list callables,
                     orbitals[i](x,y,z) evaluates the i-th occupied orbital
    
    Returns
    =======
    density       :  callable, density(x,y,z) evaluates the total electron density
    """
    
    def density(x,y,z):
        rho = 0.0*x
        for orb in orbitals:
            # orbitals are doubly occupied
            rho += 2*abs(orb(x,y,z))**2
        return rho
    
    return density

def total_charge(atomlist, rho, verbose=1):
    """
    compute the total charge of the molecule as the difference between the
    sum of the nuclear charges (positive) and the number of electrons obtained
    by integrating the electron density numerically.

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    rho          :  callable, rho(x,y,z) evaluates the electron density on a grid,
                    None if only the nuclear charge should be calculated

    Optional
    ========
    verbose      :  if > 0, print charges

    Returns
    =======
    qnuc, qelec, qtot:  nuclear (positive), electronic (negative) and total charge
    """
    # nuclear charge
    qnuc = 0.0
    for Z,pos in atomlist:
        qnuc += Z

    if rho is None:
        qelec = 0.0

    else:
        # electronic charge

        # integrate electron density
        nelec = integral(atomlist, rho)
        # electronic charge
        qelec = -nelec

    # total charge
    qtot = qnuc + qelec

    if verbose > 0:
        print ""
        print "  Charges"
        print "  ======="
        print "  nuclear charge     qnuc  = %+e" % qnuc
        print "  electronic charge  qelec = %+e   (-1) x integrated electron density" % qelec
        print "  total charge       qtot  = %+e" % qtot
        print ""
    
    return qnuc, qelec, qtot

def center_of_charge(atomlist, rho,
                     verbose=1):
    """
    computes the centers of the nuclear charges and electronic 
    density:

                            /
         <r>      = 1/nelec | r * rho  (r)  dV
            elec            /       elec

                       
         <r>      = sum  Z  R   / sum  Z
            nuc        i  i  i       i  i

    The total charge density is defined as

         chrg  (r)    = - rho  (r)  +  sum  Z  delta(r - R )
             elec+nuc        elec         i  i            i

    so that the total center of charge takes the form
 
                        /                       / /
         <r>         =  | r * chrg   (r)   dV  /  | chrg  (r)  dV
            elec+nuc    /         elec+nuc    /   /     elec+nuc
    

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    rho          :  callable, rho(x,y,z) evaluates the electron density on a grid.
                    `rho` can be set to None if only the nuclear center of charge should be calculated.

    Optional
    ========
    verbose      :  if > 0, write additional output

    Returns
    =======
    center_nuc   :  center of nuclear charge, numpy array of shape (3,)
    center_elec  :  center of electronic charge density or None, if
                    the electron density is zero
    center_tot   :  center of total charge or None if the position
                    is undefined because the total charge (elec.+nuc.)
                    is zero
    """
    # charges
    qnuc, qelec, qtot = total_charge(atomlist, rho, verbose=0)
    # number of electrons
    nelec = -qelec
    
    # center of nuclear charge
    center_nuc = np.zeros(3)
    for Z,pos in atomlist:
        center_nuc += Z * np.array(pos)
    center_nuc /= qnuc

    if rho is None:
        # No electrons, so we only calculate the nuclear center
        return center_nuc, None, None
    
    # compute center of electronic charge density (!)
    # integrand for
    #       /                 /  /
    # <r> = | rho(r) r d^3r  /   | rho(r) d^3r
    #       /               /    /
    def integrand(x,y,z, xyz=0):
        r = [x,y,z]
        return rho(x,y,z)/nelec * r[xyz]

    if abs(nelec) > 1.0e-1:
        center_elec = np.zeros(3)
        for xyz in [0,1,2]:
            # integrals for <x>, <y> and <z>
            center_elec[xyz] = integral(atomlist,
                                        lambda x,y,z: integrand(x,y,z, xyz=xyz))
    else:
        # There are no electrons in this molecule,
        # so we cannot define the center of charge
        center_elec = None
            
    # center of total charge
    if abs(qtot) > 1.0e-1:
        # The total charge should be an integer, however because of
        # numerical integration qtot may deviate from 0, even if
        # qnuc - nelec = 0
        center_tot = (qelec * center_elec + qtot * center_nuc)/qtot
    else:
        # For neutral molecules the center of charge is undefined
        # because we would have to divide by 0.
        center_tot = None

    if verbose > 0:
        print " "
        print "  Center of Charge"
        print "  ================"
        print "  center of nuclear charge density  =  %+e  %+e  %+e" % tuple(center_nuc)
        if center_elec is None:
            print "  center of electronic charge is undefined because nelec = %e" % nelec
        else:
            print "  center of elec. charge density    =  %+e  %+e  %+e" % tuple(center_elec)
        if center_tot is None:
            print "  center of total charge is undefined because qtot = qelec + qnuc = %+e " % qtot
        else:
            print "  center of total charge            =  %+e  %+e  %+e" % tuple(center_tot)
        print " "
        
    return center_nuc, center_elec, center_tot
    

def effective_potential_func(atomlist, rho, xc,
                             nelec=np.inf,
                             nuclear=True):
    """
    create function for evaluating the effective Kohn-Sham potential

                            /  rho(r')
        Veff(r) = Vnuc(r) + | --------- dr'  +  Vxc[rho](r)
                            /  |r - r'|

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    rho          :  callable, rho(x,y,z) evaluates the electron density on a grid
    xc           :  instance of libXCFunctional or None

    Optional
    ========
    nelec        :  number of electrons, for nelec < 2, we don't have
                    to consider the electron-electron interaction (Veff=Vnuc)
    nuclear      :  if set to False, the nuclear potential is excluded from 
                    the effective potential, so Veff = Vcoul + Vxc

    Returns
    =======
    veff         :  callable, veff(x,y,z) evaluates effective potential on any grid
    """
    if xc is None:
        if nelec < 2:
            print "NOTE: One-electron systems are treated exactly (without xc-potential)!"
        elif nelec == 2:
            print "NOTE: Two-electron systems are treated with Hartree-Fock theory!"
        else:
            print "NOTE: Since no xc-functional was specified, many electron systems"
            print "      are treated with Hartree theory (no exchange, no correlation)!"
    else:
        print "NOTE: Many-electron systems are treated with density functional theory!"
        
    # Bring geometry data into a form understood by the module MolecularIntegrals
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    Nat = len(atomic_numbers)

    if nelec > 1 or (not (xc is None)):
        # solve Poisson equation to get classical Coulomb potential
        # due to electron density
        #
        coulomb_potential = multicenter_poisson(rho,
                                                atomic_coordinates, atomic_numbers,
                                                radial_grid_factor=settings.radial_grid_factor,
                                                lebedev_order=settings.lebedev_order)

    if not (xc is None):
        # GGA xc-functionals need the gradient squared of the density
        #          __     2
        # sigma = (\/ rho)
        sigma = multicenter_gradient2(rho,
                                      atomic_coordinates, atomic_numbers,
                                      radial_grid_factor=settings.radial_grid_factor,
                                      lebedev_order=settings.lebedev_order)

    if not nuclear:
        print "nuclear potential Vnuc is not included in effective potential, Veff = Vcoul + Vxc!"
        
    def effective_potential(x,y,z):
        if nuclear:
            # nuclear attraction
            # Vnuc is also called the external potential.
            Vnuc = 0.0*x        
            for i in range(0, Nat):
                Zi = atomic_numbers[i]
                Ri = atomic_coordinates[:,i]
                
                Vnuc += - Zi / np.sqrt( (x-Ri[0])**2 + (y-Ri[1])**2 + (z-Ri[2])**2 )
        else:
            Vnuc = 0.0*x    

        if xc is None:
            # If no xc-functional is provided, 1 and 2 and many-electron
            # systems are treated differently.
            if nelec < 2:
                # One-electron systems such as H2+ can be solved
                # exactly, Vcoul = Vxc = 0.
                Veff = Vnuc
            elif nelec == 2:
                # For 2-electron systems with closed-shell ground state
                # such as H2 or H3+, we use Hartree-Fock theory, where
                # Coulomb and exchange energies cancel partly.
                #
                Vcoul = coulomb_potential(x,y,z)
                Veff = Vnuc + 0.5*Vcoul
            else:
                # In the absence of an xc-functional, the electron-electron
                # interaction consists only of the Hartree term.
                Vcoul = coulomb_potential(x,y,z)
                Veff = Vnuc + Vcoul
        else:
            # If an xc-functional is provided, all systems are treated the same.
            # The Coulomb interaction consists of the electron-electron repulsion
            # (Vcoul) and exchange and correlation (Vxc), even if there is only
            # a single electron. 
            # Vcoul is also called the Hartree potential.
            Vcoul = coulomb_potential(x,y,z)
            # Input for exchange-correlation functional:
            #   density
            n = rho(x,y,z)
            #   gradient of density, (grad rho)^2
            s = sigma(x,y,z)
            #s = 0.0*x
            # The wrapper around libxc expects 1d arrays as input.
            #   exchange part
            Vx = xc.func_x.vxc(n.flatten(), s.flatten())
            #   correlation part
            Vc = xc.func_c.vxc(n.flatten(), s.flatten())
            # 
            Vxc = np.reshape(Vx + Vc, n.shape)

            Veff = Vnuc + Vcoul + Vxc

        return Veff

    return effective_potential

def nuclear_potential_func(atomlist):
    """
    create function for evaluating the nuclear potential

                           Z(I)
     V    (r) =  sum  - ----------
      nuc           I    |r-R(I)|

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions

    Returns
    =======
    vnuc         :  callable, veff(x,y,z) evaluates nuclear potential on any grid
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    Nat = len(atomic_numbers)
    
    def nuclear_potential(x,y,z):
        # nuclear attraction
        # Vnuc is also called the external potential.
        Vnuc = 0.0*x        
        for i in range(0, Nat):
            Zi = atomic_numbers[i]
            Ri = atomic_coordinates[:,i]
            
            Vnuc += - Zi / np.sqrt( (x-Ri[0])**2 + (y-Ri[1])**2 + (z-Ri[2])**2 )
        return Vnuc

    return nuclear_potential


def laplacian_func(atomlist, phi):
    """
    create function for evaluating the Laplacian of phi,
      __2
      \/  phi

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    phi          :  callable, phi(x,y,z) evaluates the wavefunctions 

    Returns
    =======
    lap          :  callable, lap(x,y,z) evaluates the Laplacian of phi
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # compute Laplacian on multicenter grid
    lap = multicenter_laplacian(phi,
                                atomic_coordinates, atomic_numbers,
                                radial_grid_factor=settings.radial_grid_factor,
                                lebedev_order=settings.lebedev_order)

    return lap
    
def energy(atomlist, bra, ket, potential):
    """
    compute matrix element of electronic energy

                     /    *        __2 
       <bra|H|ket> = | bra  ( -1/2 \/  + V(r) ) ket dr
                     /

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    bra, ket     :  callables, bra(x,y,z) and ket(x,y,z) evaluate
                    the wavefunctions 
    potential    :  callable, potential(x,y,z) evaluates effective potential

    Returns
    =======
    en           : float, energy matrix element
    """
    # compute Laplacian of ket
    lap_ket = laplacian_func(atomlist, ket)
    
    def energy_density(x,y,z):
        # evaluate bra and ket on grid
        bra_xyz = bra(x,y,z)
        ket_xyz = ket(x,y,z)
        # apply KS Hamiltonian to ket
        Hket = -0.5 * lap_ket(x,y,z) + potential(x,y,z)*ket_xyz
        return bra_xyz * Hket

    # integrate energy density
    en = integral(atomlist, energy_density)
    
    return en


def energy_ee(atomlist, bra, ket, potential_ee):
    """
    compute matrix element of electronic energy

                     /    *        __2 
       <bra|H|ket> = | bra  ( -1/2 \/  + V(r) ) ket dr
                     /

    Unlike in `energy(...)` the kinetic energy and the nuclear attraction
    energy are computed together to avoid numerical errors.

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    bra, ket     :  callables, bra(x,y,z) and ket(x,y,z) evaluate
                    the wavefunctions 
    potential_ee :  callable potential_ee(x,y,z) evaluates the effective potential 
                    WITHOUT the nuclear potential, only Vcoul + Vxc. The nuclear potential
                    is constructed from the atomic coordinates and numbers.

    Returns
    =======
    en           : float, energy matrix element
    """
    # apply Hamiltonian to ket H\ket
    ham_ket = residual_ee_func(atomlist, ket, potential_ee, 0.0)
    
    def energy_density(x,y,z):
        # evaluate bra and H.ket on grid
        bra_xyz = bra(x,y,z)
        Hket_xyz = ham_ket(x,y,z)
        #  bra(x,y,z) H ket(x,y,z)
        return bra_xyz * Hket_xyz

    # integrate energy density
    en = integral(atomlist, energy_density)
    
    return en

def residual_func(atomlist, phi, potential, e):
    """
    create function for evaluating the "residual" (see Ref. [3])
                __2
       R = -1/2 \/  phi(r)  +  (V-e) phi(r)

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    phi          :  callable, phi(x,y,z) evaluates the approximate wavefunction
    potential    :  callable, potential(x,y,z) evaluates effective potential
    e            :  float, approximate energy of phi

    Returns
    =======
    residual     :  callable, residual(x,y,z) evaluates R on a grid
    """
    # compute Laplacian of phi
    lap_phi = laplacian_func(atomlist, phi)

    def residual(x,y,z):
        R = -0.5 * lap_phi(x,y,z) + (potential(x,y,z) - e)*phi(x,y,z)
        return R

    return residual

def residual_ee_func(atomlist, phi, potential_ee, e):
    """
    create function for evaluating the "residual" (see Ref. [3])
                __2
       R = -1/2 \/  phi(r)  +  (V-e) phi(r)

    Unlike in `residual_func(...)` the kinetic energy and the nuclear attraction
    energy are computed together to avoid numerical errors.

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    phi          :  callable, phi(x,y,z) evaluates the approximate wavefunction
    potential_ee :  callable potential_ee(x,y,z) evaluates the effective potential 
                    WITHOUT the nuclear potential, only Vcoul + Vxc. The nuclear potential
                    is constructed from the atomic coordinates and numbers.
    e            :  float, approximate energy of phi

    Returns
    =======
    residual     :  callable, residual(x,y,z) evaluates R on a grid
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    
    residual = multicenter_residual(phi, potential_ee, e,
                                    atomic_coordinates, atomic_numbers,
                                    radial_grid_factor=settings.radial_grid_factor,
                                    lebedev_order=settings.lebedev_order)
    
    return residual


def spherical_average_residual_func(atom, residual):
    """
    average the square of the residual R = (H-E)phi over 
    all angles so that only a function of the radius remains

               /
      avg(r) = | dOmega |R(r,Omega)|^2 
               /

    Parameters
    ==========
    atom        :   tuple (Zat,(X,Y,Z)) with the coordinates X,Y,Z of
                    the center around which the residual should be spherically averaged.
                    The atomic number Zat is used to select the radial grid.
    residual    :   callable, residual(x,y,z) evaluates the deviation from an exact solution

    Returns
    =======
    avg         :   callable, avg(r) evaluates the average error at radius r
    """
    def residual2(x,y,z):
        return residual(x,y,z)**2

    avg = spherical_average_func(atom, residual2,
                                 radial_grid_factor=settings.radial_grid_factor,
                                 lebedev_order=settings.lebedev_order)

    return avg


def energy_correction(atomlist, residual, phi,
                      method="R2"):
    """
    estimate the energy correction from the approximate solution phi 
    and the residual R.

    PROBLEM: 
    According to eqn. (11) in Ref. [3] the energy correction should
    estimated as:
                /            /
      delta_e = | R phi dV = | phi (H-e) phi dV
                /            /
    However, with this definition we will always get delta_e = 0 
    if e is the expectation value <phi|H|phi>.

    SOLUTION:
    Starting from Schroedinger's equation

      H (phi + dphi) = (e + de) (phi + dphi)
      (H-e)( phi + dphi) = de (phi + dphi)

    applying (H-e) again,

      (H-e)^2 (phi + dphi) = (de)^2 (phi + dphi)

    and multiplying from the left by (phi+dphi)^*, itegrating and 
    using the orthogonality of the eigenfunctions we get
                          2        2
      ||(H-e) (phi+dphi)||   = (de)

    After neglecting terms of order dphi, we get a usable approximation
    for the energy correction:

      delta_e = ||(H-e)phi|| = norm of the residual

    So
                          /
      delta_e = +/- sqrt{ | R^2 dV }
                          /

    Since we are interested in the solution with the lowest energy we
    choose delta_e = - ||R||. 
    """
    if method == "Becke":
        # energy correction as proposed in Becke's paper

        def integrand(x,y,z):
            return residual(x,y,z) * phi(x,y,z)

        delta_e = integral(atomlist, integrand)

    else:
        # my energy correction
        
        def integrand(x,y,z):
            return abs(residual(x,y,z))**2

        nrm2_residual = integral(atomlist, integrand)

        delta_e = -np.sqrt(nrm2_residual)
    
    return delta_e


def source_func(delta_e, phi, residual):
    """ 
    create function for evaluating the source term 

         delta_e phi - R

    in the inhomogeneous Schroedinger equation (see eqn. (8) in Ref. [3])
    """
    
    def source(x,y,z):
        return delta_e * phi(x,y,z) - residual(x,y,z)

    return source

def orbital_correction(atomlist, potential, source,  e, delta_e):
    """
    determine the orbital correction  delta_phi  by solving the
    inhomogeneous Schroedinger equation 
           __2
      -1/2 \/   (delta_phi) + (V - e - delta_e) delta_phi = delta_e phi - R

    approximately for a potential that is spherically averaged
    around each atom.

    PROBLEM:
    The inhomogeneous Schroedinger equation

         H (phi + dphi) = (e + de) (phi + dphi)
    or
         (H-e-de) dphi = -(H-e)phi + de phi

    has the trivial solution dphi = -phi, which is not what we want.
    How can we avoid this?

    SOLUTION:
    One way consists in modifying the inhomogeneous Schroedinger equation
    in such a way, that the dphi = -phi is not a solution anymore. This
    can be achieved by removing de from the left hand side of the equation
    and solving 

         (H-e) dphi = -(H-e) phi + de phi

    instead.
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # Left hand side is (H-e-0*delta_e) dphi = (H-e) dphi
    delta_phi = multicenter_inhomogeneous_schroedinger(potential, source,  e+0*delta_e,
                                                       atomic_coordinates, atomic_numbers,
                                                       radial_grid_factor=settings.radial_grid_factor,
                                                       lebedev_order=settings.lebedev_order)
    return delta_phi

def inhomogeneous_schroedinger(atomlist, potential, source, E):
    """
    solve the inhomogeneous Schroedinger equation

           __2
      -1/2 \/   phi + (V(r) - E) phi = S(r)

    for the wavefunction (correction) phi. 

    In order to decouple different angular momenta the equation is not 
    solved for the full potential V but for a simplified potential

            (sph)                   /
           V       =  sum  1/(4 pi) | (w * V )(r-R ) dOmega
                         I          /   I         I         I

    which is the spherical average of the potential around each atom I.
    R_I is the position of the atom and w_I is the weight function 
    of the fuzzy Voronoi decomposition.


    Parameters
    ==========
    atomlist     :  list of tuples (Zat,(x,y,z)) with the atomic
                    positions that define the multicenter grid
    potential    :  callable, potential(x,y,z) evaluates the V on
                    on a grid
    source       :  callable, source(x,y,z) evaluates the source term 
                    S on the right hand side
    E            :  float, energy of solution

    Returns
    =======
    phi          :  callable, wavefunction phi(x,y,z)
    
    """
    # Translate molecular geometry into the format understood
    # by the multicenter integration code
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    print "Schroedinger equation ..."
    dphi = multicenter_inhomogeneous_schroedinger(potential, source,  E,
                                                  atomic_coordinates, atomic_numbers,
                                                  radial_grid_factor=settings.radial_grid_factor,
                                                  lebedev_order=settings.lebedev_order)

    return dphi

def variational_mixture(atomlist, phi, delta_phi, potential):
    """
    Since the orbital correction is only approximate due to the 
    spherical averaging, for the improved wavefunction we make
    the ansatz
              
          psi(r)  = a phi(r) + b delta_phi(r)

    where the coefficients a and b are optimized variationally.
    I understand this to mean that we solve the generalized 
    eigenvalue problem

      ( Haa Hab ) (a)     ( Saa Sab )
      (         ) ( ) = E (         )
      ( Hba Hbb ) (b)     ( Sab Sbb )

    where Haa = <phi|H|phi>, Hab = <phi|H|delta_phi>, etc.

    Returns
    =======
    a,b      :  floats, mixing coefficients
    """
    # First we need some matrix elemens
    Saa = overlap(atomlist, phi,phi)
    Sab = overlap(atomlist, phi,delta_phi)
    Sbb = overlap(atomlist, delta_phi, delta_phi)

    Haa = energy(atomlist, phi, phi, potential)
    Hab = energy(atomlist, phi, delta_phi, potential)
    Hbb = energy(atomlist, delta_phi, delta_phi, potential)

    H = np.array([
        [Haa, Hab],
        [Hab, Hbb]])
    S = np.array([
        [Saa, Sab],
        [Sab, Sbb]])

    eigvals, eigvecs = sla.eigh(H,S)

    # Choose the solution with the lower energy
    a,b = eigvecs[:,0]

    """
    # Of the two possible solutions choose the one where
    # abs(a) > abs(b) so that the b*delta_phi is a small
    # correction to a*phi.
    
    # The chosen solution has the smaller ratio between |b| and |a|
    ratio = abs(eigvecs[1,:]) / abs(eigvecs[0,:])
    isel = np.argmin(ratio)
    a,b = eigvecs[:,isel]
    """

    # Sign of an eigenvector is arbitrary. We choose
    # the sign such that `a` is always positive.
    if a < 0:
        a *= -1
        b *= -1

    print "Eigenenergies= %s" % eigvals
    #print "selected solution isel= %d" % (isel+1)
    print "Coefficients for mixing old orbital and orbital correction"
    print "   phi^(new) = a*phi^(old) + b*delta_phi "
    print "a = %s" % a
    print "b = %s" % b
        
    return a,b

def add_two_functions(atomlist, f,g, a,b):
    """
    add two functions 
    
       h(r) = a*f(r) + b*g(r)
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    h = multicenter_operation([f,g], lambda fs: a*fs[0]+b*fs[1],
                              atomic_coordinates, atomic_numbers,
                              radial_grid_factor=settings.radial_grid_factor,
                              lebedev_order=settings.lebedev_order)
    return h

def coulomb_energy(atomlist, rho):
    """
    compute classical electron-electron repulsion energy 

                   /  /   rho(r1) rho(r2)
       Ecoul = 1/2 |d1|d2 ---------------
                   /  /      |r1-r2|

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    rho          :  callable, rho(x,y,z) evaluates the electron density on a grid
    
    Returns
    =======
    Ecoul        :  float, electrostatic interaction of density with itself
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # solve Poisson equation to get classical Coulomb potential
    # due to electron density
    coulomb_potential = multicenter_poisson(rho,
                                            atomic_coordinates, atomic_numbers,
                                            radial_grid_factor=settings.radial_grid_factor,
                                            lebedev_order=settings.lebedev_order)

    # integrate Vcoul(r)*rho(r)
    def energy_density(x,y,z):
        return coulomb_potential(x,y,z) * rho(x,y,z)
    
    Ecoul = 0.5 * integral(atomlist, energy_density)
    
    return Ecoul


def exchange_correlation_energy(atomlist, rho, xc):
    """
    compute the exchange-correlation energy for density 'rho'
                    /
         E  [rho] = | e  [rho](r) rho(r) dr
          xc        /  xc

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    rho          :  callable, rho(x,y,z) evaluates the electron density on a grid
    xc           :  instance of libXCFunctional or None
    
    Returns
    =======
    Exc        :  float, exchange-correlation energy

    """
    if xc is None:
        return 0.0
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # GGA xc-functionals need the gradient squared of the density
    #          __     2
    # sigma = (\/ rho)
    sigma = multicenter_gradient2(rho,
                                  atomic_coordinates, atomic_numbers,
                                  radial_grid_factor=settings.radial_grid_factor,
                                  lebedev_order=settings.lebedev_order)

    def energy_density(x,y,z):
        # compute the energy density  e_xc[rho](x,y,z) * rho(x,y,z)
        
        # Input for exchange-correlation functional:
        #   density
        n = rho(x,y,z)
        #   gradient of density, (grad rho)^2
        # NOTE: The computation of sigma = (grad rho)^2 is tricky,
        #       probably because at the nuclei the density has cusps, which
        #       makes numerical differentiation very difficult.
        s = sigma(x,y,z)
        #s = 0.0*x
        # The wrapper around libxc expects 1d arrays as input.
        #   exchange part
        ex = xc.func_x.exc(n.flatten(), s.flatten())
        #   correlation part
        ec = xc.func_c.exc(n.flatten(), s.flatten())
        # 
        exc = np.reshape(ex + ec, n.shape)

        return exc * n

    Exc = integral(atomlist, energy_density)

    return Exc

    
def total_dft_energy(atomlist, orbitals, xc,
                     nelec=0):
    """
    get the total DFT energy
                                   
      E_tot =  - 1/2 sum_i <i|T|i>             kinetic energy

                  /
               +  | rho(r) V(r)                nuclear attraction
                  /         nuc

                     /  /   rho(1) rho(2)
               + 1/2 |d1|d2 -------------      classical electron-electron repulsion
                     /  /       r12

               + Exc[rho]                      xc-energy

               + Erep                          nuclear-nuclear repulsion

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    orbitals     :  doubly occupied orbitals, list of callables, 
                    orbitals[i](x,y,z) evaluates the wavefunction of the i-th orbital
    xc           :  instance of libXCFunctional or None
    
    Optional
    ========
    nelec        :  number of electrons, if nelec==1, the interaction of a single electron
                    with itself is excluded
    
    Returns
    =======
    Etot         :  float, total energy
    """   
    # kinetic energy of orbitals
    Ekin = 0.0
    # electron-nuclear attraction 
    Enuc = 0.0
    for orb in orbitals:
        # Each orbital is double occupied and
        # spin-up and spin-down orbitals have the same integrals
        Ekin += 2*kinetic(atomlist, orb, orb)
        Enuc += 2*nuclear(atomlist, orb, orb)

    if xc is None:
        # If no xc-functional is provided 1 electron systems are treated
        # exactly, 2 electron system with HF theory and many-electron systems
        # with Hartree theory (neither exchange nor correlation)
        if nelec == 1:
            # For one-electron systems there is no electron-electron interaction
            Ecoul = 0.0
            Exc = 0.0
        
            Eelec = Ekin + Enuc
        elif nelec == 2:
            # electron-electron interaction
            rho = density_func(orbitals)
            # classical Coulomb energy
            Ecoul = coulomb_energy(atomlist, rho)
            #
            # Since rho = 2*|phi|^2, the Coulomb energy contains the interaction
            # of each electron with itself:
            #
            #              /  /   (2*phi(r1)^2) (2*phi(r2)^2)
            #  Ecoul = 1/2 |d1|d2 --------------------------- = 2 * (phi,phi|phi,phi) 
            #              /  /        |r1 - r2|
            #
            #        = 2 * J_11 = 2 * K_11
            
            # For two-electron closed-shell systems, the
            # self-interaction energy cancels with the exchange energy:
            #   2 * J_11 - K_11 = J_11
            Exc = -0.5 * Ecoul
            
            # Eelec = Ekin + Enuc + 0.5*Ecoul
            Eelec = Ekin + Enuc + Ecoul + Exc
        else:
            # no exchange or correlation
            Eelec = Ekin + Enuc + Ecoul
    else:
        # electron-electron interaction
        rho = density_func(orbitals)
        # classical Coulomb energy
        Ecoul = coulomb_energy(atomlist, rho)
        # exchange and correlation energy
        Exc = exchange_correlation_energy(atomlist, rho, xc)
        # total electronic energy
        Eelec = Ekin + Enuc + Ecoul + Exc

    # nuclear repulsion
    Erep = nuclear_repulsion(atomlist)
        
    Etot = Eelec + Erep

    print "   DFT Energies"
    print "   ------------"
    print "kinetic energy                 Ekin  = %e Hartree" % Ekin
    print "nuclear attraction energy      Enuc  = %e Hartree" % Enuc
    print "electron-electron repulsion    Ecoul = %e Hartree" % Ecoul
    print "exchange-correlation energy    Exc   = %e Hartree" % Exc
    print ""
    print "total electronic energy        Eelec = %e Hartree" % Eelec
    print "nuclear repulsion energy       Erep  = %e Hartree" % Erep
    print "virial -V/T = %8.5f" % (-(Etot-Ekin)/Ekin)
    print ""
    print "total energy                   Etot  = %e Hartree" % Etot
        
    return Etot

def overlap_matrix(atomlist, orbitals):
    """
    compute overlap matrix between a set of orbitals
    """
    norb = len(orbitals)
    # overlap matrix
    S = np.zeros((norb,norb), dtype=complex)
    for i in range(0, norb):
        for j in range(i, norb):
            S[i,j] = overlap(atomlist, orbitals[i], orbitals[j])
            S[j,i] = S[i,j]

    #
    print "overlap matrix S"
    labels = ["orb.  %3.1d" % (i+1) for i in range(0, norb)]
    print annotated_matrix(S, labels, labels)

    return S

def hamiltonian_matrix(atomlist, orbitals, potential):
    """
    compute matrix elements of Hamiltonian (Kohn-Sham) between a set of orbitals
    """
    norb = len(orbitals)
    H = np.zeros((norb,norb))
    for i in range(0, norb):
        for j in range(i, norb):
            H[i,j] = energy(atomlist, orbitals[i], orbitals[j], potential)
            # assume that Hamiltonian is symmetric
            H[j,i] = H[i,j]

    #
    print "Hamiltonian matrix H"
    labels = ["orb.  %3.1d" % (i+1) for i in range(0, norb)]
    print annotated_matrix(H, labels, labels)

    return H

def hamiltonian_matrix_ee(atomlist, orbitals, potential_ee):
    """
    compute matrix elements of Hamiltonian (Kohn-Sham) between a set of orbitals

    The potential contains only the electron-electron interaction but not the
    interaction between nuclei and elelctrons which is automatically added.
    """
    norb = len(orbitals)
    H = np.zeros((norb,norb))
    for i in range(0, norb):
        for j in range(i, norb):
            H[i,j] = energy_ee(atomlist, orbitals[i], orbitals[j], potential_ee)
            # assume that Hamiltonian is symmetric
            H[j,i] = H[i,j]

    #
    print "Hamiltonian matrix H"
    labels = ["orb.  %3.1d" % (i+1) for i in range(0, norb)]
    print annotated_matrix(H, labels, labels)

    return H

def orbital_transformation(atomlist, orbitals, X):
    """
    perform a linear orbital transformation, the new orbitals phi'
    are a linear combination of the old orbitals phi


      phi' = sum  X    phi
         i      j  j,i    j

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    orbitals     :  old orbitals, list of callables, 
                    orbitals[i](x,y,z) evaluates the wavefunction of the i-th orbital
                    on a grid
    X            :  matrix defining linear transformation, 
                    with shape (norb,norb)

    Returns
    =======
    orbitals_transf  :  new transformed orbitals, list of callables
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    norb = len(orbitals)

    orbitals_transf = []
    for i in range(0, norb):
        # i-th orbital in the transformed set
        phi_transf = multicenter_operation(orbitals,
                                          lambda orbs: sum( [X[j,i] * orbs[j] for j in range(0, norb)] ),
                                          atomic_coordinates, atomic_numbers,
                                          radial_grid_factor=settings.radial_grid_factor,
                                          lebedev_order=settings.lebedev_order)
        orbitals_transf.append( phi_transf )

    return orbitals_transf

def solve_matrix_schroedinger(atomlist, orbitals, potential):
    """
    solve the generalized eigenvalue problem

       H.C = E.S.C

    for the molecular orbital coefficients C. The Hamiltonian and
    overlap matrix H and S are computed in the basis set defined by the 
    list of `orbitals`. 

    Parameters
    ==========
    atomlist  :  list of tuples (Z,[x,y,z]) with atomic numbers
                 and positions    
    orbitals  :  orbitals, list of callables 
    potential :  callable, potential(x,y,z) evaluates the effective potential 

    Returns
    =======
    orbe     :  energies of occupied orbitals
    orbitals :  occupied orbitals, list of callables 
    """
    S = overlap_matrix(atomlist, orbitals)
    H = hamiltonian_matrix(atomlist, orbitals, potential)

    # solve generalized eigenvalue problem for eigenvalues E[i]
    # and eigenvectors C[j,i]
    E, C = sla.eigh(H,S)

    # wavefunctions of solutions are linear combinations of input orbitals:
    #   orbitals_eig[i](x,y,z) = sum_j C[j,i] orbitals[j](x,y,z)
    orbitals_eig = orbital_transformation(atomlist, orbitals, C) 

    # show coefficients of new orbitals MO' in the basis of the old orbitals MO .
    print "  Orbital coefficients"
    print "  ===================="
    print annotated_matrix(np.vstack((E, C)),\
                ["en. (a.u.)", "-"] + ["MO %s" % (i+1) for i in range(0, len(E))],\
                                      ["MO'%s" % (i+1) for i in range(0, len(E))], colwidth=10)
    
    return E, orbitals_eig

def solve_matrix_schroedinger_ee(atomlist, orbitals, potential_ee):
    """
    solve the generalized eigenvalue problem

       H.C = E.S.C

    for the molecular orbital coefficients C. The Hamiltonian and
    overlap matrix H and S are computed in the basis set defined by the 
    list of `orbitals`. 

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions    
    orbitals     :  orbitals, list of callables 
    potential_ee :  callable, potential_ee(x,y,z) evaluates the effective potential 
                    WITHOUT the nuclear attraction potential, which is added automatically

    Returns
    =======
    orbe     :  energies of occupied orbitals
    orbitals :  occupied orbitals, list of callables 
    """
    S = overlap_matrix(atomlist, orbitals)
    H = hamiltonian_matrix_ee(atomlist, orbitals, potential_ee)

    # solve generalized eigenvalue problem for eigenvalues E[i]
    # and eigenvectors C[j,i]
    E, C = sla.eigh(H,S)

    # wavefunctions of solutions are linear combinations of input orbitals:
    #   orbitals_eig[i](x,y,z) = sum_j C[j,i] orbitals[j](x,y,z)
    orbitals_eig = orbital_transformation(atomlist, orbitals, C) 

    # show coefficients of new orbitals MO' in the basis of the old orbitals MO .
    print "  Orbital coefficients"
    print "  ===================="
    print annotated_matrix(np.vstack((E, C)),\
                ["en. (a.u.)", "-"] + ["MO %s" % (i+1) for i in range(0, len(E))],\
                                      ["MO'%s" % (i+1) for i in range(0, len(E))], colwidth=10)
    
    return E, orbitals_eig


def orthogonalize_orbitals(atomlist, orbitals):
    """
    orthogonalize a set of orbitals using Loewdin's symmetric transformation

      phi' = sum  X    phi
         i      j  j,i    j

    with X = S^{-1/2}

    see Szabo & Ostlund section 3.4.5 
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    norb = len(orbitals)

    # perform Loewdin orthogonalization, X = S^{-1/2}
    print "BEFORE ORTHOGONALIZATION"
    S = overlap_matrix(atomlist, orbitals)
    X = sla.sqrtm(sla.inv(S))
    
    orbitals_ortho = orbital_transformation(atomlist, orbitals, X)

    # check orthogonalization
    print "AFTER ORTHOGONALIZATION"
    S = overlap_matrix(atomlist, orbitals_ortho)
    
    return orbitals_ortho

#
# tight-binding DFT provides the initial guesses for
# the orbitals, which are refined using Becke's algorithm
# of Ref. [3]
#

class BasissetFreeDFT(object):
    def __init__(self, atomlist, xc, charge=0):
        """
        Parameters
        ==========
        atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                        and positions
        xc           :  exchange-correlation functional, instance of libXCFunctional or None

        Optional
        ========
        charge       :  total charge of molecule
        """
        self.atomlist = atomlist
        self.xc = xc
        if not (xc is None):
            # Show some information about the XC functional
            desc_x = xc.func_x.description()
            desc_c = xc.func_c.description()
            print "XC-Functional    "
            print "-------------  \n"
            print desc_x
            print desc_c
            """
            ### I think this issue has been fixed, we can compute (grad n)^2 now.
            #
            # Only LDA functionals are allowed so far
            if not (("XC_FAMILY_LDA" in desc_x) and ("XC_FAMILY_LDA" in desc_c)):
                msg  = "\n\nERROR: Only LDA functionals can be used so far!\n"
                msg += "  This is because the density gradient sigma = (grad n)^2 \n"
                msg += "  needed for GGA functionals has delta-like peaks at \n"
                msg += "  the atomic positions due to the cusp condition. I don't know \n"
                msg += "  yet how to remove this discontinuity of the density gradient\n"
                msg += "  at the cusps. \n"
                raise ValueError(msg)
            """
        self.charge = charge
    def getOrbitalGuess(self):
        """
        guess occupied valence orbitals by a DFTB calculation, 
        core orbitals are taken from the confined pseudoatoms.

        Returns
        =======
        orbitals  :   list of doubly occupied DFTB orbitals (core + valence),
                      orbitals[i](x,y,z) evaluates the i-th orbital
                      on a grid
        """
        from DFTB.PES import PotentialEnergySurfaces
        from DFTB import XYZ
        from DFTB.BasisSets import AtomicBasisSet, LCAOBasisSet
        # initialize DFTB calculator
        pes = PotentialEnergySurfaces(self.atomlist, charge=self.charge)

        # basis set of atomic core orbitals
        basisAO_core = AtomicBasisSet(self.atomlist, orbital_set="core")
        core_orbitals = basisAO_core.bfs

        # The guess for the valence orbitals is taken from
        # a DFTB calculation, if there are occupied valence orbitals.
        if pes.tddftb.dftb2.Nelec_val > 0:
            # compute ground state electronic structure at DFTB level
            x = XYZ.atomlist2vector(self.atomlist)
            pes.getEnergy_S0(x)
            # Kohn-Sham energies and MO coefficients
            orbe = pes.tddftb.dftb2.getKSEnergies()
            orbs = pes.tddftb.dftb2.getKSCoefficients()
            
            # indices of occupied orbitals
            f = pes.tddftb.dftb2.getOccupation()
            occupied_orbs = np.where(f > 0.1)[0]
            # energies and MOs of occupied orbitals
            orbe_occ = orbe[occupied_orbs]
            orbs_occ = orbs[:,occupied_orbs]

            # basis set of atomic valence orbitals 
            basisAO_val = AtomicBasisSet(self.atomlist, orbital_set="valence")
            # basis set consisting of the occupied molecular orbitals
            basisMO_val = LCAOBasisSet(basisAO_val.bfs, orbs_occ)
            valence_orbitals = basisMO_val.bfs

        else:
            # There are no occupied valence orbitals, this happens
            # for instance in the Li^+ cation.
            valence_orbitals = []

        print "number of core orbitals    : %d" % len(core_orbitals)
        print "number of valence orbitals : %d" % len(valence_orbitals)
        
        orbitals = core_orbitals + valence_orbitals
        
        return orbitals

    def solveKohnSham(self, max_iter=300, thresh=1.0e-4):
        """
        solve the Kohn-Sham equations iteratively on the multicenter grid

        Suppose we have an approximate solution phi with energy e
        that differs from the sought-after wavefuntion psi by dphi.

               H psi = E psi

        After making the ansatz

          psi = phi + dphi
          E   = e + de

        the Schroedinger equation turns into an inhomogeneous equation for the
        orbital correction dphi:

           (H - e - de) dphi = -(H-e) phi + de phi       (*)

        On the right hand side we define the residual R = (H-e) phi that tells
        us how much the solution phi deviates from fulfilling the Schroedinger eqn.

        In addition to the orbital correction, the inhomogeneous equation (*) 
        has the solution dphi = -phi, provided that phi is very close to the
        true solution:

              dphi = -phi  =>  (H - e)(-phi) + de phi = de phi - (H-e)(-phi)
                           =>  (H - e) phi = 0
                           =>  phi = psi

        The solution dphi = -phi is equivalent to the solution dphi = 0, which
        indicates convergence. Therefore we need to stop the orbital correction
        cycle if either
           (1) |R| = |(H-e) phi| = 0 
        or 
           (2) |dphi| = 0


        Optional
        ========
        max_iter    :   maximum number of SCF cycles, 
                        the solution of the Kohn-Sham equations and the self-consistency 
                        of the density are intertwined
        thresh      :   threshold for SCF convergence, the iteration is stopped if the
                        sum of the orbital corrections drops below this threshold.

        Returns
        =======
        Etot        :   total DFT energy
        orbitals    :   doubly occupied orbitals, list of callables, orbital[i](x,y,z)
        energies    :   list of Kohn-Sham orbital energies
        """
        print "initial orbital guess from DFTB calculation"
        orbitals = self.getOrbitalGuess()
        
        norb = len(orbitals)
        # all orbitals are doubly occupied
        nelec = 2*norb

        orbital_energies = np.zeros(norb)

        print "solve Kohn-Sham equations on a multicenter grid"

        print_grid_summary(self.atomlist,
                           settings.lebedev_order, settings.radial_grid_factor)
        
        for k in range(0, max_iter):
            print "Iteration k= %3.1d" % k
                       
            rho = density_func(orbitals)
            veff_new = effective_potential_func(self.atomlist, rho, self.xc, nelec=nelec)
            if k > 0:
                if total_error < 1.0:
                    # The density is only updated if the orbitals have more
                    # or less stabilized.
                    print "updating effective potential veff[rho^(new)]"

                    #
                    #  (next)        (old)        (new)
                    # V       = 2/3 V      + 1/3 V
                    veff = add_two_functions(self.atomlist,
                                             veff, veff_new, 2.0/3.0, 1.0/3.0)
            else:
                veff = veff_new
            
            total_error = 0.0

            orbital_corrections = []
            for i,phi in enumerate(orbitals):
                # compute corrections for each orbital
            
                # energy expectation value <phi|H|phi>
                e = energy(self.atomlist, phi, phi, veff)
                
                residual = residual_func(self.atomlist, phi, veff, e)
                delta_e = energy_correction(self.atomlist, residual, phi, method="Becke")
                source = source_func(delta_e, phi, residual)

                #
                nrm2_residual = overlap(self.atomlist, residual, residual)
                print "residual norm for orbital %d  |(H-E)phi|^2= %e" % (i+1, nrm2_residual)

                if nrm2_residual < thresh**2:
                    # If the residual norm is converged we don't have to
                    # compute the orbital correction for this orbital. 
                    print "Residual norm for orbital %d CONVERGED" % (i+1)
                    continue
                
                # orbital correction
                delta_phi = orbital_correction(self.atomlist, veff, source,  e, delta_e)
                nrm2_delta_phi = overlap(self.atomlist, delta_phi, delta_phi)
                total_error += nrm2_delta_phi
                
                print "    orbital i= %d    energy E= %e   delta E= %+e   |dphi|^2= %e" % (i+1, e, delta_e, nrm2_delta_phi)
                
                orbital_corrections.append(delta_phi)

                """
                ### DEBUG
                import matplotlib.pyplot as plt
                r = np.linspace(0.00001, 10.0, 5000)
                
                plt.plot(r, veff(r,0*r,0*r), ls="-.", lw=2, label=r"$v_{eff}$")

                plt.plot(r, residual(r,0*r,0*r), label=r"residual $(H-\epsilon)\phi$")
                plt.plot(r, source(r,0*r,0*r), label=r"source $\Delta \epsilon \phi - (H-E)\phi$")
                # wavefunctions
                plt.plot(r, phi(r,0*r,0*r), label=r"guess $\phi$")
                plt.plot(r, delta_phi(r,0*r,0*r), label=r"correction $\Delta \phi$")
                plt.plot(r, phi(r,0*r,0*r) + delta_phi(r,0*r,0*r), label=r"$\phi + \Delta \phi$")
                plt.plot(r, phi(r,0*r,0*r) - delta_phi(r,0*r,0*r), label=r"$\phi - \Delta \phi$")
                
                plt.legend()
                plt.show()
                ###
                """
                
            print "k= %3.1d  error = %e (threshold = %e)" % (k, total_error, thresh)

            # Since the Kohn-Sham equations are solved separately for
            # each orbital, the orthogonality may be lost, if one orbital
            # collapses to a lower-energy solution. 


            # Now we solve the Kohn-Sham equations in the basis of the current
            # orbitals and orbitals corrections {phi,delta phi} (at most 2*norb basis functions)
            # to restore orthogonality and improve the orbitals.
            orbital_energies, orbitals = solve_matrix_schroedinger(self.atomlist,
                                                                   orbitals+orbital_corrections, veff)
            # Only occupied orbitals are needed
            orbitals = orbitals[:norb]
            orbital_energies = orbital_energies[:norb]

            Etot = total_dft_energy(self.atomlist, orbitals, self.xc,
                                    nelec=nelec)
            print "total DFT energy Etot= %e" % Etot
            
            if total_error < thresh:
                print "SCF cycle CONVERGED"
                break
        else:
            raise RuntimeError("SCF did not converge in '%s' iterations!" % max_iter)

        return Etot, orbitals, orbital_energies

    def solveKohnSham2(self, max_iter=300, thresh=1.0e-4):
        """
        solve the Kohn-Sham equations iteratively on the multicenter grid

        Optional
        ========
        max_iter    :   maximum number of SCF cycles, 
                        the solution of the Kohn-Sham equations and the self-consistency 
                        of the density are intertwined
        thresh      :   threshold for SCF convergence, the iteration is stopped if the
                        sum of the orbital corrections drops below this threshold.

        Returns
        =======
        Etot        :   total DFT energy
        orbitals    :   doubly occupied orbitals, list of callables, orbital[i](x,y,z)
        energies    :   list of Kohn-Sham orbital energies
        """
        print "initial orbital guess from DFTB calculation"
        orbitals = self.getOrbitalGuess()
        
        norb = len(orbitals)
        # all orbitals are doubly occupied
        nelec = 2*norb

        orbital_energies = np.zeros(norb)

        print "solve Kohn-Sham equations on a multicenter grid"
        
        for k in range(0, max_iter):
            print "Iteration k= %3.1d" % k
                       
            rho = density_func(orbitals)
            veff_new = effective_potential_func(self.atomlist, rho, self.xc, nelec=nelec)
            if k > 0:
                print "updating effective potential veff[rho^(new)]"

                #
                #  (next)        (old)        (new)
                # V       = 2/3 V      + 1/3 V
                veff = add_two_functions(self.atomlist,
                                         veff, veff_new, 2.0/3.0, 1.0/3.0)
            else:
                veff = veff_new
            
            total_error = 0.0

            for i,phi in enumerate(orbitals):
                # compute corrections for each orbital
            
                # energy expectation value <phi|H|phi>
                e = energy(self.atomlist, phi, phi, veff)
                orbital_energies[i] = e
                
                residual = residual_func(self.atomlist, phi, veff, e)
                delta_e = energy_correction(self.atomlist, residual, phi)
                source = source_func(delta_e, phi, residual)
                
                # orbital correction
                delta_phi = orbital_correction(self.atomlist, veff, source,  e, delta_e)
                nrm2_delta_phi = overlap(self.atomlist, delta_phi, delta_phi)
                total_error += nrm2_delta_phi

                # update orbital by mixing with orbital correction
                a,b = variational_mixture(self.atomlist, phi, delta_phi, veff)
                # improved orbital
                #    (new)        (old)              
                # phi      = a phi        +  b delta_phi
                #
                phi = add_two_functions(self.atomlist, phi, delta_phi, a, b)
                orbitals[i] = phi

                # check normalization
                norm2 = overlap(self.atomlist, phi, phi)
                print "<phi|phi>= %e" % norm2
                

                
                print "    orbital i= %d    energy E= %e   delta E= %+e   |dphi|^2= %e" % (i+1, e, delta_e, nrm2_delta_phi)
                            
            print "k= %3.1d  error = %e (threshold = %e)" % (k, total_error, thresh)

            # Since the Kohn-Sham equations are solved separately for
            # each orbital, the orthogonality may be lost, if one orbital
            # collapses to a lower-energy solution. 

            #orbitals = orthogonalize_orbitals(self.atomlist, orbitals)
                        
            Etot = total_dft_energy(self.atomlist, orbitals, self.xc,
                                    nelec=nelec)
            print "total DFT energy Etot= %e" % Etot
            
            if total_error < thresh:
                print "Converged"
                break
        else:
            raise RuntimeError("SCF did not converge in '%s' iterations!" % max_iter)

        return Etot, orbitals, orbital_energies

    def solveKohnSham_new(self, max_iter=300, thresh=1.0e-4):
        """
        solve the Kohn-Sham equations iteratively on the multicenter grid

        Suppose we have an approximate solution phi with energy e
        that differs from the sought-after wavefuntion psi by dphi.

               H psi = E psi

        After making the ansatz

          psi = phi + dphi
          E   = e + de

        the Schroedinger equation turns into an inhomogeneous equation for the
        orbital correction dphi:

           (H - e - de) dphi = -(H-e) phi + de phi       (*)

        On the right hand side we define the residual R = (H-e) phi that tells
        us how much the solution phi deviates from fulfilling the Schroedinger eqn.

        In addition to the orbital correction, the inhomogeneous equation (*) 
        has the solution dphi = -phi, provided that phi is very close to the
        true solution:

              dphi = -phi  =>  (H - e)(-phi) + de phi = de phi - (H-e)(-phi)
                           =>  (H - e) phi = 0
                           =>  phi = psi

        The solution dphi = -phi is equivalent to the solution dphi = 0, which
        indicates convergence. Therefore we need to stop the orbital correction
        cycle if either
           (1) |R| = |(H-e) phi| = 0 
        or 
           (2) |dphi| = 0


        Optional
        ========
        max_iter    :   maximum number of SCF cycles, 
                        the solution of the Kohn-Sham equations and the self-consistency 
                        of the density are intertwined
        thresh      :   threshold for SCF convergence, the iteration is stopped if the
                        sum of the orbital corrections drops below this threshold.

        Returns
        =======
        Etot        :   total DFT energy
        orbitals    :   doubly occupied orbitals, list of callables, orbital[i](x,y,z)
        energies    :   list of Kohn-Sham orbital energies
        """
        print "initial orbital guess from DFTB calculation"
        orbitals = self.getOrbitalGuess()
        
        norb = len(orbitals)
        # all orbitals are doubly occupied
        nelec = 2*norb

        orbital_energies = np.zeros(norb)

        # attraction potential between nuclei and electrons
        nuclear_potential = nuclear_potential_func(self.atomlist)
        
        print "solve Kohn-Sham equations on a multicenter grid"

        print_grid_summary(self.atomlist,
                           settings.lebedev_order, settings.radial_grid_factor)

        for k in range(0, max_iter):
            print "Iteration k= %3.1d" % k
                       
            rho = density_func(orbitals)
            # effective Kohn-Sham potential without nuclear attraction
            # (only electron-electron interaction)
            veff_ee_new = effective_potential_func(self.atomlist, rho, self.xc,
                                                   nelec=nelec,
                                                   nuclear=False)
            if k > 0:
                if total_error < 1.0:
                    # The density is only updated if the orbitals have more
                    # or less stabilized.
                    print "updating effective potential veff[rho^(new)]"

                    #
                    #  (next)        (old)        (new)
                    # V       = 2/3 V      + 1/3 V
                    veff_ee = add_two_functions(self.atomlist,
                                                veff_ee, veff_ee_new, 2.0/3.0, 1.0/3.0)
            else:
                veff_ee = veff_ee_new

            # effective potential (electron-electron interaction + nuclei-electrons)
            def veff(x,y,z):
                return veff_ee(x,y,z) + nuclear_potential(x,y,z)
                
            total_error = 0.0

            orbital_corrections = []
            for i,phi in enumerate(orbitals):
                # compute corrections for each orbital
            
                # energy expectation value <phi|H|phi>
                e = energy_ee(self.atomlist, phi, phi, veff_ee)
                
                residual = residual_ee_func(self.atomlist, phi, veff_ee, e)
                delta_e = energy_correction(self.atomlist, residual, phi, method="Becke")
                source = source_func(delta_e, phi, residual)

                #
                nrm2_residual = overlap(self.atomlist, residual, residual)
                print "residual norm for orbital %d  |(H-E)phi|^2= %e" % (i+1, nrm2_residual)

                if nrm2_residual < thresh**2:
                    # If the residual norm is converged we don't have to
                    # compute the orbital correction for this orbital. 
                    print "Residual norm for orbital %d CONVERGED" % (i+1)
                    continue
                
                # orbital correction
                delta_phi = orbital_correction(self.atomlist, veff, source,  e, delta_e)
                nrm2_delta_phi = overlap(self.atomlist, delta_phi, delta_phi)
                total_error += nrm2_delta_phi
                
                print "    orbital i= %d    energy E= %e   delta E= %+e   |dphi|^2= %e" % (i+1, e, delta_e, nrm2_delta_phi)
                
                orbital_corrections.append(delta_phi)

                """
                ### DEBUG
                import matplotlib.pyplot as plt
                r = np.linspace(0.00001, 10.0, 5000)
                
                plt.plot(r, veff(r,0*r,0*r), ls="-.", lw=2, label=r"$v_{eff}$")

                plt.plot(r, residual(r,0*r,0*r), label=r"residual $(H-\epsilon)\phi$")
                plt.plot(r, source(r,0*r,0*r), label=r"source $\Delta \epsilon \phi - (H-E)\phi$")
                # wavefunctions
                plt.plot(r, phi(r,0*r,0*r), label=r"guess $\phi$")
                plt.plot(r, delta_phi(r,0*r,0*r), label=r"correction $\Delta \phi$")
                plt.plot(r, phi(r,0*r,0*r) + delta_phi(r,0*r,0*r), label=r"$\phi + \Delta \phi$")
                plt.plot(r, phi(r,0*r,0*r) - delta_phi(r,0*r,0*r), label=r"$\phi - \Delta \phi$")
                
                plt.legend()
                plt.show()
                ###
                """
                
            print "k= %3.1d  error = %e (threshold = %e)" % (k, total_error, thresh)

            # Since the Kohn-Sham equations are solved separately for
            # each orbital, the orthogonality may be lost, if one orbital
            # collapses to a lower-energy solution. 


            # Now we solve the Kohn-Sham equations in the basis of the current
            # orbitals and orbitals corrections {phi,delta phi} (at most 2*norb basis functions)
            # to restore orthogonality and improve the orbitals.
            orbital_energies, orbitals = solve_matrix_schroedinger_ee(self.atomlist,
                                                                      orbitals+orbital_corrections, veff_ee)
            # Only occupied orbitals are needed
            orbitals = orbitals[:norb]
            orbital_energies = orbital_energies[:norb]

            Etot = total_dft_energy(self.atomlist, orbitals, self.xc,
                                    nelec=nelec)
            print "total DFT energy Etot= %e" % Etot
            
            if total_error < thresh:
                print "SCF cycle CONVERGED"
                break
        else:
            raise RuntimeError("SCF did not converge in '%s' iterations!" % max_iter)

        return Etot, orbitals, orbital_energies

    
    def _writeIteration(self, k, atomlist, center,
                        E, charge, l, m, delta, 
                        phi0, residual0, phi, residual):
        """
        
        """
        print "compute radial wavefunctions and residuals for plotting..."
        # The radial grid is chosen for the atom with the largest
        # atomic number
        Zmax = max([Z for (Z,pos) in atomlist])
        # The origin of the radial grid is placed at `center`
        atom_center = (Zmax, center)

        # radial part of Coulomb wave without phase shift
        phi0_rad = radial_wave_func(atomlist, phi0, l, m,
                                    center=center)

        # spherical average of residual for phi0    R0 = (H-E)phi0
        residual0_avg = spherical_average_residual_func(atom_center, residual0)

        # radial part of normalized scattering solution
        phi_rad = radial_wave_func(atomlist, phi, l, m,
                                   center=center)

        # spherical average of residual for phi     R  = (H-E)phi
        residual_avg  = spherical_average_residual_func(atom_center, residual)
        
        # shifted regular coulomb function
        Cf_shift = regular_coulomb_func(E, charge, l, m, delta,
                                            center=center)
    
        # radial part of shifted Coulomb wave
        Cf_shift_rad = radial_wave_func(atomlist, Cf_shift, l, m,
                                        center=center)

        # save radial wavefunctions and spherically averaged residual
        print "# ITERATION k= %d RADIAL_WAVEFUNCTIONS" % k
        print "# Asymptotic wavefunction:"
        print "#   charge            Z= %+d" % charge
        print "#   energy            E= %e  k= %e" % (E, k)
        print "#   angular momentum  l= %d m= %+d" % (l, m)
        print "#   phase shift       delta= %e rad" % delta
        print "# "
        print "#     R/bohr                 Coulomb               Coulomb           radial wavefunction   spherical avg. residual"
        print "#                                                  shifted                R_{l,m}(r)           <|(H-E)phi|^2>"
        import sys
        # write table to console
        r = np.linspace(1.0e-3, 100, 1000)
        data = np.vstack((r, phi0_rad(r), Cf_shift_rad(r), phi_rad(r), residual_avg(r))).transpose()
        np.savetxt(sys.stdout, data, fmt="    %+e    ")
        print "# END"

        import matplotlib.pyplot as plt

        plt.title("Iteration k= %d" % k)
        # plot radial wavefunction phi, Coulomb wave phi0
        # and residuals
        r = np.linspace(1.0e-3, 100.0, 20000)
            
        plt.xlabel("r / bohr")

        # radial (l,m) component of wavefunctions
        plt.plot(r, phi0_rad(r), label=r"initial $\phi^{(0)}_{l=%d,m=%+d}(r)$" % (l,m))
        plt.plot(r, phi_rad(r), label=r"current $\phi_{l=%d,m=%+d}(r)$" % (l,m))
        
        # residuals (spherically averaged)
        plt.plot(r, residual0_avg(r), label=r"residual $\int d\Omega (H-E)\phi_0$")
        plt.plot(r, residual_avg(r), label=r"residual $\int d\Omega (H-E)\phi$")
            
        plt.legend()
        #plt.savefig("/tmp/iteration_%4.4d.png" % k)
        plt.show()
    
    def solveScatteringProblem(self, rho, homo, E, l, m,
                               nmax=5,
                               debug=1):
        """
        solve the Schroedinger equation (for 1 electron), the Hartree-Fock equation (for 2 electrons) 
        or the Kohn-Sham equations (for > 2 electrons) for a continuum orbital with
        energy E and asymptotic angular momentum (l,m) using Becke's multicenter grids

        Parameters
        ==========
        rho          :  callable rho(x,y,z), evaluates the electron density,
                        or None for systems with only 1 electron
        homo         :  callable homo(x,y,z), evaluates the wavefunction of the
                        HOMO, or None for systems with only 1 electrons
        E            :  energy of continuum orbital, E=1/2 k^2 > 0
        l,m          :  angular quantum numbers of asymptotic wavefunction
                        -l <= m <=l

        Optional
        ========
        nmax         :  maximum number of iterations employed in inverting
                        the multicenter spherical averaging operation S.V = U. 
                        nmax is the number of times the iteration

                           V_{n} = (1 - S) V_{n-1}      with V_0 = V

                        is performed. The desired inverse is obtained by summing
                        this series
                                  nmax
                           U = sum     v_{n}
                                  n=0

        debug        :  for debug > 0, radial wavefunctions and residuals are
                        plotted at the end

        Returns
        =======
        delta        :  phase shift (in radians)
        phi          :  callable phi(x,y,z), evaluates the continuum wavefunction
                        on a grid
        """
        # The continuum orbital is specified by its energy and asymptotic
        # angular momentum (E,l,m)
    
        # length of wave vectors
        k = np.sqrt(2*E)

        print " "
        print "  Asymptotic continuum wavefunction"
        print "  ================================="
        print "  energy           E= %e Hartree  ( %e eV )" % (E, E*AtomicData.hartree_to_eV)
        print "  wavevector       k= %e a.u." % k
        print "  angular moment   l= %d m= %+d" % (l,m)
        print " "

        atomlist = self.atomlist
        atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
        # show number of radial and angular points in multicenter grid
        print_grid_summary(atomlist,
                           settings.lebedev_order, settings.radial_grid_factor)
        
        print "total charges..."
        qnuc, qelec, qtot = total_charge(atomlist, rho)
        # number of electrons from integration
        nelec = int(np.round(-qelec))
    
        # We have to fix the center of the coordinate system relative
        # to which the asymptotic Coulomb waves are defined. The origin
        # R0 should be chosen such that the effective potential at large
        # distances appraches
        #                           charge
        #     V(r) --->  V0(r) = - -------
        #                           |r-R0|
        #        
        
        # In a neutral molecule the asymptotic behaviour of the effective
        # Kohn-Sham potential
        #
        #     V    (r)   ---> - 1/r               for r -> oo
        #      eff
        # is not the result of electrostatics (because the monopole of the
        # total charge distribution vanishes) but of the asymptotic behaviour
        # of the exchange functional
        #
        #     V        [rho](r)  --->  -1/r       for r -> oo
        #      exchange 
        
        #
        # The origin r0 such that
        #
        #     V_exchange ----> -1/|r-r0|
        #
        # is taken as the center of electronic charge in the HOMO
        #
        #    r0 = <HOMO|r|HOMO>
        # 
        #
        print ""
        print "  Origin"
        print "  ======"
        
        if nelec > 1:
            # density of the highest occupied molecular orbital
            def rho_HOMO(x,y,z):
                return homo(x,y,z)**2
        
            print "integrals over charge density of HOMO orbital..."
            center_nuc, center_homo, center_tot = center_of_charge(atomlist, rho_HOMO, verbose=1)

            print "  The center of the electron density in the HOMO"
            print "  is chosen as the the origin."
            center = center_homo
        else:
            center_nuc, dummy, dummy = center_of_charge(atomlist, None)

            print "  As there are no other electrons, the center of the nuclear charge"
            print "  is chosen as the origin."
            center = center_nuc

        print "  The origin of the coordinate system, relative to which"
        print "  asymptotic solutions are defined, is placed at"
        print ""
        print "    center = %+e  %+e  %+e" % tuple(center)
        print ""

        # Asymptotically the effective potential is Coulomb-like -(Q+1)/r . It's strength consists
        # of an electrostatic part (Q=0 for a neutral molecule, Q=+1 for a cation) and an
        # exchange part (+1, since v_xc[rho] -> -1/r for r -> oo)
        # if there are other electrons present.
        charge = qtot
        if nelec > 1:
            # add charge felt due to exchange
            charge = qtot + 1
        print " electrostatic charge  qtot = %+e" % qtot
        print " asymptotic charge          = %+e" % charge
    
        print "effective potential..."
        veff = effective_potential_func(atomlist, rho, self.xc, nelec=nelec)

        """
        #### TEST
        # We place an additional dummy atom at the center, which increases
        # the resolution of the grid but does not change the nuclear potential
        
        # list of real nuclei
        atomlist_nuclei = atomlist
        # central dummy atom
        dummy_atom = (int(charge), center)
        # list of real nuclei + dummy atoms
        atomlist = atomlist_nuclei + [dummy_atom]
        print " "
        print "  A dummy atom with charge Z=%+d is placed at the center." % (int(charge))
        # show size of grid used for dummy atom
        print_grid_summary([dummy_atom],
                           settings.lebedev_order, settings.radial_grid_factor)
        ####
        """

        # solve S.U = Veff
        norm_rem, u = solve_spherical_potential(atomlist, veff, nmax=nmax)
        
        # reconstruct V from U
        #   V^(rec) = S.U
        v_rec, dummy = spherical_remainder(atomlist, u)

        # compute reconstruction error |Veff - V^(rec)|

        # The resolution of the grid is temporarily reduced for
        # computing the reconstruction error. This is faster and
        # ensures that no points lie at the singularity of the
        # Coulomb potential
        rfac, Lmax = settings.radial_grid_factor, settings.lebedev_order
        settings.radial_grid_factor, settings.lebedev_order = 5, 21    
        error = np.sqrt( integral(atomlist, lambda x,y,z: (veff(x,y,z) - v_rec(x,y,z))**2 ) )
        # reset original grid resolution
        settings.radial_grid_factor, settings.lebedev_order = rfac, Lmax
        
        print "reconstruction error: |Veff-V^(rec))| = %e" % error

        # initial guess for the continuum wavefunction
        
        # asymptotically correct solution for V0 = -charge/r 
        phi0 = regular_coulomb_func(E, charge, l, m, 0.0,
                                    center=center)
        
        # origin of Coulomb waves and V0
        x0,y0,z0 = center
        def v0(x,y,z):
            # shift coordinates
            x,y,z = x-x0,y-y0,z-z0
            # distance from origin
            r = np.sqrt(x*x+y*y+z*z)
            return -charge/r
        
        # total effective potential is decomposed into
        # asymptotic potential and sort-range potential:
        #  Veff = V0 + V1
        def v1(x,y,z):
            return veff(x,y,z) - v0(x,y,z)

        # right-hand side of inhomogeneous Schroedinger equation
        def source(x,y,z):
            return -v1(x,y,z) * phi0(x,y,z)
        
        ##### Visualization ##########
        if debug > 1:
            # show cuts along one axis
            import matplotlib.pyplot as plt
            r = np.linspace(-5.0,5.0, 5000)

            # cut along z-axis
            x = 0*r
            y = 0*r
            z = r
            
            plt.xlabel("x / bohr")
            plt.ylabel("potential / Hartree")
            
            plt.plot(r, veff(x,y,z), color="black", ls="-", lw=2, label=r"original $V$")
            plt.plot(r, u(x,y,z), color="red", ls="-.", lw=2, label=r"inverted $U = \hat{S}^{-1} V$")
            plt.plot(r, v_rec(x,y,z), color="orange", ls="--", lw=2, label=r"reconstructed $V^{(rec)} = \hat{S} U$")
            plt.plot(r, veff(x,y,z)-v_rec(x,y,z), color="green", label=r"deviation $V - V^{(rec)}$")

            plt.legend()
            plt.show()
        ################################
            
        #
        #  solve  (T + U - E) dphi = - V1 phi0
        #
        # for orbital correction dphi
        
        print "Schroedinger equation..."
        dphi = multicenter_inhomogeneous_schroedinger(u, source,  E,
                                                      atomic_coordinates, atomic_numbers,
                                                      radial_grid_factor=settings.radial_grid_factor,
                                                      lebedev_order=settings.lebedev_order)

        # Combine asymptotically correct solution with correction 
        #   phi = phi0 + dphi
        phi = add_two_functions(atomlist, phi0, dphi, 1.0, 1.0)

        print "normalization and phase shifts..."
        delta, phi = phaseshift(atomlist, phi, E, charge, l, m,
                                center=center)
        
        # The continuum orbital should be orthogonal to the bound
        # orbitals belonging to the same Hamiltonian. I think this
        # should come out correctly by default.
        if nelec > 1:
            olap_hc = overlap(atomlist, homo, phi)
            print "  overlap  <homo | continuum> = %e" % olap_hc

        ##### Visualization ##########
        if debug > 0:
            print "residuals..."
            # residual for phi0    R0 = (H-E)phi0
            residual0 = residual_func(atomlist, phi0, veff, E)
            # residual for final solution  R = (H-E)phi
            residual  = residual_func(atomlist, phi, veff, E)
            # spherical average of residual function
            # The origin is placed at a central dummy atom.
            # central dummy atom
            dummy_atom = (int(charge), center)
            residual_avg = spherical_average_residual_func(dummy_atom, residual)
                
            ##### table with radial wavefunctions and residual 
            print "radial wavefunctions..."
            # shifted regular coulomb function
            Cf_shift = regular_coulomb_func(E, charge, l, m, delta,
                                            center=center)
    
            # save radial wavefunctions and spherically averaged residual
            # radial part of Coulomb wave without phase shift
            phi0_rad = radial_wave_func(atomlist, phi0, l, m,
                                        center=center)
            # radial part of shifted Coulomb wave
            Cf_shift_rad = radial_wave_func(atomlist, Cf_shift, l, m,
                                            center=center)
            # radial part of normalized scattering solution
            phi_rad = radial_wave_func(atomlist, phi, l, m,
                                       center=center)

            print "# RADIAL_WAVEFUNCTIONS"
            print "# Asymptotic wavefunction:"
            print "#   charge            Z= %+d" % charge
            print "#   energy            E= %e  k= %e" % (E, k)
            print "#   angular momentum  l= %d m= %+d" % (l, m)
            print "#   phase shift       delta= %e rad" % delta
            print "# "
            print "#     R/bohr                 Coulomb               Coulomb           radial wavefunction   spherical avg. residual"
            print "#                                                  shifted                R_{l,m}(r)           <|(H-E)phi|^2>"
            import sys
            # write table to console
            r = np.linspace(1.0e-3, 100.0, 1000)
            data = np.vstack((r, phi0_rad(r), Cf_shift_rad(r), phi_rad(r), residual_avg(r))).transpose()
            np.savetxt(sys.stdout, data, fmt="    %+e    ")
            print "# END"

            ##### plot potentials, wavefunctions and residual and an axis 
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(3,1, sharex=True)
            r = np.linspace(-30.0, +30.0, 50000)

            # cut along x-axis
            x = 0*r
            y = 0*r
            z = r            
            
            axes[0].set_ylabel("potential / Hartree")
            axes[1].set_ylabel("wavefunction")
            axes[2].set_ylabel("residual")

            # potentials
            axes[0].plot(r, veff(x,y,z) , ls="-",  label=r"original $V$")
            axes[0].plot(r, u(x,y,z),     ls="-.", label=r"inverted $U = S^{-1} V$")
            axes[0].plot(r, v_rec(x,y,z), ls="--", label=r"reconstructed $V^{(rec)}$")
            axes[0].plot(r, source(x,y,z), ls="-.", label=r"source $-(V-V_0) \phi_0$")

            # wavefunctions
            axes[1].plot(r, phi0(x,y,z), ls="-.", label=r"$\phi_0$")
            axes[1].plot(r, phi(x,y,z),  ls="-" , label=r"$\phi$")

            ### LCAO continuum wavefunction
            for i,(Z,pos) in enumerate(atomlist):
                # Hydrogen continuum orbital on center i
                phiI = regular_coulomb_func(E, +1, l, m, 0.0,
                                            center=pos)  
                if i == 0:
                    phi_LCAO = phiI
                else:
                    phi_LCAO = add_two_functions(atomlist, phi_LCAO, phiI, 1.0, 1.0)
            axes[1].plot(r, phi_LCAO(x,y,z), ls="--", label=r"$\phi_{LCAO}$")
            
            ##################################
            
            # residuals
            axes[2].plot(r, residual0(x,y,z), ls="-.", label=r"residual $(T+V-E)\phi_0$")
            axes[2].plot(r, residual(x,y,z),  ls="-",  label=r"residual $(T+V-E)\phi$")

            ######################################
            # residual of LCAO wavefunction
            residual_LCAO  = residual_func(atomlist, phi_LCAO, veff, E)                    
            axes[2].plot(r, residual_LCAO(x,y,z), ls="--", label=r"residual $(T+V-E)\phi_{LCAO}$")
            ######################################
            
            for ax in axes:
                ax.set_xlabel("x / bohr   (y=0, z=0)")
                ax.legend()
        
            fig.tight_layout()
            
            plt.show()

            
        ############################################
        
        # Here the k-dependent factors of the normalization are included,
        # so that the normalized continuum wave function has the asymptotic form
        #
        #  phi_norm_k(r) --->  k^{-1/2} 1/r sin(k*r + ...)
        # 
        phi_norm_k = multicenter_operation([phi], lambda fs: fs[0]/np.sqrt(k), 
                                         atomic_coordinates, atomic_numbers,
                                         radial_grid_factor=settings.radial_grid_factor,
                                         lebedev_order=settings.lebedev_order)
            
        return delta, phi_norm_k

    def solveScatteringProblem_new(self, rho, E, l,m, 
                                   nmax=5, 
                                   debug=1):
        """
        NOT WORKING BECAUSE OF CONCEPTUAL MISTAKE

        solve the Schroedinger equation (for 1 electron), the Hartree-Fock equation (for 2 electrons) 
        or the Kohn-Sham equations (for > 2 electrons) for a continuum orbital with
        energy E and asymptotic angular momentum (l,m) using Becke's multicenter grids

        Parameters
        ==========
        rho          :  callable rho(x,y,z), evaluates the electron density,
                        or None for systems with only 1 electron
        E            :  energy of continuum orbital, E=1/2 k^2 > 0
        l,m          :  angular quantum numbers of the asymptotic wavefunction of interest
                        -l <= m <=l


        Optional
        ========
        nmax         :  maximum number of iterations employed in inverting
                        the multicenter spherical averaging operation S.V = U. 
                        nmax is the number of times the iteration

                           V_{n} = (1 - S) V_{n-1}      with V_0 = V

                        is performed. The desired inverse is obtained by summing
                        this series
                                  nmax
                           U = sum     v_{n}
                                  n=0

        debug        :  for debug > 0, radial wavefunctions and residuals are
                        plotted at the end

        Returns
        =======
        delta        :  phase shift (in radians)
        phi          :  callable phi(x,y,z), evaluates the continuum wavefunction
                        on a grid, the asymptotic angular momentum can be specified
                        via keywords, e.g. phi(x,y,z,l=1,m=0). If no keywords are given,
                        the defaults are l=0, m=0.
        """
        # The continuum orbital is specified by its energy and asymptotic
        # angular momentum (E,l,m)
        # 
    
        # length of wave vectors
        k = np.sqrt(2*E)

        print " "
        print "  Asymptotic continuum wavefunction"
        print "  ================================="
        print "  energy           E= %e Hartree  ( %e eV )" % (E, E*AtomicData.hartree_to_eV)
        print "  wavevector       k= %e a.u." % k
        print " "

        atomlist = self.atomlist
        atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
        # show number of radial and angular points in multicenter grid
        print_grid_summary(atomlist,
                           settings.lebedev_order, settings.radial_grid_factor)
        
        print "total charges..."
        qnuc, qelec, qtot = total_charge(atomlist, rho)
        # number of electrons from integration
        nelec = int(np.round(-qelec))
    
        # We have to fix the center of the coordinate system relative
        # to which the asymptotic Coulomb waves are defined. The origin
        # R0 should be chosen such that the effective potential at large
        # distances appraches
        #                           charge
        #     V(r) --->  V0(r) = - -------
        #                           |r-R0|
        #        
        
        # In a neutral molecule the asymptotic behaviour of the effective
        # Kohn-Sham potential
        #
        #     V    (r)   ---> - 1/r               for r -> oo
        #      eff
        # is not the result of electrostatics (because the monopole of the
        # total charge distribution vanishes) but of the asymptotic behaviour
        # of the exchange functional
        #
        #     V        [rho](r)  --->  -1/r       for r -> oo
        #      exchange 
        
        #
        # The origin r0 such that
        #
        #     V_exchange ----> -1/|r-r0|
        #
        # is taken as the center of electronic charge in the HOMO
        #
        #    r0 = <HOMO|r|HOMO>
        # 
        #
        print ""
        print "  Origin"
        print "  ======"
        
        if nelec > 1:
            # density of the highest occupied molecular orbital
            def rho_HOMO(x,y,z):
                return homo(x,y,z)**2
        
            print "integrals over charge density of HOMO orbital..."
            center_nuc, center_homo, center_tot = center_of_charge(atomlist, rho_HOMO, verbose=1)

            print "  The center of the electron density in the HOMO"
            print "  is chosen as the the origin."
            center = center_homo
        else:
            center_nuc, dummy, dummy = center_of_charge(atomlist, None)

            print "  As there are no other electrons, the center of the nuclear charge"
            print "  is chosen as the origin."
            center = center_nuc

        print "  The origin of the coordinate system relative to which"
        print "  asymptotic solutions are defined is placed at"
        print ""
        print "    center = %+e  %+e  %+e" % tuple(center)
        print ""

        # Asymptotically the effective potential is Coulomb-like -(Q+1)/r . It's strength consists
        # of an electrostatic part (Q=0 for a neutral molecule, Q=+1 for a cation) and an
        # exchange part (+1, since v_xc[rho] -> -1/r for r -> oo)
        # if there are other electrons present.
        charge = qtot
        if nelec > 1:
            # add charge felt due to exchange
            charge = qtot + 1
        print " electrostatic charge  qtot = %+e" % qtot
        print " asymptotic charge          = %+e" % charge
    
        print "effective potential..."
        veff = effective_potential_func(atomlist, rho, self.xc, nelec=nelec)

        # solve S.U = Veff
        norm_rem, u = solve_spherical_potential(atomlist, veff, nmax=nmax)
        
        # reconstruct V from U
        #   V^(rec) = S.U
        v_rec, dummy = spherical_remainder(atomlist, u)

        # compute reconstruction error |Veff - V^(rec)|

        # The resolution of the grid is temporarily reduced for
        # computing the reconstruction error. This is faster and
        # ensures that no points lie at the singularity of the
        # Coulomb potential
        rfac, Lmax = settings.radial_grid_factor, settings.lebedev_order
        settings.radial_grid_factor, settings.lebedev_order = 5, 21    
        error = np.sqrt( integral(atomlist, lambda x,y,z: (veff(x,y,z) - v_rec(x,y,z))**2 ) )
        # reset original grid resolution
        settings.radial_grid_factor, settings.lebedev_order = rfac, Lmax
        
        print "reconstruction error: |Veff-V^(rec))| = %e" % error
            
        #
        #  solve  (T + U - E) phi = 0
        #
        # for for continuum orbital phi
        
        print "continuum Schroedinger equation..."
        phi = multicenter_continuum_schroedinger(u, E, charge, 
                                                 atomic_coordinates, atomic_numbers,
                                                 radial_grid_factor=settings.radial_grid_factor,
                                                 lebedev_order=settings.lebedev_order)

        #####
        #print "normalization and phase shifts..."
        #delta, phi = phaseshift(atomlist, lambda x,y,z: phi(x,y,z,l=l,m=m),
        #                        E, charge, l, m,
        #                        center=center)
        delta = 0.0
        ####
        
        ##### Visualization ##########
        if debug > 0:
            print "residuals..."
            # regular coulomb function without phase shift
            phi0 = regular_coulomb_func(E, charge, l, m, 0.0,
                                        center=center)
            # residual for phi0    R0 = (H-E)phi0
            residual0 = residual_func(atomlist, phi0, veff, E)
            # residual for final solution  R = (H-E)phi
            residual  = residual_func(atomlist, phi, veff, E)
            # Laplacian
            lap_phi = laplacian_func(atomlist, phi)
            # spherical average of residual function
            # The origin is placed at a central dummy atom.
            # central dummy atom
            dummy_atom = (int(charge), center)
            residual_avg = spherical_average_residual_func(dummy_atom, residual)
                
            ##### table with radial wavefunctions and residual 
            print "radial wavefunctions..."
            # shifted regular coulomb function
            Cf_shift = regular_coulomb_func(E, charge, l, m, delta,
                                            center=center)
    
            # save radial wavefunctions and spherically averaged residual
            # radial part of Coulomb wave without phase shift
            phi0_rad = radial_wave_func(atomlist, phi0, l, m,
                                        center=center)
            # radial part of shifted Coulomb wave
            Cf_shift_rad = radial_wave_func(atomlist, Cf_shift, l, m,
                                            center=center)
            # radial part of normalized scattering solution
            phi_rad = radial_wave_func(atomlist, phi, l, m,
                                       center=center)

            print "# RADIAL_WAVEFUNCTIONS"
            print "# Asymptotic wavefunction:"
            print "#   charge            Z= %+d" % charge
            print "#   energy            E= %e  k= %e" % (E, k)
            print "#   angular momentum  l= %d m= %+d" % (l, m)
            print "#   phase shift       delta= %e rad" % delta
            print "# "
            print "#     R/bohr                 Coulomb               Coulomb           radial wavefunction   spherical avg. residual"
            print "#                                                  shifted                R_{l,m}(r)           <|(H-E)phi|^2>"
            import sys
            # write table to console
            r = np.linspace(1.0e-3, 100.0, 1000)
            data = np.vstack((r, phi0_rad(r), Cf_shift_rad(r), phi_rad(r), residual_avg(r))).transpose()
            np.savetxt(sys.stdout, data, fmt="    %+e    ")
            print "# END"

            ##### plot potentials, wavefunctions and residual and an axis 
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(3,1, sharex=True)
            r = np.linspace(-20.0, +20.0, 50000)

            # cut along x-axis
            x = 0*r
            y = 0*r
            z = r            
            
            axes[0].set_ylabel("potential / Hartree")
            axes[1].set_ylabel("wavefunction")
            axes[2].set_ylabel("residual")

            # potentials
            axes[0].plot(r, veff(x,y,z) , ls="-",  label=r"original $V$")
            axes[0].plot(r, u(x,y,z),     ls="-.", label=r"inverted $U = S^{-1} V$")
            axes[0].plot(r, v_rec(x,y,z), ls="--", label=r"reconstructed $V^{(rec)}$")

            # wavefunctions
            axes[1].plot(r, phi0(x,y,z), ls="-.", label=r"$\phi_0$")
            axes[1].plot(r, phi(x,y,z),  ls="-" , label=r"$\phi$")

            ### LCAO continuum wavefunction
            for i,(Z,pos) in enumerate(atomlist):
                # Hydrogen continuum orbital on center i
                phiI = regular_coulomb_func(E, +1, l, m, 0.0,
                                            center=pos)  
                if i == 0:
                    phi_LCAO = phiI
                else:
                    phi_LCAO = add_two_functions(atomlist, phi_LCAO, phiI, 1.0, 1.0)
            axes[1].plot(r, phi_LCAO(x,y,z), ls="--", label=r"$\phi_{LCAO}$")
            
            ##################################
            
            # residuals
            axes[2].plot(r, residual0(x,y,z), ls="-.", label=r"residual $(T+V-E)\phi_0$")
            axes[2].plot(r, residual(x,y,z),  ls="-",  label=r"residual $(T+V-E)\phi$")

            # kinetic vs. potential energy
            axes[2].plot(r, -0.5 * lap_phi(x,y,z), ls="--", label=r"kinetic energy $T\phi$")
            axes[2].plot(r, (veff(x,y,z)-E)*phi(x,y,z), ls="--", label=r"potential energy $(V-E)\phi$")
            
            ######################################
            # residual of LCAO wavefunction
            residual_LCAO  = residual_func(atomlist, phi_LCAO, veff, E)                    
            axes[2].plot(r, residual_LCAO(x,y,z), ls="--", label=r"residual $(T+V-E)\phi_{LCAO}$")
            ######################################
            
            for ax in axes:
                ax.set_xlabel("x / bohr   (y=0, z=0)")
                ax.legend()
        
            fig.tight_layout()
            
            plt.show()

            
        ############################################
        
        # Here the k-dependent factors of the normalization are included,
        # so that the normalized continuum wave function has the asymptotic form
        #
        #  phi_norm_k(r) --->  k^{-1/2} 1/r sin(k*r + ...)
        # 
        phi_norm_k = multicenter_operation([phi], lambda fs: fs[0]/np.sqrt(k), 
                                         atomic_coordinates, atomic_numbers,
                                         radial_grid_factor=settings.radial_grid_factor,
                                         lebedev_order=settings.lebedev_order)
            
        return delta, phi_norm_k
        

    
###############################################################################################
#
# Scattering, phase shifts
#
###############################################################################################

import mpmath

def regular_coulomb_func(E, Z, l, m, delta_l,
                         center=(0.0, 0.0, 0.0)):
    """
    create a regular Coulomb function, which is a continuum solution of the 
    Schroedinger equation with energy E for the hydrogen atom (if Z=1)

              __2
        (-1/2 \/  -  Z/r) Cf(x,y,z) = E Cf(x,y,z)


    The angular part of Cf are real (!) spherical harmonics.

    Parameters
    ==========
    E          : energy of solution, > 0
    Z          : charge of ion, > 0
    l,m        : angular momentum quantum numbers, -l <= m <= l
    delta_l    : phase shift

    Optional
    ========
    center     : tuple (x0,y0,z0) with coordinates of the origin
                 to which the Coulomb function should be shifted

    Returns
    =======
    Cf         : function, Cf(x,y,z) evaluates the Coulomb function
                 on a grid
    """
    # origin
    x0,y0,z0 = center
    
    def Cf(x,y,z):
        # shift coordinates to new origin
        x,y,z = x-x0, y-y0, z-z0
        # spherical coordinates
        r,th,ph = cartesian2spherical((x,y,z))
        
        k = np.sqrt(2*E)
        kr = k*r

        # radial part R_{E,l}(r) of wavefunction
        rhos = kr.flatten()
        fs = []
        for rho in rhos:
            fi = mpmath.coulombf(l, -Z/k, rho+delta_l)
            fs.append(complex(fi).real)
        wfn_radial = np.reshape(np.array(fs), kr.shape)/r 
        
        # angular part of wavefunction, Y^(real)_{l,m}(th,ph)
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,ll,mm in sph_it:
            #if ll == l and mm == m:
            #    wfn_angular = Ylm
            #    break
            if ll == l and mm == m:
                # real spherical harmonics
                if m < 0:
                    Ylm_real = -np.sqrt(2.0) * Ylm.imag
                elif m > 0:
                    Ylm_real =  np.sqrt(2.0) * (-1)**m * Ylm.real
                else:
                    # m == 0
                    Ylm_real = Ylm.real
                    
                wfn_angular = Ylm_real
                break

        # combine radial and angular parts
        wfn = wfn_radial * wfn_angular
        
        return wfn

    return Cf

def irregular_coulomb_func(E, Z, l, m, delta_l,
                         center=(0.0, 0.0, 0.0)):
    """
    create a irregular Coulomb function, which is a continuum solution of the 
    Schroedinger equation with energy E for the hydrogen atom (if Z=1)

              __2
        (-1/2 \/  -  Z/r) Cg(x,y,z) = E Cf(x,y,z)


    The angular part of Cg are real (!) spherical harmonics.

    Parameters
    ==========
    E          : energy of solution, > 0
    Z          : charge of ion, > 0
    l,m        : angular momentum quantum numbers, -l <= m <= l
    delta_l    : phase shift

    Optional
    ========
    center     : tuple (x0,y0,z0) with coordinates of the origin
                 to which the Coulomb function should be shifted

    Returns
    =======
    Cg         : function, Cg(x,y,z) evaluates the Coulomb function
                 on a grid
    """
    # origin
    x0,y0,z0 = center
    
    def Cg(x,y,z):
        # shift coordinates to new origin
        x,y,z = x-x0, y-y0, z-z0
        # spherical coordinates
        r,th,ph = cartesian2spherical((x,y,z))
        
        k = np.sqrt(2*E)
        kr = k*r

        # radial part R_{E,l}(r) of wavefunction
        rhos = kr.flatten()
        fs = []
        for rho in rhos:
            fi = mpmath.coulombg(l, -Z/k, rho+delta_l)
            fs.append(complex(fi).real)
        wfn_radial = np.reshape(np.array(fs), kr.shape)/r 
        
        # angular part of wavefunction, Y^(real)_{l,m}(th,ph)
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,ll,mm in sph_it:
            #if ll == l and mm == m:
            #    wfn_angular = Ylm
            #    break
            if ll == l and mm == m:
                # real spherical harmonics
                if m < 0:
                    Ylm_real = -np.sqrt(2.0) * Ylm.imag
                elif m > 0:
                    Ylm_real =  np.sqrt(2.0) * (-1)**m * Ylm.real
                else:
                    # m == 0
                    Ylm_real = Ylm.real
                    
                wfn_angular = Ylm_real
                break

        # combine radial and angular parts
        wfn = wfn_radial * wfn_angular
        
        return wfn

    return Cg


def radial_wave_func(atomlist, wfn, l, m, center=(0.0, 0.0, 0.0)):
    """
    compute the radial wavefunction R_{l,m}(r) by projecting the wavefunction
    wfn(x,y,z) onto the real spherical harmonic Y_{l,m}(th,ph)
                /
      R  (r)  = | Y   (Omega) wfn(r,Omega) dOmega
       l,m      /  l,m

    Parameters
    ==========
    atomlist   :  list of tuples (Zat,[xi,yi,zi]) with atomic numbers
                  and positions
    wfn        :  callable, wfn(x,y,z) evaluates the wavefunction on a grid
    l,m        :  integers, -l <= m <= l, angular quantum numbers of the radial 
                  wavefunction that should be projected out

    Optional
    ========
    center     : tuple (x0,y0,z0) with coordinates of the origin
                 around which the spherical wave expansion should 
                 be performed

    Returns
    =======
    wfn_rad    :  callable, wfn_rad(r) evaluates the radial wavefunction R_{l,m}(r)
    """
    # The radial grid is chosen for the atom with the largest atomic number
    Zmax = max([Z for (Z,pos) in atomlist])
    # The origin of the radial grid is placed at `center`
    atom_center = (Zmax, center)

    wfn_rad = radial_component_func(atom_center, wfn, l,m,
                                    radial_grid_factor=settings.radial_grid_factor,
                                    lebedev_order=settings.lebedev_order)                                    

    return wfn_rad

def phaseshift_lstsq(atomlist, phi, E, Z, l, m, rmin, rmax, Npts,
                     center=(0.0, 0.0, 0.0),
                     debug=0):
    """
    determine the scale factor and phaseshift of a continuum orbital `phi` relative
    to the regular Coulomb wavefunction

    Parameters
    ==========
    atomlist   :  list of tuples (Z,[x,y,z]) with atomic numbers
                  and positions
    phi        :  continuum orbital with unknown phase shift
    E          :  energy of continuum state
    Z          :  nuclear charge, Z=+1 for hydrogen
    l,m        :  angular quantum numbers, -l <= m <= l
    rmin,rmax  :  upper and lower bound for selecting points that
                  contribute to the matching error
    Npts       :  number of samples

    Optional
    ========
    center     :  tuple (x0,y0,z0) with position of the center of charge
    debug      :  if > 0, plot f(delta) and its root delta0

    Returns
    =======
    scale      :  scale factor, 1/scale * phi(x,y,z) has asymptotically the same
                  amplitude as the Coulomb wavefunction
    delta      :  phase shift in radians
    """
    # select radial points
    r = np.linspace(rmin, rmax, Npts)
    # weights
    w = np.ones(Npts)/float(Npts)
    # evaluate radial wavefunction at sampling points
    wfn_rad = radial_wave_func(atomlist, phi, l, m, center=center)
    # R = 1/r u
    # values of radial wavefunction at the sampling points
    u = r * wfn_rad(r)
    
    # length of wave vectors
    k = np.sqrt(2.0*E)

    # Coulomb function with phase shift delta
    def coulf(r, delta):
        cf = np.zeros(r.shape)
        for i,ri in enumerate(r):
            cfi = mpmath.coulombf(l, -Z/k, k*ri + delta)
            cf[i] = float(complex(cfi).real)
        return cf

    # The objective function to minimize is the
    # sum of squared deviations between the radial wavefunction and
    # the shifted Coulomb function
    #        2                                        *                 2
    #   error (delta) = sum  weight(i) ( r(i) R(i) - s  Cf(r(i),delta) )
    #                     i
    #
    # s*(delta) is the optimal scale factor for a given phase shift delta.
    #
    def error2(delta):
        # values of shifted Coulomb function at sampling points
        c = coulf(r, delta)
        # optimal scale factor for given delta
        scale = np.sum(w * u * c)/np.sum(w * c**2)
        # sum of squared errors
        err2 = np.sum(w * (u - scale * c)**2)
        return err2
    
    # find phase shift delta the minimizes the sum of squared errors
    delta_opt = optimize.fminbound(error2, 0.0, 2.0*np.pi)
    # find scale factor for optimal phase shift
    c_opt = coulf(r, delta_opt)
    scale_opt = np.sum(w * u * c_opt)/np.sum(w * c_opt**2)    

    # 
    err2 = np.sum(w * (u - scale_opt * coulf(r, delta_opt))**2)
    err = np.sqrt(err2)
    print "average deviation after matching to Coulomb function err= %e" % err

    if debug > 0:
        # Where does the minimum of error2(delta) lie?
        # Have we found the correct minimum
        import matplotlib.pyplot as plt
        ds = np.linspace(0.0, 2*np.pi, 1000)
        plt.xlabel(r"phase shift $\delta$")
        plt.ylabel(r"matching error$^2$")
        error2_array = [error2(d) for d in ds]
        plt.plot(ds, error2_array, label=r"$\epsilon(\delta)^2 = \vert r R(r) - s Cf(r,\delta) \vert^2$")
        plt.plot([delta_opt], [err2], "o", label=r"$\epsilon(\delta^*)$")
        plt.legend()
        plt.show()

    
    return delta_opt, scale_opt

def phaseshift_2pts(atomlist, phi, E, Z, l, m, r1, r2,
                    center=(0.0, 0.0, 0.0),
                    debug=0):
    """
    determine the scale factor and phaseshift of a continuum orbital `phi` relative
    to the regular Coulomb wavefunction from matching at two points

    Parameters
    ==========
    atomlist   :  list of tuples (Z,[x,y,z]) with atomic numbers
                  and positions
    phi        :  continuum orbital with unknown phase shift
    E          :  energy of continuum state
    Z          :  nuclear charge, Z=+1 for hydrogen
    l,m        :  angular quantum numbers, -l <= m <= l
    r1, r2     :  points at which the radial functions of the continuum
                  orbital and the Coulomb wave are matched

    Optional
    ========
    center     :  tuple (x0,y0,z0) with position of the center of charge
    debug      :  if > 0, plot f(delta) and its root delta0

    Returns
    =======
    scale      :  scale factor, 1/scale * phi(x,y,z) has asymptotically the same
                  amplitude as the Coulomb wavefunction
    delta      :  phase shift in radians
    """
    print "continuum and Coulomb wave are matched at the points"
    print "  r1= %e  bohr" % r1
    print "  r2= %e  bohr" % r2
    # integrate out the angular part Y_{l,m} of the wavefunction
    wfn_rad = radial_wave_func(atomlist, phi, l, m, center=center)
    # evaluate radial wavefunction at the two sampling points
    #  R = 1/r u
    u1 = r1 * wfn_rad(r1)
    u2 = r2 * wfn_rad(r2)

    # length of wave vectors
    k = np.sqrt(2.0*E)

    # Coulomb function with phase shift delta
    def coulf(r, delta):
        cf = mpmath.coulombf(l, -Z/k, k*r + delta)
        cf = float(complex(cf).real)
        return cf

    # The wavefunctions should agree at two points
    #
    #    (1)   u(r1) = s * c(r1;delta)
    #    (2)   u(r2) = s * c(r2;delta)
    #
    # We can solve the first equation for s = u(r1)/c(r1;delta),
    # and substitute this into the second equation. The phase shift
    # delta is then a root of
    #
    #    f(delta) = u(r2) * c(r1;delta) - u(r1) * c(r2;delta)
    #
    # The scale factor is determined for the root delta0 either from
    # eqn. (1) or (2) as
    #
    #    scale = u(r1)/c(r1;delta_root)
    # or
    #    scale = u(r2)/c(f2;delta_root)

    def f(delta):
        return u2 * coulf(r1,delta) - u1 * coulf(r2,delta)

    # find root of f(delta) with Newton's method
    delta_guess = 0.0
    
    try:
        delta_root = optimize.newton(f, delta_guess)
    except RuntimeError as e:
        msg  = "\nHINT:  You probably have to increase the resolution of\n"
        msg += "       the multicenter grid. If the computed radial wavefunction\n"
        msg += "       is not accurate enough at large radii it can happen that f(delta)\n"
        msg += "       has no roots at all!\n\n"
        print msg
        raise e

    # compute scale factor
    c1 = coulf(r1,delta_root)
    c2 = coulf(r1,delta_root)
    if (abs(c1) > 0.0):
        scale = u1/c1
    else:
        scale = u2/c2
        
    print "Optimized phase shift  delta= %e" % delta_root
    print "  f(delta)= %e" % f(delta_root)
    print "scale factor           scale= %e" % scale

    if debug > 0:
        # Where are the roots of f(delta)?
        import matplotlib.pyplot as plt
        ds = np.linspace(0.0, 2*np.pi, 1000)
        plt.xlabel(r"phase shift $\delta$")
        plt.ylabel(r"$f(\delta)$")
        fs = [f(d) for d in ds]
        # 
        plt.plot(ds, fs, label=r"$f(\delta) = u(r_2) c(r_1;\delta) - u(r_1) c(r_2;\delta)$")
        # highlight position of optimized phase shift
        plt.plot([delta_root], [f(delta_root)], "o", label=r"$f(\delta_0)$")
        # x-axis
        plt.plot([0.0,2.0*np.pi],[0.0,0.0], color="red")
        
        plt.legend()
        plt.show()
        
    return scale, delta_root

    
def phaseshift2overlaps(E, Z, l,m, delta_l, rmin, rmax):
    """
    compute the overlaps of two Coulomb functions over the interval [rmin,rmax],
    where one or both Coulomb functions have a phase shift delta_l.

    It is assumed that rmin is large enough so that the asymptotic
    expansion of the regular Coulomb function 

         F(rho) = sin(rho - eta log(2*rho) - l pi/2 + sigma_l)

    is warranted. For this expansion the integrals
                         /rho_max
      I_0d(delta) =  1/k |        F(rho) F(rho+delta) drho
                         /rho_min
    and
                         /rho_max
      I_dd(delta) =  1/k |        F(rho+delta) F(rho+delta) drho
                         /rho_min
    
    can be performed analytically (with Mathematica).

    Parameters
    ==========
    E          :  energy of continuum state
    Z          :  nuclear charge, Z=+1 for hydrogen
    l,m        :  angular quantum numbers, -l <= m <= l
    delta_l    :  phase shift in radians
    rmin,rmax  :  upper and lower bound for integration over
                  radial coordinate

    Returns
    =======
    I_0d       :  float, overlap of unshifted and shifted Coulomb function
                  on the interval [rmin,rmax]
                                /rmax
                  I_0d(delta) = |     dr   F(kr) F(kr+delta)
                                /rmin
    I_dd       :  float, overlap of shifted Coulomb function with itself
                  on the interval [rmin,rmax]
                                /rmax
                  I_dd(delta) = |     dr   F(kr+delta) F(kr+delta)
                                /rmin
    """
    k = np.sqrt(2*E)
    eta = -Z/k
    # Coulomb phase shift, sigma_l = arg Gamma(l+1+i*eta)
    g = mpmath.gammainc(l+1.0+1.0j*eta)
    sigma_l = np.angle( complex(g) )

    def indefinite_integral_0d(rho):
        #             /rho
        # I_0d(rho) = | F(rho') F(rho'+delta) drho'
        #             /0
        a = 1.0j/8.0 * (-1.0)**l * np.exp(-1.0j*delta_l + np.pi*eta - 2.0j*sigma_l)
        g1 = np.exp(2.0j*(delta_l + 2.0*sigma_l)) * complex(mpmath.gammainc(1.0-2.0j*eta, -2.0j*rho))
        g2 = complex(mpmath.gammainc(1.0+2.0j*eta, 2.0j*rho))
        
        olap = 0.5*rho*np.cos(delta_l) + a * (g1-g2)
        
        return olap.real

    def indefinite_integral_dd(rho):
        #             /rho
        # I_dd(rho) = | F(rho'+delta) F(rho'+delta) drho'
        #             /0
        a = 1.0j/8.0 * (-1.0)**l * np.exp(-2.0j*(delta_l + sigma_l) + np.pi*eta)
        g1 = np.exp(4.0j*(delta_l + sigma_l)) * complex(mpmath.gammainc(1.0-2.0j*eta, -2.0j*rho))
        g2 = complex(mpmath.gammainc(1.0+2.0j*eta, 2.0j*rho))
        
        olap = 0.5*rho + a * (g1-g2)
        
        return olap.real

    # definite overlap integral over the interval [rmin,rmax]
    # after substituting rho = k*r
    #
    #  /rmax                                     /k*rmax
    #  |      F(rho(r)) F(rho(r)+delta) dr = 1/k |       F(rho) F(rho+delta) drho
    #  /rmin                                     /k*rmin
    #
    #                                      = 1/k [I_0d(k*rmax) - I_0d(k*rmax)]
    I_0d = 1.0/k * (indefinite_integral_0d(k*rmax) - indefinite_integral_0d(k*rmin))
    # similary for I_dd
    I_dd = 1.0/k * (indefinite_integral_dd(k*rmax) - indefinite_integral_dd(k*rmin))

    return I_0d, I_dd

def overlaps2phaseshift(E, Z, l, m, olap_0d, olap_dd, rmin, rmax):
    """
    determine scale factor and phaseshift of continuum orbital relative to the 
    Coulomb wavefunction from the overlaps on a radial interval [rmin,rmax].

    Asymptotically the radial part of the continuum orbital behaves as

      1/r scale * C(r+delta)

    The parameters `scale` and `delta` are determined such that the following
    two equations are satisfied:

      scale * I_0d(delta)   = olap_0d
      scale^2 * I_dd(delta) = olap_dd

    (see explanations in `phaseshift2overlaps()` above)

    Parameters
    ==========
    E          :  energy of continuum state
    Z          :  nuclear charge, Z=+1 for hydrogen
    l,m        :  angular quantum numbers, -l <= m <= l
    olap_0d    :  overlap between the Coulomb function and the continuum
                  orbital as computed from numerical integration
    olap_dd    :  overlap of the continuum orbital with itself as computed
                  from numerical integration
    rmin,rmax  :  upper and lower bound for integration over
                  radial coordinate

    Returns
    =======
    scale      :  scale factor
    delta      :  phase shift
    """
    k = np.sqrt(2.0*E)
    # Finding the scale factor and phase shift that satisfy simultaneously
    #       scale * I_0d(delta) = olap_0d
    #  and
    #       scale^2 * I_dd(delta) = olap_dd
    # is equivalent to finding the root delta0 of
    # 
    #                          2                                         2 
    #       f(delta0) = olap_0d  * I_dd(delta0) -  olap_dd * I_0d(delta0)
    #
    # The scale factor is then determined as
    #                  olap_0d                              olap_dd
    #       scale = ------------       or   scale = sqrt(------------)
    #               I_0d(delta0)                         I_dd(delta0)
    #
    def f(delta):
        I_0d, I_dd = phaseshift2overlaps(E, Z, l, m, delta, rmin, rmax)
        return  olap_0d**2 * I_dd - olap_dd * I_0d**2 

    debug = 0
    if debug > 0:
        # Where does f(delta) cross the abscissa?
        import matplotlib.pyplot as plt
        deltas = np.linspace(0.0, 2.0*np.pi, 1000)
        fs = [f(d) for d in deltas]
        I0ds = [phaseshift2overlaps(E, Z, l, m, d, rmin, rmax)[0] for d in deltas]
        Idds = [phaseshift2overlaps(E, Z, l, m, d, rmin, rmax)[1] for d in deltas]
        plt.xlabel(r"$\delta$")
        # 
        plt.plot(deltas, fs, label=r"$f(\delta)$")
        plt.plot(deltas, I0ds, label=r"$I_{0d}$")
        plt.plot(deltas, Idds, label=r"$I_{dd}$")
        # show x-axis
        plt.plot([0,2.0*np.pi], [0,0], lw=2)
        plt.legend()
        plt.show()

    #
    # Since approximately
    #     I_0d(delta) ~ 1/2 k (rmax-rmin) cos(delta)
    # and
    #     I_dd(delta) ~ 1/2 k (rmax-rmin)
    # an initial guess for the phase shift is
    #                                       olap_0d^2
    #     delta = arccos(+/- sqrt( ------------------------- ) )
    #                              olap_dd 1/2 k (rmax-rmin)
    # 
    cos = np.sqrt( olap_0d**2 /(olap_dd * 0.5 * k*(rmax-rmin)) )
    assert -1.0 <= cos <= 1.0
    delta_guess = np.arccos(cos)
    print "Initial guess for phase shift  delta = %e" % delta_guess

    # find root of f(delta) with Newton's method
    try:
        delta_root = optimize.newton(f, delta_guess)
    except RuntimeError as e:
        msg  = "\nHINT:  You probably have to increase the resolution of\n"
        msg += "       the multicenter grid. If the computed radial wavefunction\n"
        msg += "       is not accurate enough at large radii it can happen that f(delta)\n"
        msg += "       has no roots at all!\n\n"
        print msg
        raise e
        
    print "Optimized phase shift  delta = %e" % delta_root

    # If f(delta)=0 then there should be 3 other roots in the interval [0,2*pi], which
    # are equivalent:
    #     f(pi-delta)=0
    #     f(pi+delta)=0
    # and f(2*pi-delta)=0
    # An additional phase shift of pi is compensated by a sign change in the scale factor,
    # so (delta,scale) and (delta+pi,-scale) give the same wavefunction.
    # Let's find the other roots and compare with their predicted positions.
    delta_root1 = delta_root
    delta_root2 = optimize.newton(f, np.pi-delta_root)
    delta_root3 = optimize.newton(f, np.pi+delta_root)
    delta_root4 = optimize.newton(f, 2*np.pi-delta_root)
    #
    print " "
    print "All roots of f(delta)"
    print "====================="
    print " root 1 delta1 = %e                    " % (delta_root1)
    print " root 2 delta2 = %e     pi-delta   = %e" % (delta_root2, np.pi-delta_root1)
    print " root 3 delta3 = %e     pi+delta   = %e" % (delta_root3, np.pi+delta_root1)
    print " root 4 delta4 = %e     2*pi-delta = %e" % (delta_root4, 2*np.pi-delta_root1)
    print " "
    
    # compute scale factor
    I_0d_root, I_dd_root = phaseshift2overlaps(E, Z, l, m, delta_root, rmin, rmax)
    if (abs(I_0d_root) > 0.0):
        scale = olap_0d / I_0d_root
    else:
        scale = np.sqrt( olap_dd / I_dd_root )
    
    return scale, delta_root

def phaseshift_olap(atomlist, phi, E, Z, l, m, rmin, rmax,
                    center=(0.0, 0.0, 0.0)):
    """
    determine the scale factor and phaseshift of a continuum orbital `phi` relative
    to the regular Coulomb wavefunction

    Parameters
    ==========
    atomlist   :  list of tuples (Z,[x,y,z]) with atomic numbers
                  and positions
    phi        :  continuum orbital with unknown phase shift
    E          :  energy of continuum state
    Z          :  nuclear charge, Z=+1 for hydrogen
    l,m        :  angular quantum numbers, -l <= m <= l
    rmin,rmax  :  upper and lower bound for integration over
                  radial coordinate

    Optional
    ========
    center     :  tuple (x0,y0,z0) with position of the center of charge

    Returns
    =======
    scale      :  scale factor, 1/scale * phi(x,y,z) has asymptotically the same
                  amplitude as the Coulomb wavefunction
    delta      :  phase shift in radians
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # Asymptotically phi(x,y,z) approaches a shifted Coulomb function
    # with unknown scale factor and phaseshift delta:
    #
    #  phi(x,y,z; delta) ---->  scale * Cf(x,y,z; delta)
    #
    # Since we know the overlap of two Coulomb waves as a function of
    # the phase shift, we can deduce the scale factor and the phase shift
    # of the orbital `phi` from the integrals
    #
    #          rmax    /rmax       /
    #  <Cf|phi>     =  |    r^2 dr | dOmega  Cf(x,y,z) phi(x,y,z; delta)
    #          rmin    /rmin       /
    # 
    #                  /rmax       /
    #    ---->  scale  |    r^2 dr | dOmega  Cf(x,y,z) Cf(x,y,z; delta)
    #                  /rmin       /
    #
    #         = scale * I_0d(delta)
    #
    # and 
    # 
    #          rmax    /rmax       /
    # <phi|phi>      = |    r^2 dr | dOmega  phi(x,y,z; delta) phi(x,y,z; delta)
    #          rmin     /rmin       /
    #
    #                  /rmax       /
    #   ----> scale^2  |    r^2 dr | dOmega  Cf(x,y,z; delta) Cf(x,y,z; delta)
    #                  /rmin       /
    #
    #       = scale^2 * I_dd(delta)
    #
    
    # The window selects points with radii in the interval [rmin,rmax]

    # regular coulomb function
    Cf = regular_coulomb_func(E, Z, l, m, 0.0, center=center)

    # origin
    x0,y0,z0 = center
    
    def window(x,y,z):
        # shift coordinates to new origin
        x,y,z = x-x0, y-y0, z-z0
        # distance from origin
        r = np.sqrt(x*x+y*y+z*z)
        w = 0.0*r
        w[(rmin <= r) & (r <= rmax)] = 1.0
        return w

    # To check whether the grid is fine enough, we compute the volume of
    # the window numerically and compare with exact expression:
    #             /rmax        /
    #    volume = |     r^2 dr | dOmega = 4/3 pi (rmax^3 - rmin^3)
    #             /rmin        /
    shell_volume = integral(atomlist, window)

    shell_volume_exact = 4.0/3.0 * np.pi * (rmax**3 - rmin**3)
    print "numerical volume integral = %e" % shell_volume
    print "exact volume integral     = %e" % shell_volume_exact
    print "relative error in volume  = %e %%" % (abs(shell_volume - shell_volume_exact)/shell_volume_exact * 100)
    
    # overlap between Coulomb function and continuum orbital
    def integrand_0d(x,y,z):
        # Our wavefunctions are real so we don't need complex conjugation in the bra.
        return Cf(x,y,z) * phi(x,y,z) * window(x,y,z)

    olap_0d = integral(atomlist, integrand_0d)

    print "overlap <Coulomb|continuum> over radial iterval [%s,%s] = %s" % (rmin, rmax, olap_0d)

    # overlap of continuum orbital with itself
    def integrand_dd(x,y,z):
        # Our wavefunctions are real so we don't need complex conjugation in the bra.
        return phi(x,y,z)**2 * window(x,y,z)

    olap_dd = integral(atomlist, integrand_dd)

    print "overlap <continuum|continuum> over radial interval [%s,%s] = %s" % (rmin, rmax, olap_dd)

    # Solve the system of equations
    #    <Couldomb|continuum>  = scale   * I_0d(delta)
    #    <continuum|continuum> = scale^2 * I_dd(delta)
    # for the scale factor `scale` and the phase shift `delta`.
    scale, delta = overlaps2phaseshift(E, Z, l, m, olap_0d, olap_dd, rmin, rmax)

    return scale, delta

def phaseshift(atomlist, phi, E, charge, l, m,
               center=(0.0, 0.0, 0.0),
               method="2pts"):
    """
    Determine the normalization factor and phaseshift for a continuum
    orbital `phi` such that asymptotically it is matched to the
    amplitude and phase of a Coulomb wave with energy E and angular
    momentum quantum numbers (l,m).

    Parameters
    ==========
    atomlist   :  list of tuples (Z,[x,y,z]) with atomic numbers
                  and positions
    phi        :  callable, phi(x,y,z) evaluates the unnormalized
                  continuum wavefunction
    E          :  energy of continuum state
    charge     :  charge Z defining asymptotic Coulomb potential,
                    V0(r) = -Z/|r-center|
    l,m        :  angular quantum numbers of asymptotic solution, 
                  -l <= m <= l
    
    Optional
    ========
    center     :  tuple (x0,y0,z0), origin of Coulomb wave
    method     :  "2pts", "olap" or "lstsq", determines which algorithm 
                  should be used for matching the continuum solution
                  to the Coulomb wave.

    Returns
    =======
    delta      :  phase shift (in radians)
    phi_norm   :  callable phi_norm(x,y,z), continuum orbital which asymptotically
                  has the same amplitude as a Coulomb wave with phase shift `delta`
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    
    # length of wave vector
    k = np.sqrt(2*E)
    wavelength = 2.0*np.pi / k
    print "wavelength = %e" % wavelength

    if method == "olap":
        # The phase shift is determined by integrating over
        # the radial coordinate in the interval [rmin,rmax].
        # On the one hand rmin < rmax should be chosen large enough,
        # so that the continuum orbital approaches its asymptotic form,
        # on the other hand rmax should be small enough that the accuracy
        # of the solution due to the sparse r-grid is still high enough.
        # A compromise has to be struck depending on the size of the radial grid.
        
        # We integrate over several periods, but not more than 30 bohr.
        rmin = 70.0 
        rmax = rmin + max(10*wavelength, 30.0)
            
        scale, delta = phaseshift_olap(atomlist, phi, E, charge, l, m, rmin, rmax,
                                       center=center)

    elif method == "lstsq":
        # The phase shift is determined by matching the radial wavefunction
        # to a shifted and scaled Coulomb function at a number of radial
        # sampling points drawn from the interval [rmin, rmax].
        # On the one hand rmin < rmax should be chosen large enough,
        # so that the continuum orbital approaches its asymptotic form,
        # on the other hand rmax should be small enough that the accuracy
        # of the solution due to the sparse r-grid is still high enough.
        # A compromise has to be struck depending on the size of the radial grid.
        
        # The matching points are spread over several periods,
        # but not more than 30 bohr.
        rmin = 70.0 
        rmax = rmin + max(10*wavelength, 30.0)
        Npts = 100

        # determine phase shift and scaling factor by a least square
        # fit the the regular Coulomb function
        scale, delta = phaseshift_lstsq(atomlist, phi, E, charge, l, m, rmin, rmax, Npts,
                                        center=center)

    elif method == "2pts":
        # The phase shift is determined from matching the radial wavefunction
        # to a shifted and scaled Coulomb function at two points r1 and r2.
        r1 = 50.0
        r2 = r1 + 0.01*wavelength
        
        scale, delta = phaseshift_2pts(atomlist, phi, E, charge, l, m, r1, r2,
                                       center=center)

    else:
        raise ValueError("method for phase matching should be 'olap', 'lstsq' or '2pts'")
        
        
    print "scale factor (relative to Coulomb wave) = %s" % scale
    print "phase shift (relative to Coulomb wave) = %e " % delta

    # normalize wavefunction, so that 1/scale phi(x,y,z) approaches
    # asymptotically a phase-shifted Coulomb wave
    phi_norm = multicenter_operation([phi], lambda fs: fs[0]/scale, 
                                     atomic_coordinates, atomic_numbers,
                                     radial_grid_factor=settings.radial_grid_factor,
                                     lebedev_order=settings.lebedev_order)

    return delta, phi_norm

######################################################################################
#
# Operator S for multicenter spherical averaging
# ----------------------------------------------
# The operator for averaing a function f spherically around the a single center R_i is
# given by
#                     /2pi   /pi
#  avg (f) = 1/(4 pi) | dph  |  sin(th ) dth  f(r ,th ,ph )
#     i               /0   i /0       i     i    i   i   i
#
# where (r_i,th_i,ph_i) are the spherical coordinates of the vector r - R_i .
#
# The multicenter spherical averaging operator is
#
#  S(f) = sum   avg (w f)
#            i     i  i 
#
# where w_i(r) is the weight function associated with the fuzzy Voronoi cell i.
#
# The operator S has the following properties:
#   1)  S is a linear, bounded operator with norm |S| <= 1
#   2)  |1-S| <= 1
#   3)  range(S) depends on the number and positions of the centers
#
######################################################################################

def spherical_remainder(atomlist, potential):
    """
    split the potential into two parts: the spherical average
    around each atom V^(sph) and the remainder V^(rem) = V - V^(sph)
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    spherical, remainder = multicenter_spherical_remainder(potential,
                                                   atomic_coordinates, atomic_numbers,
                                                   radial_grid_factor=settings.radial_grid_factor,
                                                   lebedev_order=settings.lebedev_order)

    return spherical, remainder

def solve_spherical_potential(atomlist, potential,
                              nmax=10):
    """
    given a potential function v, we solve

         S u = v

    for u, such that the multicenter spherical average of u is the desired potential v.

    To construct the inverse of S we iterate

       v  = v
        0
                                                   n
       v  = v    -  S v    = (1 - S) v    = (1 - S)  v
        n    n-1       n-1            n-1             0

    This series should converge to 0:

       v  ----> 0     for  n --> oo
        n

    The original potential may be reconstructed from the elements v_n:
             oo                 oo            oo         n
      v = sum     S v   =  S sum    v  = S sum    (1 - S)   v
             n=0     n          n=0  n        n=0

    which implies that 
             oo         n          oo  
      u = sum    (1 - S)  v  =  sum    v
             n=0                   n=0  n

    The inverse of S on the domain D can be found as the von Neumann series 
    (geometric series) provided the |1-S| < 1. We can only prove |1-S| <= 1, so
    we need to restrict the domain on which we can hope to find an inverse.
    The domain D is defined as the set of functions u for which (1-S)^n v converges.

       -1        oo         n          1           -1
      S    =  sum    (1 - S)   = -------------  = S
                 n=0             1 - ( 1 - S )


    Parameters
    ----------
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions, which  define the multicenter grid
    potential    :  callable, potential(x,y,z) evaluates v on the right hand side
                    of the equation  
                        S u = v

    Optional
    --------
    nmax         :  integer, maximum number of elements in the sequence
                    v_0,v_1,...,v_nmax which are constructed iteratively

    Returns
    -------
    norm_rem     :  norm of the last element, |v_{nmax}|. This number should be small, otherwise
                    nmax has to be increased or the potential v lies outside the function
                    space for which Su=v can be solved.
    u            :  callable, u(x,y,z) evaluates the solution of
                        S u = v
    """
    print "invert multicenter spherical averaging by solving  S.u = v ..."
    # By applying the operator (1-S) repeatedly to v, we construct the series
    # v_n and at the same time build up the solution u.

    print " "
    print "         Convergence of Remainder                                  "
    print " "
    print "            v_0 = v                                                "
    print "            v_n = (1-S) v_{n-1}                                    "
    print " "
    print "  Iteration           norm of remainder                  ratio     "
    print "      n                      |v_n|                |v_{n}|/|v_{n-1}|"
    print "  -----------------------------------------------------------------"
    
    # u = v0 + ...
    u = potential
    # |v_n| is undefined unless at least one iteration has been done.
    norm_remainder = None
    for n in range(0, nmax):
        #  potential = v_n
        #  spherical = S v_n
        #  remainder = v_{n+1} = (1 - S) v_n
        spherical, remainder = spherical_remainder(atomlist, potential)
        # |v_{n+1}| should converge to 0   for n --> oo
        norm_remainder = np.sqrt(integral(atomlist, lambda x,y,z: remainder(x,y,z)**2))

        # Print table that shows how the series |v_n| converges or doesn't converge.
        if n > 0:
            # Since ||(1-S)f|| <= ||f||, we should always have |v_{n}| <= |v_{n-1}|.
            # If this is not the case, this is a sign that the resolution of the grid
            # is not large enough.
            ratio = norm_remainder / norm_last_remainder

            msg = ""
            if ratio > 1.0:
                # Show a warning message
                msg = "(ratio exceeds 1, increase grid resolution!)"
            print "    %4.1d               %e               %e         %s"            % (n+1, norm_remainder, ratio, msg)

            if ratio > 1.0:
                # If the ratio exceeds 1, there is no reason to continue the iteration,
                # because we cannot expect to improve the potential futher, if the error
                # from the numerical integration is as large as the deviation we are trying
                # to correct.
                break
        else:
            print "    %4.1d               %e"            % (n+1, norm_remainder)

        # u += v_{n+1} 
        u = add_two_functions(atomlist, u, remainder, 1.0, 1.0)
            
        # prepare for next iteration
        # v_{n+1} --> v_{n}
        potential = remainder
        # keep |v_{n}|
        norm_last_remainder = norm_remainder
        
    print " "
    
    return norm_remainder, u

##############################################################
#
# Projection Method
#
##############################################################

def imaginary_time_propagation(atomlist, phi0, potential, e, t, n):
    """
    The solutions of the continuum Schroedinger equation

         (H - E) phi = 0

    are the same as those of
                2 
         (H - E)   phi = 0

    where E > 0 is given. `phi` can be obtained by propagation in imaginary
    time. The operator exp(-(H-E)^2 t) projects out the eigenfunctions of R^2 = (H-E)^2
    with eigenvalue 0, since all other eigenfunctions have eigenvalues > 0
    and decay.

          -(H-E)^2 t        t -->        0*t
         e           phi0 ----------->  e   phi  +  exponentially decaying terms

    The exponential is approximated by
                                         t   n
         exp(-(H-E)^2 t) = (1 - (H-E)^2 --- )
                                         n

    The initial guess phi0 should have the correct asymptotic behaviour.

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    phi0         :  callable, phi0(x,y,z) evaluates the initial guess for the wavefunction
    potential    :  callable, potential(x,y,z) evaluates effective potential V
    e            :  float, energy of continuum wavefunction
    t            :  float > 0, time of propagation
    n            :  integer, number of steps

    Returns
    =======
    phi          :  callable, phi(x,y,z) evaluates the final wavefunction
    """
    phi = phi0
    # time step
    dt = t/n
    
    for i in range(0, n):
        # compute (H-E) phi  = R phi
        residual  = residual_func(atomlist, phi, potential, e)
        # compute (H-E)^2 phi = R(R(phi))
        residual2 = residual_func(atomlist, residual, potential, e)
        # deviation from exact solution
        #  error^2 = ||(H-E)^2 phi||
        ### DEBUG
        import matplotlib.pyplot as plt
        r = np.linspace(-1.0, 1.0, 10000)
        x = 0*r
        y = 0*r
        z = r
        print residual2(x,y,z)
        plt.cla()
        plt.clf()
        plt.plot(r, residual(x,y,z), label="$(H-E) \phi$")
        plt.plot(r, residual2(x,y,z), label="$(H-E)^2 \phi$")
        plt.legend()
        plt.show()
        ###
        error = np.sqrt( integral(atomlist, residual2) )
        print " n= %d   |(H-E)^2 phi(%d)|= %e" % (i,i,error)
        #                       2
        # phi     =  ( 1 - (H-E)  dt ) phi
        #    n+1                          n
        phi = add_two_functions(atomlist, phi, residual2, 1.0, -dt)

    return phi


###############################################################################################
#
# TESTING
#
################################################################################################
from DFTB.MolecularIntegrals.Ints1e import overlap, nuclear_repulsion

def test_hydrogen_molecular_ion():
    """
    ground state energy of the hydrogen molecular ion (HMI, H2+)
    using Becke's basis-set-free solver for the Schroedinger equation.
    
    The initial guess for the 1\sigma_g orbital is the normalized sum
    of two 1s hydrogen orbitals centered on each proton.

    The numerically exact eigenenergies of the HMI for a bond length 
    of rHH = 2 bohr can be obtained with the program 
    'DFTB/Scattering/hydrogen_molecular_ion.py':

    1\sigma_g  Energy: -1.10263535012 Hartree
    2\sigma_g  Energy: -0.360864925127 Hartree
    3\sigma_g  Energy: -0.177680825727 Hartree
    1\sigma_u  Energy: -0.667534428815 Hartree
    2\sigma_u  Energy: -0.255413235739 Hartree
    3\sigma_u  Energy: -0.13730603945 Hartree
    1\pi_g     Energy: -0.226699702323 Hartree
    2\pi_g     Energy: -0.126703197559 Hartree
    1\delta_g  Energy: -0.212732490978 Hartree
    2\delta_g  Energy: -0.121011034108 Hartree
    
    """
    #xc = XCFunctionals.libXCFunctional("lda_x", "lda_c_xalpha")
    xc = XCFunctionals.XC_None()
    
    # bond length of H2+ (in bohr)
    rHH = 2.0
    atomlist = [(1, (0.0, 0.0, 0.0)),
                (1, (0.0, 0.0, rHH))]
    
    # initial guess for orbital
    def hydrogen_1s(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        psi = 2.0/np.sqrt(4.0*np.pi) * np.exp(-r)
        return psi
    
    def phi_1sg_unnorm(x,y,z):
        """1sigma_g without normalization"""
        psiA = hydrogen_1s(x,y,z)
        psiB = hydrogen_1s(x,y,z-rHH)
        psiApB = psiA + psiB
        return psiApB

    # normalize 1sigma_g orbital
    s = overlap(atomlist, phi_1sg_unnorm, phi_1sg_unnorm)

    def phi_1sg(x,y,z):
        """normalized 1sigma_g orbital"""
        psiA = hydrogen_1s(x,y,z)
        psiB = hydrogen_1s(x,y,z-rHH)
        psiApB = (psiA + psiB)/np.sqrt(s)
        return psiApB
        
    # initial guess
    orbitals = [phi_1sg]

    max_iter = 10
    thresh = 1.0e-4

    ####
    # plot cut through orbital along the z-axis
    import matplotlib.pyplot as plt

    Npts = 100
    zero = np.zeros(Npts)
    r = np.linspace(-2*rHH,3*rHH, Npts)

    fig, ax = plt.subplots(1,3)
    ####

    Enuc = nuclear_repulsion(atomlist)
    print "nuclear repulsion energy = %e" % Enuc
    
    for k in range(0, max_iter):
        rho = density_func(orbitals)
        veff_new = effective_potential_func(atomlist, rho, xc, nelec=1)
        if k > 0:
            #  (next)        (old)        (new)
            # V       = 2/3 V      + 1/3 V
            veff = add_two_functions(atomlist,
                                     veff, veff_new, 2.0/3.0, 1.0/3.0)
        else:
            veff = veff_new

        total_error = 0.0
        for i,phi in enumerate(orbitals):

            ####
            # cut along x-axis
            ax[0].set_xlabel("x / bohr")
            ax[0].set_ylabel(r"$\phi(x,0,0)$")
            linex, = ax[0].plot(r, phi(r,zero,zero), label=r"$\phi$   iter. %d" % k)

            # cut along y-axis
            ax[1].set_xlabel("y / bohr")
            ax[1].set_ylabel(r"$\phi(0,y,0)$")
            liney, = ax[1].plot(r, phi(zero,r,zero), label=r"$\phi$   iter. %d" % k)

            # cut along z-axis
            ax[2].set_xlabel("z / bohr")
            ax[2].set_ylabel(r"$\phi(0,0,z)$")
            linez, = ax[2].plot(r, phi(zero,zero,r), label=r"$\phi$   iter. %d" % k)

            #plt.show()
            
            ####

            # compute corrections for each orbital
            
            # energy expectation value <phi|H|phi>
            e = energy(atomlist, phi, phi, veff)

            residual = residual_func(atomlist, phi, veff, e)
            delta_e = energy_correction(atomlist, residual, phi)
            source = source_func(delta_e, phi, residual)

            print "    orbital= %d    energy E= %e   delta E= %e" % (i+1, e, delta_e)
            
            # orbital correction
            delta_phi = orbital_correction(atomlist, veff, source,  e, delta_e)
            total_error += overlap(atomlist, delta_phi, delta_phi)

            ####
            # cut along x-axis
            ax[0].plot(r, delta_phi(r,zero,zero), ls="-.", color=linex.get_color(), label=r"$\Delta \phi$ iter. %d" % k)
            ax[0].legend()

            # cut along y-axis
            ax[1].plot(r, delta_phi(zero,r,zero), ls="-.", color=liney.get_color(), label=r"$\Delta \phi$ iter. %d" % k)
            ax[1].legend()

            # cut along z-axis
            ax[2].plot(r, delta_phi(zero,zero,r), ls="-.", color=linez.get_color(), label=r"$\Delta \phi$ iter. %d" % k)
            ax[2].legend()

            #plt.show()
            
            ####
            
            # update orbital by mixing with orbital correction
            a,b = variational_mixture(atomlist, phi, delta_phi, veff)
            # improved orbital
            #    (new)        (old)              
            # phi      = a phi        +  b delta_phi
            #
            phi = add_two_functions(atomlist, phi, delta_phi, a, b)

            # check normalization
            norm2 = overlap(atomlist, phi, phi)
            print "<phi|phi>= %e" % norm2
                
            orbitals[i] = phi
            
        print " %3.1d   error = %e (threshold = %e)" % (k, total_error, thresh)
        
        if total_error < thresh:
            print "Converged"
            break

    print "final ground state energy E= %e Hartree" % e
    
    plt.show()

def test_h3_ion():
    """
    Reproduce the H3+ ground state calculation from Ref. [3].

    The final Hartree-Fock energy should be -1.30040 Hartree.
    """
    # equilateral H3+ geometry with a bond length of rHH=1.6405 bohr
    atomlist = [
        (1, ( 0.623730953331291,  1.301681603473471, 0.000000000000000)),
        (1, (-0.210837994662163, -0.110669251391040, 0.000000000000000)),
        (1, ( 1.429589880195322, -0.127247798591583, 0.000000000000000))]
    charge = +1
    # Hartree-Fock theory 
    xc = None

    # choose resolution of multicenter grids
    settings.radial_grid_factor = 3      # controls size of radial grid 
    settings.lebedev_order = 23          # controls size of angular grid
    
    RDFT = BasissetFreeDFT(atomlist, xc, charge=charge)
    RDFT.solveKohnSham()

def test_lithium_cation():
    """
    Try to reproduct the total HF energy of the Li^+ cation 
    which only has a closed 1s shell.
    """
    # Li^+ atom
    atomlist = [(3, (0.0, 0.0, 0.0))]
    charge = +1

    for rfac in [3, 6, 9, 12, 15, 18, 21, 24, 27, 30]:
        # choose resolution of multicenter grids
        settings.radial_grid_factor = rfac      # controls size of radial grid  
        settings.lebedev_order = 25           # controls size of angular grid
        # show number of radial and angular points in multicenter grid
        grid_sizes = print_grid_summary(atomlist,
                                        settings.lebedev_order, settings.radial_grid_factor)
    
        # 1s core orbitals for Li+^ atom
        RDFT = BasissetFreeDFT(atomlist, None, charge=charge)
        Etot, orbitals, orbital_energies = RDFT.solveKohnSham()

        Nr,Nang = grid_sizes[0]
        print "GRID= %d x %d        total energy = %e" % (Nr,Nang, Etot)

        
    
def test_h2o():
    import warnings
    warnings.filterwarnings("error")
    
    # experimental geometry of water
    #  r(OH) = 0.958 Ang, angle(H-O-H) = 104.4776 degrees
    atomlist = [
        (8, (0.000000000000000,  0.000000000000000, -0.222540557483415)),
        (1, (0.000000000000000, +1.431214118579765,  0.886071388908105)),
        (1, (0.000000000000000, -1.431214118579765,  0.886071388908105))]
    charge = 0    
    # LDA functional
    #xc = XCFunctionals.libXCFunctional('lda_x', 'lda_c_vwn_3')
    xc = XCFunctionals.libXCFunctional('gga_x_pbe', 'gga_c_pbe')
    
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 10      # controls size of radial grid 
    settings.lebedev_order = 23          # controls size of angular grid
    
    RDFT = BasissetFreeDFT(atomlist, xc, charge=charge)
    #RDFT.solveKohnSham()
    RDFT.solveKohnSham_new()

def test_lithium_dimer():
    
    # experimental bond length of Li2
    #  r(Li-Li) = 2.6729 Ang = 5.051048594872604 bohr
    rLiLi = 5.051048594872604
    atomlist = [
        (3, (0.0, 0.0, -0.5*rLiLi)),
        (3, (0.0, 0.0, +0.5*rLiLi))]
    charge = 0    
    # LDA functional
    xc = XCFunctionals.libXCFunctional('lda_x', 'lda_c_vwn_3')

    # choose resolution of multicenter grids
    settings.radial_grid_factor = 3      # controls size of radial grid 
    settings.lebedev_order = 23          # controls size of angular grid
    
    RDFT = BasissetFreeDFT(atomlist, xc, charge=charge)
    RDFT.solveKohnSham()

def test_h2():
    
    # experimental bond length of H2 (from NIST)
    #  r(H-H) = 0.74144 Ang = 1.401118436971957 bohr
    rHH = 1.401118436971957
    atomlist = [
        (1, (0.0, 0.0, -0.5*rHH)),
        (1, (0.0, 0.0, +0.5*rHH))]
    charge = 0    
    # LDA functional
    xc = XCFunctionals.libXCFunctional('lda_x', 'lda_c_vwn_3')

    # choose resolution of multicenter grids
    settings.radial_grid_factor = 3      # controls size of radial grid 
    settings.lebedev_order = 23          # controls size of angular grid
    
    RDFT = BasissetFreeDFT(atomlist, xc, charge=charge)
    RDFT.solveKohnSham2()
    
def test_coulomb_waves():
    """
    check that Coulomb waves phi(r) satisfy Schroedinger's equation
    for the hydrogen atom,
              __2
        (-1/2 \/  -  1/r) phi = E phi

    by plotting the residual   R(r) = (H-E) phi(r).

    At some point r > 20.0 the residual R(r) will deviate from 0 because
    the multicenter grid becomes very sparse for large r, so that the 
    numerical derivatives are not accurate.
    """    
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 40      # controls size of radial grid 
    settings.lebedev_order = 23          # controls size of angular grid
    # proton at the origin
    atomlist = [(1, (0.0, 0.0, 0.0))]

    # atomic number
    Z = 1
    # energy of the continuum state
    E = 1.0
    # wave vector
    k = np.sqrt(2*E)
    # angular quantum numbers of asymptotic continuum orbital
    # l=1, m=+1 corresponds to a px-orbital
    l = 1      
    m = +1
    # no phase shift 
    delta_l = 0.0

    # define regular coulomb function
    Cf = regular_coulomb_func(E, Z, l, m, 0.0)

    # There are no bound occupied orbitals, so the density is rho=0
    rho = density_func([])
    # potential energy, just nuclear attraction  -1/r
    potential = effective_potential_func(atomlist, rho, None, nelec=1)
    # residual R(r) = (H-E) Cf(r) 
    residual = residual_func(atomlist, Cf, potential, E)    

    # plot regular Coulomb wavefunction and residual along x-axis

    import matplotlib.pyplot as plt
    
    r = np.linspace(1.0e-3, 50.0, 1000)    
    plt.xlabel("r / bohr")
    
    plt.plot(r[r>0.5], potential(r[r>0.5],0*r[r>0.5],0*r[r>0.5]), label="V(r) = -1/r")
    plt.plot(r, Cf(r,0*r,0*r), label=r"regular Coulomb function C_f(r,0,0)")
    plt.plot(r, residual(r,0*r,0*r), label=r"residual $(H-E) C_f(r)$")
    plt.legend()
    
    plt.show()

def test_scattering_solution():
    """
    solve the Schroedinger equation for a continuum state 

    Assume that the Hamiltonian H = H0 + V can be split into a part 
                   __2
         H0 = -1/2 \/  + V0(r)

    for which the solutions phi0(r) are known, (H0-E)phi0 = 0, 
    and a perturbation V1(r) with V1 ---> 0 exponentially fast 
    for r ---> oo.

    The solutions for H consist of phi0 plus a correction dphi.
    With this ansatz the Schroedinger equation

         (H0 + V1 - E)(phi0 + dphi) = 0

    may be converted into an inhomogeneous Schroedinger equation:

         (H0 + V1 - E)dphi = -V1 phi0     (*)

    which has to be solved for dphi, given the source term -V1 phi0.

    As an example we take the potential of the hydrogen atom
     
         V0 = - 1/r 

    and add the perturbation 
         
         V1(r) = -exp(-r^2)

    The inhomogeneous Schroedinger equation (*) is solved on Becke's 
    multicenter grid numerically with the boundary condition 

         dphi(r->oo) = 0

    This is not the correct boundary condition for continuum states,
    but the solution is still correct (How come?).
    """        
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 120      # controls size of radial grid 
    settings.lebedev_order = 41          # controls size of angular grid
    
    # hydrogen atom
    atomlist = [(1, (0.0, 0.0, 0.0))]
    Z = 1

    # show number of radial and angular points in multicenter grid
    print_grid_summary(atomlist,
                       settings.lebedev_order, settings.radial_grid_factor)
    
    # The continuum orbital is specified by its energy and asymptotic
    # angular momentum (E,l,m)
    
    # energy of continuum orbital
    E = 1.0  
    # length of wave vectors
    k = np.sqrt(2*E)
    # angular quantum numbers of asymptotic solution  (px-orbital)
    l = 1
    m = +1

    # regular coulomb function
    Cf = regular_coulomb_func(E, Z, l, m, 0.0)
    
    # There are no bound orbitals, so the density of H+ is rho=0
    rho = density_func([])
    
    # V0(r) = -1/r 
    v0 = effective_potential_func(atomlist, rho, None, nelec=1)

    # perturbation V1(r) = -exp(-r*^2)
    def v1(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return -np.exp(-r**2)

    # combined potential V0 + V1
    def potential(x,y,z):
        return v0(x,y,z) + v1(x,y,z)

    # asymptotically correct solution for H0
    phi0 = Cf

    # right-hand side of inhomogeneous Schroedinger equation
    def source(x,y,z):
        return -v1(x,y,z) * phi0(x,y,z)

    #
    #  solve  (H0 + V1 - E) dphi = - V1 phi0
    # for orbital correction dphi

    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    
    dphi = multicenter_inhomogeneous_schroedinger(potential, source,  E,
                                                  atomic_coordinates, atomic_numbers,
                                                  radial_grid_factor=settings.radial_grid_factor,
                                                  lebedev_order=settings.lebedev_order)

    # Combine asymptotically correct solution with correction 
    #   phi = phi0 + dphi
    phi = add_two_functions(atomlist, phi0, dphi, 1.0, 1.0)

    # residual for phi0    R0 = (H-E)phi0
    residual0 = residual_func(atomlist, phi0, potential, E)
    # residual for final solution  R = (H-E)phi
    residual  = residual_func(atomlist, phi, potential, E)

    # spherical average of residual function
    residual_avg = spherical_average_residual_func(atomlist[0], residual)

    # Since phi() should be a Coulomb wave, its phase shift should
    # be zero. The phase shift is determined by integrating over
    # the radial coordinate in the interval [rmin,rmax].
    # On the one hand rmin < rmax should be chosen large enough,
    # so that the continuum orbital approaches its asymptotic form,
    # on the other hand rmax should be small enough that the accuracy
    # of the solution due to the sparse r-grid is still high enough.
    # A compromise has to be struck depending on the size of the radial grid.

    # We integrate over one period, but not more than 10 bohr.
    wavelength = min(2.0*np.pi / k, 10.0)
    print "wavelength = %e" % wavelength
    rmin = 10.0 
    rmax = rmin + wavelength

    scale, delta = phaseshift_olap(atomlist, phi, E, Z, l, m, rmin, rmax)
    print "scale factor (relative to Coulomb wave) = %s" % scale
    print "phase shift (relative to Coulomb wave) = %e " % delta
    
    # normalize wavefunction, so that 1/scale phi(x,y,z) approaches
    # asymptotically a phase-shifted Coulomb wave
    phi_norm = multicenter_operation([phi], lambda fs: fs[0]/scale, 
                                     atomic_coordinates, atomic_numbers,
                                     radial_grid_factor=settings.radial_grid_factor,
                                     lebedev_order=settings.lebedev_order)

    # shifted regular coulomb function
    Cf_shift = regular_coulomb_func(E, Z, l, m, delta)
    
    # plot wavefunction, potential and residual along x-axis
    r = np.linspace(1.0e-3, 50.0, 1000)

    import matplotlib.pyplot as plt
    plt.ylim((-0.55, 0.25))
    plt.xlabel("r / bohr")

    plt.plot(r, phi0(r,0*r,0*r), label=r"$\phi_0$")
    plt.plot(r, dphi(r,0*r,0*r), label=r"$\Delta \phi$")
    plt.plot(r, phi(r,0*r,0*r), label=r"$\phi = \phi_0 + \Delta \phi$")
    plt.plot(r, phi_norm(r,0*r,0*r), label=r"normalized $\phi$")
    plt.plot(r, Cf_shift(r,0*r,0*r), ls="-.", label=r"shifted Coulomb")

    # 
    plt.plot(r, potential(r,0*r,0*r), ls="-.", label=r"potential $V(r) = V_0(r) + V_1(r)$")
    # 
    plt.plot(r, residual0(r,0*r,0*r), label=r"residual $(H-E)\phi_0$")
    plt.plot(r, residual( r,0*r,0*r), label=r"residual $(H-E)\phi$")
    plt.plot(r, residual_avg(r), label=r"residual $\int d\Omega (H-E)\phi$")
    plt.legend()
    
    plt.show()

def test_phaseshift_overlap():
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 60      # controls size of radial grid 
    settings.lebedev_order = 41           # controls size of angular grid
    
    # hydrogen atom
    atomlist = [(1, (0.0, 0.0, 0.0))]
    Z = 1

    # show number of radial and angular points in multicenter grid
    print_grid_summary(atomlist,
                       settings.lebedev_order, settings.radial_grid_factor)
    
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    
    # The continuum orbital is specified by its energy and asymptotic
    # angular momentum (E,l,m)
    
    # energy of continuum orbital
    E = 1.0
    # length of wave vectors
    k = np.sqrt(2*E)
    # angular quantum numbers of asymptotic solution 
    l = 1
    m = +1

    # regular coulomb function
    Cf = regular_coulomb_func(E, Z, l, m, 0.0)

    # phase-shifted Coulomb function
    delta_l = 1.5 #1.0
    Cf_shift = regular_coulomb_func(E, Z, l, m, delta_l)

    # A normalized solution of the Schroedinger equation on a finite
    # grid will not conform to the normalization for continuum states.
    # So usually at large radii the continuum orbital will differ by a scale factor
    # from the shifted Coulomb wavefunction.
    scale = -4.0
    phi = multicenter_operation([Cf_shift], lambda fs: scale*fs[0], 
                                atomic_coordinates, atomic_numbers,
                                radial_grid_factor=settings.radial_grid_factor,
                                lebedev_order=settings.lebedev_order)
    
    # We integrate
    #               /rmax       /
    # I(delta_l) =  |    r^2 dr | dOmega  Cf(x,y,z) Cf_shift(x,y,z; delta_l)
    #               /rmin       /
    # to get the overlap between two Coulomb functions as a function of the
    # phase shift delta_l. This should allow us to deduce the phase shift of
    # continuum orbital from the overlap with the Coulomb function.

    # We integrate over one period, but not more than 10 bohr.
    wavelength = min(2.0*np.pi / k, 10.0)
    print "wavelength = %e" % wavelength
    rmin = 10.0 
    rmax = rmin + wavelength
    
    # The window selects points with radii in the interval [rmin,rmax]
    def window(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        w = 0.0*r
        w[(rmin <= r) & (r <= rmax)] = 1.0
        return w

    def integrand_0d(x,y,z):
        # Our wavefunctions are real so we don't need complex conjugation in the bra.
        return Cf(x,y,z) * phi(x,y,z) * window(x,y,z)
    
    def integrand_dd(x,y,z):
        # Our wavefunctions are real so we don't need complex conjugation in the bra.
        return phi(x,y,z) * phi(x,y,z) * window(x,y,z)

    I_0d = integral(atomlist, integrand_0d)

    I_dd = integral(atomlist, integrand_dd)

    I_0d_ana, I_dd_ana = phaseshift2overlaps(E, Z, l,m, delta_l, rmin, rmax)

    print "Numerical overlaps"
    print " I_0d = %s" % I_0d
    print " I_dd = %s" % I_dd
    print "Analytical overlaps"
    print " I_0d = %s" % I_0d_ana
    print " I_dd = %s" % I_dd_ana

    # recover scaling and phas shift from overlaps
    scale0, delta0 = overlaps2phaseshift(E, Z, l, m, I_0d, I_dd, rmin, rmax)

    print " scale factor = %s" % scale0
    print " phase shift = %s" % delta0
    
    print " expected scale factor = %s" % scale
    print " expected phase shift  = %s" % delta_l

def test_asymptotic_regular_coulomb():
    """
    compare regular Coulomb wavefunction graphically with
    its asymptotic form

         1/r sin(kr - eta log(2 k r) - l pi/2 + sigma_l)

    where sigma_l is the Coulomb phase shift and eta = -Z/k.
    """
    # charge of core
    Z = +1
    # energy of continuum orbial
    E = 0.1
    k = np.sqrt(2*E)
    # asymptotic angular momentum quantum numbers, (pz-orbital)
    l = 1
    m = 0
    
    # regular coulomb function
    Cf = regular_coulomb_func(E, Z, l, m, 0.0)

    eta = -Z/k
    # Coulomb phase shift, sigma_l = arg Gamma(l+1+i*eta)
    g = mpmath.gammainc(l+1.0+1.0j*eta)
    sigma_l = np.angle( complex(g) )

    def asymptotic(x,y,z):
        # spherical coordinates
        r,th,ph = cartesian2spherical((x,y,z))
        # radial part of asymptotic wavefunction
        rho = k*r
        wfn_radial = np.sin(rho - eta*np.log(2*rho) - l*np.pi/2.0 + sigma_l)/r
        # angular part of wavefunction, Y^(real)_{l,m}(th,ph)
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,ll,mm in sph_it:
            if ll == l and mm == m:
                # real spherical harmonics
                if m < 0:
                    Ylm_real = -np.sqrt(2.0) * Ylm.imag
                elif m > 0:
                    Ylm_real =  np.sqrt(2.0) * (-1)**m * Ylm.real
                else:
                    # m == 0
                    Ylm_real = Ylm.real
                    
                wfn_angular = Ylm_real
                break

        # combine radial and angular parts
        wfn = wfn_radial * wfn_angular

        return wfn

    import matplotlib.pyplot as plt
    plt.xlabel("z / bohr")

    r = np.linspace(0.003, 1000.0, 10000)
    plt.plot(r, Cf(0*r,0*r,r), label=r"Coulomb wavefunction")
    plt.plot(r, asymptotic(0*r,0*r,r), label=r"$\sin(kr - \eta \log(2 k r) - l/2 \pi + \sigma_l)/r$")
    plt.legend()

    plt.show()

def test_improved_asymptotic_regular_coulomb():
    """
    The regular Coulomb wavefunction is approximated asymptotically
    by the expression (see eqns. (31) and (32) in Ref. [1])

       Cf(rho) ----->  |L| sin(rho - eta*log(2*rho) - 1/2 pi l + sigma_l + arg(L))

    where

         L = 2F0(l+1+i*eta,-l+i*eta, 1/(2*i*rho))
    
    is related to the generalized hypergeometric function 2F0, that
    may be expressed as a power series

    References
    ----------
    [1] A. Burgess, "The Determination of Phases and Amplitudes of Wave Functions",
        Proc. Phys. Soc. 1963, vol. 81
    """
    from scipy import special
    
    Z = +1
    l = 1
    k = 0.5

    eta = -Z/k

    def poch(a,n):
        """
        rising factorial or Pochhammer symbol
         (a)_0 = 1
         (a)_n = a (a+1) (a+2) ... (a+n-1)
        """
        p = 1
        for m in range(0, n):
            p *= a+m
        return p

    """
    # test Pochhammer symbols
    for n in range(0, 10):
        print "%s   %s" % (special.poch(0.3, n), poch(0.3, n))
    """
 
    def coulombf(rho):
        cf = np.zeros(rho.shape)
        for i,rho_i in enumerate(rho):
            cf[i] = float(mpmath.coulombf(l, -Z/k, rho_i))
        return cf

    def asymptotic(rho, nmax=0):
        # Coulomb phase shift, sigma_l = arg Gamma(l+1+i*eta)
        g = mpmath.gammainc(l+1.0+1.0j*eta)
        sigma_l = np.angle( complex(g) )

        phi0 = rho - eta*np.log(2*rho) - 0.5*np.pi*l + sigma_l

        # correction terms
        a1 = l+1+1.0j*eta
        a2 = -l+1.0j*eta
        z = 1.0/(2.0j*rho)

        lamb = 0.0
        for n in range(0, nmax+1):
            lamb += poch(a1,n) * poch(a2,n) / special.factorial(n) * z**n
        
        return abs(lamb) * np.sin(phi0 + np.angle(lamb))

    import matplotlib.pyplot as plt
    plt.xlabel("r / bohr")
    plt.ylim((-10.0,10.0))

    r = np.linspace(0.003, 1000.0, 10000)
    rho = k*r
    for nmax in range(0, 5):
        plt.plot(r, asymptotic(rho,nmax=nmax), label=r"$|\Lambda_{%d}|\sin(kr - \eta \log(2 k r) - l/2 \pi + \sigma_l + \arg(\Lambda_{%d}))/r$" % (nmax, nmax))
    plt.plot(r, coulombf(rho), ls="-.", label=r"Coulomb wavefunction")

    plt.legend()

    plt.show()


    # Coulomb phase shift, sigma_l = arg Gamma(l+1+i*eta)
    g = mpmath.gammainc(l+1.0+1.0j*eta)
    sigma_l = np.angle( complex(g) )

    phi0 = rho - eta*np.log(2*rho) - 0.5*np.pi*l + sigma_l

    # correction terms
    a1 = l+1+1.0j*eta
    a2 = -l+1.0j*eta
    z = 1.0/(2.0j*rho)
    
    lamb = 0.0
    for n in range(0, nmax+1):
        lamb += poch(a1,n) * poch(a2,n) / special.factorial(n) * z**n

    delta = np.arcsin(coulombf(rho)/abs(lamb)) - phi0 - np.angle(lamb)
    delta = delta.real % (np.pi)

    plt.xlabel("r / bohr")
    plt.ylim((-10.0,10.0))
    
    plt.plot(r, delta, label=r"$\delta$ mod $\pi$")
    
    plt.legend()
    plt.show()

    
def test1_lithium_scattering_solution():
    """
    compute continuum orbital in the electrostatic potential of the Li^+ core
    """
    # Li^+ atom
    #atomlist = [(3, (0.0, 0.0, 0.0))]
    atomlist = [(3, np.random.rand(3))]

    # choose resolution of multicenter grids for bound orbitals
    settings.radial_grid_factor = 20      # controls size of radial grid  
    settings.lebedev_order = 25          # controls size of angular grid
    # 1s core orbitals for Li+^ atom
    RDFT = BasissetFreeDFT(atomlist, None, charge=+1)
    #bound_orbitals = RDFT.getOrbitalGuess()
    Etot, bound_orbitals, orbital_energies = RDFT.solveKohnSham()

    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 20 #120      # controls size of radial grid  
    settings.lebedev_order = 25 #41          # controls size of angular grid
    # show number of radial and angular points in multicenter grid
    print_grid_summary(atomlist,
                       settings.lebedev_order, settings.radial_grid_factor)
    
    print "electron density..."
    # electron density of two electrons in the 1s core orbital
    rho = density_func(bound_orbitals)

    print "total charges..."
    qnuc, qelec, qtot = total_charge(atomlist, rho)
    # number of electrons from integration
    nelec = int(np.round(-qelec))
    
    # We have to fix the center of the coordinate system relative
    # to which the asymptotic Coulomb waves are defined. The origin
    # R0 should be chosen such that the effective potential at large
    # distances appraches
    #                           charge
    #     V(r) --->  V0(r) = - -------
    #                           |r-R0|
    #        

    # In a neutral molecule the asymptotic behaviour of the effective
    # Kohn-Sham potential
    #
    #     V    (r)   ---> - 1/r               for r -> oo
    #      eff
    # is not the result of electrostatics (because the monopole of the
    # total charge distribution vanishes) but of the asymptotic behaviour
    # of the exchange functional
    #
    #     V        [rho](r)  --->  -1/r       for r -> oo
    #      exchange 

    #
    # The origin r0 such that
    #
    #     V_exchange ----> -1/|r-r0|
    #
    # is taken as the center of electronic charge in the HOMO
    #
    #    r0 = <HOMO|r|HOMO>
    # 
    #
    print ""
    print "  Origin"
    print "  ======"
    
    if nelec > 0:
        homo = bound_orbitals[-1]
        # density of the highest occupied molecular orbital
        def rho_HOMO(x,y,z):
            return homo(x,y,z)**2
        
        print "integrals over charge density of HOMO orbital..."
        center_nuc, center_homo, center_tot = center_of_charge(atomlist, rho_HOMO, verbose=0)

        print "  The center of the electron density in the HOMO"
        print "  is chosen as the the origin."
        center = center_homo
    else:
        center_nuc, dummy, dummy = center_of_charge(atomlist, None)

        print "  As there are no electrons, the center of the nuclear charge"
        print "  is chosen as the origin."
        center = center_nuc

    print "  The origin of the coordinate system relative to which"
    print "  asymptotic solutions are defined is placed at"
    print ""
    print "    center = %+e  %+e  %+e" % tuple(center)
    print ""
    
    print "NOTE: Two-electron systems are treated with Hartree-Fock theory."
    print "      For Li+ the effective HF potential goes as -2/r for large r!"
    print "      Therefore the continuum electron sees a core with charge +2."
    xc = None
    # THIS IS WEIRD, so the continuum electron sees Li^(2+) instead of Li^+???
    # It looks strange, but makes perfect sense:
    # The effective Hartree-Fock potential for Li^+
    #   Veff = -3/r + 0.5 * Vcoul(r)
    # approaches
    #     V0 = -2/r
    # asymptotically, since for a closed-shell 2-electron system Vcoul(r) -> 2/r.
    # Therefore we have to match the solution to a Coulomb wave with charge +2.
    # If we would use a different effective potential for bound (-2/r) and continuum (-1/r)
    # states, the states would not be solutions of the same Schroedinger equation,
    # which seems inconsistent.

    # Asymptotically the effective potential is Coulomb-like -(Q+1)/r . It's strength consists
    # of an electrostatic part (Q=0 for a neutral molecule, Q=+1 for a cation) and an
    # exchange part (+1, since v_xc[rho] -> -1/r for r -> oo)
    # if there are other electrons present.
    charge = qtot
    if nelec > 0:
        # add charge felt due to exchange
        charge = qtot + 1
    print " electrostatic charge  qtot = %+e" % qtot
    print " asymptotic charge          = %+e" % charge
    
    print "effective potential..."
    # potential energy for Li nucleus and 2 core electrons
    potential = effective_potential_func(atomlist, rho, xc, nelec=2)

    # origin of Coulomb waves and V0
    x0,y0,z0 = center
    def v0(x,y,z):
        # shift coordinates
        x,y,z = x-x0,y-y0,z-z0
        # distance from origin
        r = np.sqrt(x*x+y*y+z*z)
        return -charge/r

    def v1(x,y,z):
        return potential(x,y,z) - v0(x,y,z)
    
    # The continuum orbital is specified by its energy and asymptotic
    # angular momentum (E,l,m)
    
    # energy of continuum orbital
    E = 1.0  
    # length of wave vectors
    k = np.sqrt(2*E)
#    # angular quantum numbers of asymptotic solution  (s-orbital)
#    l = 0
#    m = 0
    # angular quantum numbers of asymptotic solution  (px-orbital)
    l = 1
    m = +1

    print " "
    print "  Asymptotic continuum wavefunction"
    print "  ================================="
    print "  energy           E= %e Hartree  ( %e eV )" % (E, E*AtomicData.hartree_to_eV)
    print "  wavevector       k= %e a.u." % k
    print "  angular moment   l= %d m= %+d" % (l,m)
    print " "
    
    # asymptotically correct solution for V0 = -1/r (hydrogen)
    Cf = regular_coulomb_func(E, charge, l, m, 0.0,
                              center=center)
    phi0 = Cf
    
    # right-hand side of inhomogeneous Schroedinger equation
    def source(x,y,z):
        return -v1(x,y,z) * phi0(x,y,z)

    #
    #  solve  (H0 + V1 - E) dphi = - V1 phi0
    # for orbital correction dphi

    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    print "Schroedinger equation ..."
    dphi = multicenter_inhomogeneous_schroedinger(potential, source,  E,
                                                  atomic_coordinates, atomic_numbers,
                                                  radial_grid_factor=settings.radial_grid_factor,
                                                  lebedev_order=settings.lebedev_order)

    # Combine asymptotically correct solution with correction 
    #   phi = phi0 + dphi
    phi = add_two_functions(atomlist, phi0, dphi, 1.0, 1.0)

    print "residuals (H-E)phi ..."
    # residual for phi0    R0 = (H-E)phi0
    residual0 = residual_func(atomlist, phi0, potential, E)
    # residual for final solution  R = (H-E)phi
    residual  = residual_func(atomlist, phi, potential, E)

    # spherical average of residual function

    # The radial grid is chosen for the atom with the largest atomic number
    Zmax = max([Z for (Z,pos) in atomlist])
    # The origin of the radial grid is placed at `center`
    atom_center = (Zmax, center)
    residual_avg = spherical_average_residual_func(atom_center, residual)
    """
    ### DEBUG
    # show radial wavefunction
    wfn_rad = radial_wave_func(atomlist, phi, l, m)
    import matplotlib.pyplot as plt
    r = np.linspace(1.0e-3, 200.0, 4000)
    plt.plot(r, wfn_rad(r), label=r"$R_{%d,%d}(r)$" % (l,m))
    plt.show()
    ###
    """

    print "phase shifts..."
    # Which algorithm should be used for matching the continuum solution
    # to the Coulomb wave to determine the phase shift and scaling factor?
    phase_matching = "2pts" # "2pts" "olap" "lstsq"  
    
    if phase_matching == "olap":
        # The phase shift is determined by integrating over
        # the radial coordinate in the interval [rmin,rmax].
        # On the one hand rmin < rmax should be chosen large enough,
        # so that the continuum orbital approaches its asymptotic form,
        # on the other hand rmax should be small enough that the accuracy
        # of the solution due to the sparse r-grid is still high enough.
        # A compromise has to be struck depending on the size of the radial grid.

        # We integrate over several periods, but not more than 30 bohr.
        wavelength = 2.0*np.pi / k
        print "wavelength = %e" % wavelength
        rmin = 70.0 
        rmax = rmin + max(10*wavelength, 30.0)
        
        scale, delta = phaseshift_olap(atomlist, phi, E, charge, l, m, rmin, rmax,
                                       center=center)

    elif phase_matching == "lstsq":
        # The phase shift is determined by matching the radial wavefunction
        # to a shifted and scaled Coulomb function at a number of radial
        # sampling points drawn from the interval [rmin, rmax].
        # On the one hand rmin < rmax should be chosen large enough,
        # so that the continuum orbital approaches its asymptotic form,
        # on the other hand rmax should be small enough that the accuracy
        # of the solution due to the sparse r-grid is still high enough.
        # A compromise has to be struck depending on the size of the radial grid.
        
        # The matching points are spread over several periods,
        # but not more than 30 bohr.
        wavelength = 2.0*np.pi / k
        print "wavelength = %e" % wavelength
        rmin = 70.0 
        rmax = rmin + max(10*wavelength, 30.0)
        Npts = 100

        # determine phase shift and scaling factor by a least square
        # fit the the regular Coulomb function
        scale, delta = phaseshift_lstsq(atomlist, phi, E, charge, l, m, rmin, rmax, Npts,
                                        center=center)

    elif phase_matching == "2pts":
        # The phase shift is determined from matching the radial wavefunction
        # to a shifted and scaled Coulomb function at two points r1 and r2.
        wavelength = 2.0*np.pi / k
        print "wavelength = %e" % wavelength

        r1 = 50.0
        r2 = r1 + 0.001*wavelength
        
        scale, delta = phaseshift_2pts(atomlist, phi, E, charge, l, m, r1, r2,
                                       center=center)

        
        
    print "scale factor (relative to Coulomb wave) = %s" % scale
    print "phase shift (relative to Coulomb wave) = %e " % delta

    # normalize wavefunction, so that 1/scale phi(x,y,z) approaches
    # asymptotically a phase-shifted Coulomb wave
    phi_norm = multicenter_operation([phi], lambda fs: fs[0]/scale, 
                                     atomic_coordinates, atomic_numbers,
                                     radial_grid_factor=settings.radial_grid_factor,
                                     lebedev_order=settings.lebedev_order)

    # The continuum orbital should be orthogonal to the bound
    # orbitals belonging to the same Hamiltonian. I think this
    # should come out correctly by default.
    print " "
    print "  Overlaps between bound orbitals and continuum orbital"
    print "  ====================================================="
    for ib,bound_orbital in enumerate(bound_orbitals):
        olap_bc = overlap(atomlist, bound_orbital, phi_norm)
        print "  <bound %d| continuum> = %e" % (ib+1, olap_bc)
    print ""

    # shifted regular coulomb function
    Cf_shift = regular_coulomb_func(E, charge, l, m, delta,
                                    center=center)
    
    # save radial wavefunctions and spherically averaged residual
    # radial part of Coulomb wave without phase shift
    phi0_rad = radial_wave_func(atomlist, phi0, l, m,
                                center=center)
    # radial part of shifted Coulomb wave
    Cf_shift_rad = radial_wave_func(atomlist, Cf_shift, l, m,
                                    center=center)
    # radial part of normalized scattering solution
    phi_norm_rad = radial_wave_func(atomlist, phi_norm, l, m,
                                    center=center)

    print "# RADIAL_WAVEFUNCTIONS"
    print "# Asymptotic wavefunction:"
    print "#   charge            Z= %+d" % charge
    print "#   energy            E= %e  k= %e" % (E, k)
    print "#   angular momentum  l= %d m= %+d" % (l, m)
    print "#   phase shift       delta= %e rad" % delta
    print "# "
    print "#     R/bohr                 Coulomb               Coulomb           radial wavefunction   spherical avg. residual"
    print "#                                                  shifted                R_{l,m}(r)           <|(H-E)phi|^2>"
    import sys
    # write table to console
    r = np.linspace(1.0e-3, 100, 1000)
    data = np.vstack((r, phi0_rad(r), Cf_shift_rad(r), phi_norm_rad(r), residual_avg(r))).transpose()
    np.savetxt(sys.stdout, data, fmt="    %+e    ")
    print "# END"


    
    # plot wavefunction, potential and residual along x-axis
    r = np.linspace(1.0e-3, 200.0, 4000)

    # save cut of wavefunction along x-axis to file
    data = np.vstack((r, phi0(r,0*r,0*r), phi_norm(r,0*r,0*r), residual_avg(r))).transpose()
    fh = open("/tmp/li+_wavefunction.dat", "w")
    print>>fh, "# phase shift (relative to Coulomb wave) = %e" % delta
    print>>fh, "# R / bohr      Coulomb         phi      <(H-E)phi>"
    np.savetxt(fh, data, fmt="%+e")
    fh.close()
    
    import matplotlib.pyplot as plt
    plt.ylim((-0.55, 0.25))
    plt.xlabel("r / bohr")

    # wavefunctions
    plt.plot(r, phi0(r,0*r,0*r), label=r"$\phi_0$")
    plt.plot(r, dphi(r,0*r,0*r), label=r"$\Delta \phi$")
    plt.plot(r, phi(r,0*r,0*r), label=r"$\phi = \phi_0 + \Delta \phi$")
    plt.plot(r, phi_norm(r,0*r,0*r), label=r"normalized $\phi$")
    plt.plot(r, Cf_shift(r,0*r,0*r), ls="-.", label=r"shifted Coulomb")

    # potentials
    plt.plot(r, v0(r,0*r,0*r), lw=2, ls="-.", label=r"potential $V_0(r) = -%d/r$" % int(charge))
    plt.plot(r, v1(r,0*r,0*r), lw=2, ls="-.", label=r"potential $V_1(r) = V(r) - V_0(r)$")
    plt.plot(r, potential(r,0*r,0*r), lw=2, ls="-.", label=r"potential $V(r) = V_0(r) + V_1(r)$")
    
    # residuals
    plt.plot(r, residual0(r,0*r,0*r), label=r"residual $(H-E)\phi_0$")
    plt.plot(r, residual( r,0*r,0*r), label=r"residual $(H-E)\phi$")
    plt.plot(r, residual_avg(r), label=r"residual $\int d\Omega (H-E)\phi$")
    
    plt.legend()
    
    plt.show()


def test2_lithium_scattering_solution():
    """
    compute continuum orbital in the electrostatic potential of the Li^+ core
    """
    # Li^+ atom
    atomlist = [(3, np.random.rand(3))]

    # choose resolution of multicenter grids for bound orbitals
    settings.radial_grid_factor = 20      # controls size of radial grid  
    settings.lebedev_order = 25          # controls size of angular grid
    # 1s core orbitals for Li+^ atom
    RDFT = BasissetFreeDFT(atomlist, None, charge=+1)
    bound_orbitals = RDFT.getOrbitalGuess()
    #Etot, bound_orbitals, orbital_energies = RDFT.solveKohnSham()

    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 20 #120      # controls size of radial grid  
    settings.lebedev_order = 25 #41          # controls size of angular grid

    # electron density of two electrons in the 1s core orbital
    rho = density_func(bound_orbitals)
    # HOMO wavefunction
    homo = bound_orbitals[-1]

    # energy of continuum orbital
    E = 1.0  
#    # angular quantum numbers of asymptotic solution  (s-orbital)
#    l = 0
#    m = 0
    # angular quantum numbers of asymptotic solution  (px-orbital)
    l = 1
    m = +1
    
    delta, phi = RDFT.solveScatteringProblem(rho, homo, E, l, m)

    
def test1_spherical_remainder():
    """
    repeatedly subtract the multicenter spherical averaged potential
    until the remaining potential is reduced to zero.

    We iterated
                           (sph)
           V    =  V   -  V   
            n+1     n      n

    where
            (sph)                   /
           V       =  sum  1/(4 pi) | (w * V )(r-R ) dOmega
            n            I          /   I   n     I         I

    is the spherical average of the potential around each atom I 
    (R_I is the position of the atom and w_I is the weight function 
     of the fuzzy Voronoi decomposition)

    until  V  ---> 0
            n

    We also need to check that the original potential can be 
    reconstructed from the spherical parts alone:

                    (sph)
          V = sum  V 
                 n  n
    """
    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    # center of charge
    center = (0.0, 0.0, 0.0)

    
    atomlist = [(1, posH1),
                (1, posH2)]
    charge = +2
    
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 1      # controls size of radial grid  
    settings.lebedev_order = 59          # controls size of angular grid
    
    rho = None
    xc = None
    print "effective potential..."
    potential = effective_potential_func(atomlist, rho, xc, nelec=0)

    # Because v0 has a singularity at the origin, we have to add
    # a grid around the origin by placing a dummy atom there.
    atomlist_grid = [(1, center),
                     (1, posH1),
                     (1, posH2)]
    
    import matplotlib.pyplot as plt
    r = np.linspace(-5.0, +5.0, 10000)

    fig, axes = plt.subplots(1,2)

    axes[0].set_ylabel("remainders")
    axes[1].set_ylabel("spherical parts")
    
    axes[0].plot(r, potential(0.0*r, 0.0*r, r), lw=2, label=r"$V(0,0,z)$")
    axes[1].plot(r, potential(0.0*r, 0.0*r, r), lw=2, label=r"$V(0,0,z)$")

    spherical_parts = []

    nmax = 5
    for n in range(0, nmax):
        print "V_{%d} = V_{%d} - V^(sph)_{%d}" % (n+1,n,n)
        spherical, remainder = spherical_remainder(atomlist_grid, potential)
        # check that spherical and remainder add up to the original potential
        #  V = V^(sph) + (1-V^(sph))
        def error(x,y,z):
            err = abs(potential(x,y,z) - spherical(x,y,z) - remainder(x,y,z))
            return err
        
        spherical_parts.append( spherical )
        l, = axes[0].plot(r, remainder(0.0*r, 0.0*r, r), label=r"remainder $V^(rem)_{%d}(0,0,z)$" % n)
        axes[0].plot(r, error(0.0*r, 0.0*r, r), label=r"error $\vert V_{%d} - (V^(sph)_{%d} + V^(rem)_{%d}) \vert$" % (n,n,n),
                     ls="-.", color=l.get_color())
        potential = remainder

    # reconstruct full potential from spherical parts
    v_reconstructed = 0.0*r
    for n,spherical in enumerate(spherical_parts):
        v_reconstructed += spherical(0.0*r, 0.0*r, r)

        axes[1].plot(r, v_reconstructed, label=r"reconstructed $\sum_{n=0}^{%d} V^{(sph)}_n$" % n)
        
    for ax in axes:
        ax.set_xlabel("z / bohr")
        ax.legend()
        
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.5)
    
    plt.show()

def test2_spherical_remainder():
    """
    Here we check that the decomposition of a function into a spherical average around
    the different centers and the remainder becomes exact as we increase the resultion
    of the angular grid:
    
        V = V^(rem) + V^(sph)

    The reconstruction error |V - (V^(rem) + V^(sph)| is computed as a function of the 
    number of grid points in the Lebedev grid.

    For the H2+ nuclear potential and a radial grid with radial_grid_factor=3 times more
    points than in Becke's default grid for hydrogen, the following table is obtained:

     Lmax     nr. of points          reconstruction
             in angular grid             error     
       3           6                 3.492031e+01       
       5          14                 5.450564e+00       
       7          26                 1.508499e+01       
       9          38                 7.858227e+00       
      11          50                 4.171945e+00       
      13          74                 5.853172e+00       
      15          86                 5.216042e+00       
      17         110                 1.547622e+00       
      19         146                 1.765951e+00       
      21         170                 1.183767e+00       
      23         194                 7.699948e-01       
      25         230                 1.202782e+00       
      27         266                 8.486140e-01       
      29         302                 2.586290e-01       
      31         350                 8.533220e-01       
      35         434                 1.268368e-01       
      41         590                 7.664443e-02       
      47         770                 3.502057e-02       
      53         974                 2.247854e-02       
      59        1202                 9.843761e-03       
      65        1454                 6.437152e-03       
      71        1730                 2.872481e-03       
      77        2030                 1.961795e-03       
      83        2354                 1.125463e-03       
     131        5810                 8.290255e-04 

    """
    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    # center of charge
    center = (0.0, 0.0, 0.0)

    
    atomlist = [(1, posH1),
                (1, posH2)]
    charge = +2
    
    rho = None
    xc = None
    print "effective potential..."
    potential = effective_potential_func(atomlist, rho, xc, nelec=0)
    
    # add a grid around the origin by placing a dummy atom there
    atomlist_grid = [(1, center),
                     (1, posH1),
                     (1, posH2)]

    # choose radial resolution of multicenter grid
    rfac = 3                              # increase density of radial points by a factor of 3
    settings.radial_grid_factor = rfac    # controls size of radial grid
    
    # How does the angular resolution of the grids affect the reconstruction error?
    from DFTB.MolecularIntegrals.lebedev_quad_points import Lebedev_Lmax, Lebedev_Npts
    # error(Lmax)
    errors = np.zeros(len(Lebedev_Lmax))
    for i,Lmax in enumerate(Lebedev_Lmax):
        settings.lebedev_order = Lmax           # controls size of angular grid

        # decompose potential in multicenter average plus remainder
        #  V = V^(sph) + V^(rem)
        spherical, remainder = spherical_remainder(atomlist_grid, potential)

        # compute the reconstruction error
        def deviation(x,y,z):
            V = potential(x,y,z)
            Vsph = spherical(x,y,z)
            Vrem = remainder(x,y,z)
            return V - (Vsph + Vrem)

        errors[i] = np.sqrt( integral(atomlist, lambda x,y,z: deviation(x,y,z)**2) )
        
        print ""
        print "resolution of angular grid:            Lmax = %d" % Lmax
        print "reconstruction error: |V-(V^(sph)+V^(rem))| = %e" % errors[i]

    print "#"
    print "# resolution of radial grid: rfac= %d" % rfac
    print "# "
    print "#   Lmax       nr. of points            reconstruction"
    print "#             in angular grid               error     "
    for Lmax,Npts,err in zip(Lebedev_Lmax, Lebedev_Npts, errors):
        print "  %6.1d      %6.1d                 %e       " % (Lmax,Npts,err)


    import matplotlib.pyplot as plt

    plt.xlabel(r"Lebedev order $L_{max}$")
    plt.ylabel(r"reconstruction error")

    plt.semilogy(Lebedev_Lmax, errors, "o")

    plt.show()
    
    
def test_solve_spherical_potential_hmi():
    """
    check that we can find a potential u such that multicenter spherical averaging
    of u gives the molecular potential v

        S u = v

    v is the potential of the H2^+ ion

        v = -1/|r-R/2*ez| - 1/|r+R/2*ez|
    """
    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    # center of charge
    center = (0.0, 0.0, 0.0)

    # list of real atoms
    atomlist_nuclei = [(1, posH1),
                       (1, posH2)]
    charge = +2
    
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 1       # controls size of radial grid  
    settings.lebedev_order = 59          # controls size of angular grid
    
    rho = None
    xc = None
    print "effective potential..."
    potential = effective_potential_func(atomlist_nuclei, rho, xc, nelec=0)

    # We can add a grid around the origin by placing a dummy atom there.
    atomlist = [(1, center),
                (1, posH1),
                (1, posH2)]

    v = potential
    # solve S.u = v
    nmax = 5
    norm_rem, u = solve_spherical_potential(atomlist, v, nmax=nmax)
    # reconstruct v from u
    #   v_rec = S.u
    v_rec, dummy = spherical_remainder(atomlist, u)
    
    # 
    import matplotlib.pyplot as plt
    r = np.linspace(-5.0, +5.0, 10000)

    fig, axes = plt.subplots(1,3)
    fig.suptitle(r"after $n_{max} = %d$ iterations" % nmax)
    
    axes[0].set_xlabel("x / bohr")
    axes[1].set_xlabel("y / bohr")
    axes[2].set_xlabel("z / bohr")

    # plot cuts of potentials along x, y and z-axis
    xcut = (  r,0*r,0*r)  # only x-coordinate changes y=z=0
    ycut = (0*r,  r,0*r)
    zcut = (0*r,0*r,  r)
    for iax,(x,y,z) in enumerate( [xcut, ycut, zcut] ):
        axes[iax].plot(r, v(x,y,z),     lw=2,          label="$V$")
        axes[iax].plot(r, u(x,y,z),     lw=2, ls="-.", label="$U = S^{-1} V$")
        axes[iax].plot(r, v_rec(x,y,z), lw=2,          label="reconstructed $V^{(rec)} = S U$")
        
        axes[iax].legend()
        
#    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5)
    
    plt.show()

def test_solve_spherical_potential_water():
    """
    check that we can find a potential u such that multicenter spherical averaging
    of u gives the molecular potential v

        S u = v

    v is the effective Kohn-Sham potential of the water molecule

    """
    # experimental geometry of water
    #  r(OH) = 0.958 Ang, angle(H-O-H) = 104.4776 degrees
    atomlist = [
        (8, (0.000000000000000,  0.000000000000000, -0.222540557483415)),
        (1, (0.000000000000000, +1.431214118579765,  0.886071388908105)),
        (1, (0.000000000000000, -1.431214118579765,  0.886071388908105))]
    charge = 0    
    # LDA functional
    xc = XCFunctionals.libXCFunctional('lda_x', 'lda_c_vwn_3')
    
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 10       # controls size of radial grid  
    settings.lebedev_order = 65           # controls size of angular grid

    print_grid_summary(atomlist,
                       settings.lebedev_order, settings.radial_grid_factor)

    
    RDFT = BasissetFreeDFT(atomlist, xc, charge=charge)
    # doubly occupied DFTB orbitals
    orbitals = RDFT.getOrbitalGuess()
    
    # total electron density
    print "electron density..."
    rho = density_func(orbitals)

    print "total charge..."
    qnuc, qelec, qtot = total_charge(atomlist, rho)

    print "effective potential..."
    potential = effective_potential_func(atomlist, rho, xc)

    v = potential
    # solve S.u = v
    nmax = 5
    norm_rem, u = solve_spherical_potential(atomlist, v, nmax=nmax)
    # reconstruct v from u
    #   v_rec = S.u
    v_rec, dummy = spherical_remainder(atomlist, u)

    # v_missing is the difference between the original potential
    # and the reconstructed potential, its norm is the reconstruction
    # error
    #   v_missing = v - v_rec
    #
    # In order to visualize the deviation from the original
    # potential, v_missing is saved to a cube file.
    def v_missing(grid, dV):
        (x,y,z) = grid
        return v(x,y,z) - v_rec(x,y,z)

    from DFTB.Analyse import Cube
    cube_file = "/tmp/water_lda_vmissing.cube"
    Cube.function_to_cubefile(atomlist, v_missing, filename=cube_file,
                              ppb=5.0)
    print "deviation V - V_rec saved to cube file '%s'" % cube_file

    
    # save tables with cuts along the y-axis of V, V^(rec) and U
    r = np.linspace(-5.0, +5.0, 10000)
    dat_file = "/tmp/water_lda_reconstructed_potential.dat"
    fh = open(dat_file, "w")
    # z-positions of cuts along y-axis, x=0 in all cuts
    for z in [-0.222540557483415, 0.0, 0.886071388908105]:
        data = np.vstack(( r,
                           v(0*r,r,0*r+z),
                           v_rec(0*r,r,0*r+z),
                           u(0*r,r,0*r+z))).transpose()
        print>>fh, "# cut through potential along the y-axis at x=0, z=%e" % z
        print>>fh, "# "
        print>>fh, "# R/bohr    V/Hartree       V^(rec)/Hartree      U/Hartree"
        np.savetxt(fh, data, fmt="%+10.8e")
        print>>fh, " "
    print "tables with V, V^(rec)=S.U and U=S^(-1).V  written to '%s'" % dat_file
    fh.close()

def test_solve_spherical_potential_water_dummy():
    """
    check that we can find a potential u such that multicenter spherical averaging
    of u gives the molecular potential v

        S u = v

    v is the effective Kohn-Sham potential of the water molecule. 

    Additional dummy atoms are added at the origin and below and above the molecular
    plane, so as to increase the resolution of the multicenter grid in this region.

    """
    # experimental geometry of water
    #  r(OH) = 0.958 Ang, angle(H-O-H) = 104.4776 degrees
    atomlist_nuclei = [
        (8, (0.000000000000000,  0.000000000000000, -0.222540557483415)),
        (1, (0.000000000000000, +1.431214118579765,  0.886071388908105)),
        (1, (0.000000000000000, -1.431214118579765,  0.886071388908105))]
    charge = 0    
    # LDA functional
    xc = XCFunctionals.libXCFunctional('lda_x', 'lda_c_vwn_3')

    # additional dummy atoms which add grid points but no nuclear or electronic charge
    atomlist_dummy = [
        (1, ( 0.000000000000000,  0.000000000000000, 0.000000000000000)),
        (1, (+1.000000000000000,  0.000000000000000, 0.000000000000000)),
        (1, (-1.000000000000000,  0.000000000000000, 0.000000000000000)) ]

    atomlist = atomlist_nuclei + atomlist_dummy
    
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 10       # controls size of radial grid  
    settings.lebedev_order = 65           # controls size of angular grid

    print_grid_summary(atomlist,
                       settings.lebedev_order, settings.radial_grid_factor)

    
    RDFT = BasissetFreeDFT(atomlist_nuclei, xc, charge=charge)
    # doubly occupied DFTB orbitals
    orbitals = RDFT.getOrbitalGuess()
    
    # total electron density
    print "electron density..."
    rho = density_func(orbitals)

    print "total charge..."
    qnuc, qelec, qtot = total_charge(atomlist_nuclei, rho)

    print "effective potential..."
    potential = effective_potential_func(atomlist_nuclei, rho, xc)

    v = potential
    # solve S.u = v
    nmax = 5
    norm_rem, u = solve_spherical_potential(atomlist, v, nmax=nmax)
    # reconstruct v from u
    #   v_rec = S.u
    v_rec, dummy = spherical_remainder(atomlist, u)

    # v_missing is the difference between the original potential
    # and the reconstructed potential, its norm is the reconstruction
    # error
    #   v_missing = v - v_rec
    #
    # In order to visualize the deviation from the original
    # potential, v_missing is saved to a cube file.
    def v_missing(grid, dV):
        (x,y,z) = grid
        return v(x,y,z) - v_rec(x,y,z)

    from DFTB.Analyse import Cube
    cube_file = "/tmp/water+dummy_lda_vmissing.cube"
    Cube.function_to_cubefile(atomlist, v_missing, filename=cube_file,
                              ppb=5.0)
    print "deviation V - V_rec saved to cube file '%s'" % cube_file

    
    # save tables with cuts along the y-axis of V, V^(rec) and U
    r = np.linspace(-5.0, +5.0, 10000)
    dat_file = "/tmp/water+dummy_lda_reconstructed_potential.dat"
    fh = open(dat_file, "w")
    # z-positions of cuts along y-axis, x=0 in all cuts
    for z in [-0.222540557483415, 0.0, 0.886071388908105]:
        data = np.vstack(( r,
                           v(0*r,r,0*r+z),
                           v_rec(0*r,r,0*r+z),
                           u(0*r,r,0*r+z))).transpose()
        print>>fh, "# cut through potential along the y-axis at x=0, z=%e" % z
        print>>fh, "# "
        print>>fh, "# R/bohr    V/Hartree       V^(rec)/Hartree      U/Hartree"
        np.savetxt(fh, data, fmt="%+10.8e")
        print>>fh, " "
    print "tables with V, V^(rec)=S.U and U=S^(-1).V  written to '%s'" % dat_file
    fh.close()

    
def test_solver_grid_dependence_hmi():
    """
    How does the reconstruction error

         |V - V^{(rec)}| = |V - S (S^{-1} V)|

    depend on the size of the radial and angular grids?
    """
    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    # center of charge
    center = (0.0, 0.0, 0.0)

    # list of real atoms
    atomlist = [(1, posH1),
                (1, posH2)]
    charge = +2

    # choose radial resolution of multicenter grid
    rfac = 10                             # increase density of radial points by a factor of 3
    settings.radial_grid_factor = rfac    # controls size of radial grid
    
    # How does the angular resolution of the grids affect the reconstruction error?
    # The reconstruction error is computed for different grid size (rfac, Lmax).
    from DFTB.MolecularIntegrals.lebedev_quad_points import Lebedev_Lmax, Lebedev_Npts
    # error(Lmax)
    errors = np.zeros(len(Lebedev_Lmax))
    for i,Lmax in enumerate(Lebedev_Lmax):
        settings.lebedev_order = Lmax           # controls size of angular grid

        print_grid_summary(atomlist,
                           settings.lebedev_order, settings.radial_grid_factor)

        
        rho = None
        xc = None
        print "effective potential..."
        potential = effective_potential_func(atomlist, rho, xc, nelec=0)

        # solve S.u = v
        v = potential
        # number of times the iteration
        #    v_{n} = (1 - S) v_{n-1}
        # is performed.
        nmax = 5
        norm_rem, u = solve_spherical_potential(atomlist, v, nmax=nmax)
        # reconstruct v from u
        #   v_rec = S.u
        v_rec, dummy = spherical_remainder(atomlist, u)

        # reconstruction error = |V - V^(rec)| = |V - S.U| = |V - S.S^{-1}.U|
        error = np.sqrt(
            integral(atomlist, lambda x,y,z: (v(x,y,z) - v_rec(x,y,z))**2 ) )

        print "reconstruction error:  rfac= %d  Lmax= %d  error= %e" % (rfac, Lmax, error)

def test_solver_grid_dependence_water():
    """
    How does the reconstruction error

         |V - V^{(rec)}| = |V - S (S^{-1} V)|

    depend on the size of the radial and angular grids?
    """

    # experimental geometry of water
    #  r(OH) = 0.958 Ang, angle(H-O-H) = 104.4776 degrees
    atomlist = [
        (8, (0.000000000000000,  0.000000000000000, -0.222540557483415)),
        (1, (0.000000000000000, +1.431214118579765,  0.886071388908105)),
        (1, (0.000000000000000, -1.431214118579765,  0.886071388908105))]
    charge = 0    
    # LDA functional
    xc = XCFunctionals.libXCFunctional('lda_x', 'lda_c_vwn_3')

    # choose radial resolution of multicenter grid
    rfac = 10                             # increase density of radial points by a factor of 3
    settings.radial_grid_factor = rfac    # controls size of radial grid
    
    # How does the angular resolution of the grids affect the reconstruction error?
    # The reconstruction error is computed for different grid size (rfac, Lmax).
    from DFTB.MolecularIntegrals.lebedev_quad_points import Lebedev_Lmax, Lebedev_Npts
    # error(Lmax)
    errors = np.zeros(len(Lebedev_Lmax))
    for i,Lmax in enumerate(Lebedev_Lmax):
        settings.lebedev_order = Lmax           # controls size of angular grid

        print_grid_summary(atomlist,
                           settings.lebedev_order, settings.radial_grid_factor)

        
        RDFT = BasissetFreeDFT(atomlist, xc, charge=charge)
        # doubly occupied DFTB orbitals
        orbitals = RDFT.getOrbitalGuess()

        # total electron density
        print "electron density..."
        rho = density_func(orbitals)

        qnuc, qelec, qtot = total_charge(atomlist, rho)
        
        print "effective potential..."
        potential = effective_potential_func(atomlist, rho, xc)

        # solve S.u = v
        v = potential
        # number of times the iteration
        #    v_{n} = (1 - S) v_{n-1}
        # is performed.
        nmax = 5
        norm_rem, u = solve_spherical_potential(atomlist, v, nmax=nmax)
        # reconstruct v from u
        #   v_rec = S.u
        v_rec, dummy = spherical_remainder(atomlist, u)

        # reconstruction error = |V - V^(rec)| = |V - S.U| = |V - S.S^{-1}.U|
        error = np.sqrt(
            integral(atomlist, lambda x,y,z: (v(x,y,z) - v_rec(x,y,z))**2 ) )

        print "reconstruction error:  rfac= %d  Lmax= %d  error= %e" % (rfac, Lmax, error)

        
def test_iterative_inhomogeneous(): # NOT WORKING, THERE IS A CONCEPTUAL MISTAKE 
    """
    The molecular potential V is decomposed into a Coulombic potential V0, 
    which has the correct asymptotic shape plus "spherically averaged" corrections
                     oo    (sph)
        V = V  +  sum     V
             0       i=1   i

    We define the partial sums
                      n    (sph)                        n -> oo
        V  = V  +  sum    V              such that  V  ---------->  V
         n    0       i=1  i                         n

    V_0 is a long-range Coulomb potential, while all V_i^(sph) are short-range.

    The Schroedinger equation for a continuum orbital with energy E > 0 is solved
    iteratively, first for V_0, then for V_1, V_2, etc. In each step the solution of
    the previous potential, phi_{i-1}, is used to construct the source term for
    the inhomogeneous Schroedinger equation that yields the next solution.

    
    Step 0:  The solution phi_0 having the correct energy E and asymptotic angular
    momentum (l,m) for the Coulombic potential V_0 is known

         (T + V  - E) phi  = 0              V  = -charge/|r-center|       (0)
               0         0                   0

    Step 1: Starting from the previous solutino phi_0 we search for the solution of

                                                       (sph)
         (T + V  - E) phi  = 0              V  = V  + V                   (1a)
               1         1                   1    0    1

    by making the ansatz that the new solution corresponds to the old one plus a correction

         phi  = phi  + dphi
            1      0       1                                              (1b)

    Substituting this ansatz into the Schroedinger equation for phi_1 we get
                                                      (sph)
         (T + V  - E) dphi  = -(T + V  - E) phi   -  V      phi
               1          1          0         0      1        0

    which, because of eqn. (0) simplifies to an inhomogeneous Schroedinger equation for
    the correction dphi_1
                                 (sph)
         (T + V  - E) dphi  = - V      phi                                (1c)
               1          1      1        0

    Step 2: The same procedure is applied to solve
                                                        (sph)          (sph)   (sph)
         (T + V  - E) phi  = 0               V  = V  + V      =  V  + V     + V
               2         2                    2    1    2         0    1       2

    with the ansatz

          phi  = phi  + dphi
             2      1       2

    which using (1c) leads to
                                 (sph)
         (T + V  - E) dphi  = - V      phi
               2          2      2        1

    and so on ...

    Step n:
                                 (sph)             n-1
         (T + V  - E) dphi  = - V      ( phi  + sum    dphi )
               n          n      n          0      i=1     i

    PROBLEM: I think this won't work, because on the left hand side V_n needs
             to be replaced by V_n^(sph)
    """
    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    # center of charge
    center = (0.0, 0.0, 0.0)


    # list of real atoms
    atomlist_nuclei = [(1, posH1),
                       (1, posH2)]
    charge = +2

    # Because v0 has a singularity at the origin, we have to add
    # a grid around the origin by placing a dummy atom there.
    atomlist = [(1, center),
                (1, posH1),
                (1, posH2)]
    
    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 10      # controls size of radial grid  
    settings.lebedev_order = 25          # controls size of angular grid
    
    rho = None
    xc = None
    print "effective potential..."
    potential_full = effective_potential_func(atomlist_nuclei, rho, xc, nelec=0)

    # energy of continuum orbital (in Hartree)
    E = 0.2 
    # angular momentum quantum numbers
    l,m = 0,0
    
    # The initial guess for the continuum orbital is a Coulomb
    # wave centered at the origin
    phi0 = regular_coulomb_func(E, charge, l, m, 0.0,
                               center=center)
    residual0 = residual_func(atomlist, phi0, potential_full, E)
    
    # phi0 is a solution of
    #  (T + V0 - E) phi = 0
    
    # origin of Coulomb waves and V0
    x0,y0,z0 = center
    def v0(x,y,z):
        # shift coordinates
        x,y,z = x-x0,y-y0,z-z0
        # distance from origin
        r = np.sqrt(x*x+y*y+z*z)
        return -charge/r

    # short-range potential     V - V0
    v_short = add_two_functions(atomlist, potential_full, v0, 1.0, -1.0)

    potential = v0
    phi = phi0
    remainder = v_short

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,3)
    r = np.linspace(-5.0, +5.0, 10000)

    axes[0].set_ylabel("potential / Hartree")
    axes[1].set_ylabel("wavefunction")
    axes[2].set_ylabel("residual")

    axes[0].plot(r, potential_full(0.0*r, 0.0*r, r), label=r"$V(0,0,z)$")
    
    max_iter = 3
    for k in range(1, max_iter):
        print "Iteration k= %d" % k
        spherical, remainder = spherical_remainder(atomlist, remainder)
        potential = add_two_functions(atomlist, potential, spherical, 1.0, 1.0)

        source = lambda x,y,z: -spherical(x,y,z) * phi(x,y,z)
        dphi = inhomogeneous_schroedinger(atomlist, potential_full, source , E)
        phi = add_two_functions(atomlist, phi, dphi, 1.0, 1.0)

        residual = residual_func(atomlist, phi, potential_full, E)
            
        axes[0].plot(r, potential(0.0*r, 0.0*r, r), label=r"$V_{%d}(0,0,z)$" % k)
        axes[1].plot(r, phi(0.0*r,0.0*r,r), label=r"$\phi_{%d}(0,0,z)$" % k)
        axes[2].plot(r, residual(0.0*r,0.0*r,r), ls="-.", label=r"residual $(H-E)\phi_{%d}(0,0,z)$" % k)

    for ax in axes:
        ax.set_xlabel("z / bohr")
        ax.legend()

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.5)
    
    plt.show()

    
def test1_hmi_scattering():
    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    # center of charge
    center = (0.0, 0.0, 0.0)

    
    atomlist = [(1, posH1),
                (1, posH2)]
    charge = +2
    
    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 40      # controls size of radial grid  
    settings.lebedev_order = 31          # controls size of angular grid

    # energy of continuum orbital (in Hartree)
    E = 0.2 
    # angular momentum quantum numbers
    l,m = 0,0

    """
    # The initial guess for the continuum orbital is the sum of
    # two hydrogen continuum orbitals centered on each proton
    cf1 = regular_coulomb_func(E, charge, l, m, 0.0,
                               center=posH1)
    cf2 = regular_coulomb_func(E, charge, l, m, 0.0,
                               center=posH1)

    print "initial guess for continuum orbital..."
    phi0 = add_two_functions(atomlist, cf1, cf2, 1.0/np.sqrt(2), 1.0/np.sqrt(2))
    """
    # The initial guess for the continuum orbital is a Coulomb
    # wave centered at the origin
    phi0 = regular_coulomb_func(E, charge, l, m, 0.0,
                               center=center)

    rho = None
    xc = None
    print "effective potential..."
    potential = effective_potential_func(atomlist, rho, xc, nelec=0)
    
    # origin of Coulomb waves and V0
    x0,y0,z0 = center
    def v0(x,y,z):
        # shift coordinates
        x,y,z = x-x0,y-y0,z-z0
        # distance from origin
        r = np.sqrt(x*x+y*y+z*z)
        return -charge/r

    def v1(x,y,z):
        return potential(x,y,z) - v0(x,y,z)

    # right-hand side of inhomogeneous Schroedinger equation
    def source(x,y,z):
        return -v1(x,y,z) * phi0(x,y,z)

    #
    #  solve  (H0 + V1 - E) dphi = - V1 phi0
    # for orbital correction dphi

    # Because v0 has a singularity at the origin, we have to add
    # a grid around the origin by placing a dummy atom there.
    atomlist_grid = [(3, center),
                     (3, posH1),
                     (3, posH2)]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist_grid)

    print "Schroedinger equation ..."
    dphi = multicenter_inhomogeneous_schroedinger(potential, source,  E,
                                                  atomic_coordinates, atomic_numbers,
                                                  radial_grid_factor=settings.radial_grid_factor,
                                                  lebedev_order=settings.lebedev_order)

    
    
    # Combine asymptotically correct solution with correction 
    #   phi = phi0 + dphi
    phi1 = add_two_functions(atomlist_grid, phi0, dphi, 1.0, 1.0)
    # normalize continuum orbital
    delta, phi1 = phaseshift(atomlist_grid, phi1, E, charge, l, m,
                            center=center)
    
    # residual (H-E)phi0
    print "residual (H-E)phi0..."
    residual0 = residual_func(atomlist_grid, phi0, potential, E)

    # residual (H-E)phi1
    print "residual (H-E)phi1..."
    residual1 = residual_func(atomlist_grid, phi1, potential, E)

    # improve orbital
    delta_e = 0.0
    source = source_func(delta_e, phi1, residual1)
    print "orbital correction..."
    dphi = orbital_correction(atomlist_grid, potential, source,  E, delta_e)
    print "residual (H-E)dphi..."
    residual_dphi = residual_func(atomlist_grid, dphi, potential, E)

    #
    print "build matrix (H-E)^2 in basis {phi,dphi} ..."
    R2 = np.zeros((2,2))
    R2[0,0] = overlap(atomlist_grid, residual1, residual1)                                   #  <phi|(H-E)^2|phi>
    R2[0,1] = overlap(atomlist, residual1, residual_dphi)      #  <phi|(H-E)^2|delta_phi>
    R2[1,0] = R2[0,1]                                         #  all wavefunctions are real
    R2[1,1] = overlap(atomlist, residual_dphi, residual_dphi) #  <delta_phi|(H-E)^2|delta_phi>

    # diagonalize R^2
    eigvals, eigvecs = sla.eigh(R2)
    
    print "matrix R^2"
    labels = ["phi", "delta phi"]
    print annotated_matrix(R2, labels, labels)
    
    # and choose the solution with the lowest 
    # eigenvalue, I think all eigenvalues should be positive
    
    print "Eigenvalues: %s" % eigvals
    a,b = eigvecs[:,0]
    
    # Sign of an eigenvector is arbitrary. We choose
    # the sign such that `a` is always positive.
    if a < 0:
        a *= -1
        b *= -1
        
    print "Coefficients for mixing old orbital and orbital correction"
    print "   phi^(new) = a*phi^(old) + b*delta_phi "
    print "a = %s" % a
    print "b = %s" % b

    phi2 = add_two_functions(atomlist_grid, phi1, dphi, a,b)
    # normalize continuum orbital
    delta, phi2 = phaseshift(atomlist_grid, phi2, E, charge, l, m,
                            center=center)
    print "residual (H-E)phi2 ..."
    residual2 = residual_func(atomlist_grid, phi2, potential, E)

    
    # plot cuts through effective potential and residual
    import matplotlib.pyplot as plt
    
    r_small = np.linspace(-5.0, 5.0, 10000)
    r_large = np.linspace(1.0, 300.0, 20000)

    rs = [r_small, r_large]
    fig, axes = plt.subplots(1,2)

    for ax,r in zip(axes, rs):
        ax.set_xlabel("z / bohr")

        ax.plot(r, potential(0*r,0*r,r), lw=2, label=r"$v_{eff}(0,0,z)$")
        ax.plot(r, v0(0*r,0*r,r), lw=2, label=r"$v_0(0,0,z)$")
        ax.plot(r, v1(0*r,0*r,r), lw=2, label=r"$v_1(0,0,z)$")

        ax.plot(r, phi0(0*r,0*r,r), ls="-.", label=r"$\phi_0(0,0,z)$")
        ax.plot(r, phi1(0*r,0*r,r), ls="-.", label=r"$\phi_1(0,0,z)$")
        ax.plot(r, phi2(0*r,0*r,r), ls="-.", label=r"$\phi_2(0,0,z)$")

        # residuals
        ax.plot(r, residual0(0*r,0*r,r), label=r"residual $R(z) = (H-E)\phi_0$")
        ax.plot(r, residual1(0*r,0*r,r), label=r"residual $R(z) = (H-E)\phi_1$")
        ax.plot(r, residual2(0*r,0*r,r), label=r"residual $R(z) = (H-E)\phi_2$")
        
        ax.legend()
        
    plt.show()
    exit(0)
    #########################
    
    phi = phi0
    residual = residual0
    
    max_iter = 10
    for k in range(0, max_iter):
        delta_e = energy_correction(atomlist, residual, phi, method="Becke")
        print "delta_e = %e" % delta_e
        delta_e = 0.0
        source = source_func(delta_e, phi, residual)
        delta_phi = orbital_correction(atomlist, potential, source, E, delta_e)

        phi = add_two_functions(atomlist, phi, delta_phi, 1.0/np.sqrt(2), 1.0/np.sqrt(2))
        residual = residual_func(atomlist, phi, potential, E)
        nrm2_residual = overlap(atomlist, residual, residual)
        print "total error  |(H-E)phi|^2= %e" % nrm2_residual
        
        
        plt.plot(r, potential(0*r,0*r,r), label=r"$v_{eff}(0,0,z)$")
        plt.plot(r, residual0(0*r,0*r,r), label=r"residual $R(z) = (H-E)\phi_0$")
        plt.plot(r, residual(0*r,0*r,r), label=r"residual $R(z) = (H-E)\phi$")
        
        plt.ylim((-2.0, 2.0))
        plt.legend()
        plt.show()

def test2_hmi_scattering():
    import warnings
    warnings.filterwarnings("error")

    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    
    atomlist = [(1, posH1),
                (1, posH2)]
    charge = +2
    
    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 40 #60     # controls size of radial grid  
    settings.lebedev_order = 21 #41 #65          # controls size of angular grid

    RDFT = BasissetFreeDFT(atomlist, None, charge=charge)
    
    # energy of continuum orbital (in Hartree)
    E = 0.2 
    # angular momentum quantum numbers
    l,m = 0,0

    rho = None
    homo = None
    delta, phi = RDFT.solveScatteringProblem(rho, homo, E, l, m, nmax=10)
    #delta, phi = RDFT.solveScatteringProblem_new(rho, E, l,m, nmax=5)

def test_imaginary_time_propagation():
    import warnings
    warnings.filterwarnings("error")
    
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 60      # controls size of radial grid 
    settings.lebedev_order = 3          # controls size of angular grid
    
    # hydrogen atom
    atomlist = [(1, (0.0, 0.0, 0.0))]
    Z = 1

    # show number of radial and angular points in multicenter grid
    print_grid_summary(atomlist,
                       settings.lebedev_order, settings.radial_grid_factor)
    
    # The continuum orbital is specified by its energy and asymptotic
    # angular momentum (E,l,m)
    
    # energy of continuum orbital
    E = 1.0  
    # length of wave vectors
    k = np.sqrt(2*E)
    # angular quantum numbers of asymptotic solution  (px-orbital)
    l = 0
    m = 0

    # There are no bound orbitals, so the density of H+ is rho=0
    rho = density_func([])

    # soft-core Coulomb potential
    #   V(r) = -1/sqrt(r^2 + eps^2)
    eps = 1.0e-2
    def potential(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return -Z/np.sqrt(r**2+eps**2)

    # initial guess
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    print "continuum Schroedinger equation..."
    phi0 = multicenter_continuum_schroedinger(potential, E, +1, 
                                             atomic_coordinates, atomic_numbers,
                                             radial_grid_factor=settings.radial_grid_factor,
                                             lebedev_order=settings.lebedev_order)
    ####
        
    import matplotlib.pyplot as plt
    r = np.linspace(-0.01, 0.01, 50000)

    x = 0*r
    y = 0*r
    z = r

    # total propagation time
    time = 0.01
    # number of time steps for propagation from 0 to t
    nsteps = 1

    lap0 = laplacian_func(atomlist, phi0)
    
    plt.plot(r, potential(x,y,z), label="-1/sqrt(r^2+eps^2)")
    plt.plot(r, phi0(x,y,z), label="$\phi_0$")
    plt.plot(r, -0.5 * lap0(x,y,z), ls="-.", label=r"kinetic energy $T \phi_0$")
    plt.plot(r, (potential(x,y,z)-E)*phi0(x,y,z), ls="-.", label=r"potential energy $(V-E)\phi_0$")
    plt.plot(r, -0.5 * lap0(x,y,z) + (potential(x,y,z)-E)*phi0(x,y,z), ls="--", label="residual $(T+V-E)\phi_0$")
    
    plt.legend()
    plt.show()
    
    print "residual of initial guess..."
    residual0 = residual_func(atomlist, phi0, potential, E)

    plt.plot(r, residual0(x,y,z), ls="--", label=r"residual $(H-E)\phi_0$")

    lap_res = laplacian_func(atomlist, residual0)
    plt.plot(r, lap_res(x,y,z), ls="-.", label=r"$-1/2 \nabla^2$ residual")

    plt.legend()
    plt.show()
    
    phi = phi0
    for i in range(0, 3):
        print "imaginary time propagation  exp(-(H-E)^2 t) phi0 ..."
        phi = imaginary_time_propagation(atomlist, phi, potential, E, time, nsteps)

        plt.plot(r, phi(x,y,z), ls="-.", label="$\phi_{%d}$" % (i+1))

    plt.legend()
    plt.show()


def test_xc_functional():
    """
    compare exchange-correlation energy computed on the 
    multicenter grid with the values from atomic calculations
    using the same radial density.
    """
    # load nitrogen atom
    from DFTB.SlaterKoster.free_pseudo_atoms import n as atom
    
    from scipy import interpolate

    atomlist = [(atom.Z, (0.0, 0.0, 0.0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # same functional is used as for the pseudo atom calculation
    xc = XCFunctionals.libXCFunctional(atom.pseudo_orbital_x, atom.pseudo_orbital_c)

    # interpolate radial density
    rho_rad = interpolate.interp1d(atom.r[atom.r > 0.0], atom.radial_density[atom.r > 0.0],
                                   kind='cubic', fill_value=(atom.radial_density[1], 0.0), bounds_error=False)

    def rho(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        n = rho_rad(r)
        return n

    #          __     2
    # sigma = (\/ rho)
    sigma = multicenter_gradient2(rho,
                                  atomic_coordinates, atomic_numbers,
                                  radial_grid_factor=settings.radial_grid_factor,
                                  lebedev_order=settings.lebedev_order)

    import matplotlib.pyplot as plt
    r = np.linspace(-10.0, 10.0, 1000)
    x = 0*r
    y = 0*r
    z = r
    plt.plot(r, rho(x,y,z), label=r"$\rho$")
    #plt.plot(r, sigma(x,y,z), label=r"$\vert \nabla \rho \vert^2$")
    # exchange and correlation potentials
    plt.plot(r, xc.func_x.vxc(rho(x,y,z), sigma(x,y,z)), ls="--", label=r"$V_x[\rho]$")
    plt.plot(r, xc.func_c.vxc(rho(x,y,z), sigma(x,y,z)), ls="--", label=r"$V_c[\rho]$")
    # exchange and correlation energy densities
    plt.plot(r, xc.func_x.exc(rho(x,y,z), sigma(x,y,z)), ls="-.", label=r"$e_x[\rho]$")
    plt.plot(r, xc.func_c.exc(rho(x,y,z), sigma(x,y,z)), ls="-.", label=r"$e_c[\rho]$")
    
    plt.legend()
    
    plt.show()
    

    # exchange-correlation energy
    Exc = exchange_correlation_energy(atomlist, rho, xc)

    print "exchange-correlation energy (multicenter grid)   = %e" % Exc
    print "exchange-correlction energy (atomic calculation) = %e" % atom.Exc

    # classical Coulomb energy
    Ecoul = coulomb_energy(atomlist, rho)
    print "Coulomb energy (multicenter grid)   = %e" % Ecoul
    print "Coulomb energy (atomic calculation) = %e" % atom.Ecoul

    
if __name__ == "__main__":
    #test_hydrogen_molecular_ion()
    #test_h3_ion()
    #test_lithium_cation()
    #test_h2()
    #test_lithium_dimer() # not working yet
    #test_h2o()  # not working yet

    #test_coulomb_waves()
    #test_scattering_solution()
    #test_phaseshift_overlap()
    #test_asymptotic_regular_coulomb()
    #test_improved_asymptotic_regular_coulomb()
    
    #test1_lithium_scattering_solution()   # not working yet
    #test2_lithium_scattering_solution()   # not working yet

    #test1_hmi_scattering() # not working yet
    #test2_hmi_scattering() # not working yet

    #test1_spherical_remainder()
    #test2_spherical_remainder()

    #test_iterative_inhomogeneous()  # doesn't work

    #test_solve_spherical_potential_hmi()
    #test_solve_spherical_potential_water()
    #test_solve_spherical_potential_water_dummy()
    
    #test_solver_grid_dependence_hmi()
    #test_solver_grid_dependence_water()

    #test_imaginary_time_propagation()

    test_xc_functional()
