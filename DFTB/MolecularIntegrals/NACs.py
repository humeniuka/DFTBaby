#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compute non-adiabatic coupling (NAC) vectors between two electronic states 
numerically using Becke's multicenter integration scheme.

Via a sum-over-states expansion the non-adiabatic coupling vector can be rewritten as
                             <I|(d/dR Vne)|J>
   <state I| d/dR state J> = -----------------                     (1)
                                E_J - E_I

where 

         elecs.  atoms      -Za
   Vne =  sum    sum    -----------
          i       A     |R_A - r_i|

is the attraction between electrons and nuclei. By introducing the transition density
rho_IJ(r), the A-component of the derivative vector in the numerator of eqn. (1)
may be expressed as the following integral:

                                     /   R_A - r
   nac_A =  <I|(d/dR_A Vne)|J> = Z_A |  ----------- rho_IJ(r) dr 
                                     /  |R_A - r|^3

such that 

    <I|d/dR_A J> = nac_A / (E_J - E_I)

In some basis of atomic orbitals {phi_a(r)}, the transition density can be presented by the transition density matrix
 
               AOs    I,J         *
   rho_IJ(r) = sum   P      phi(r)  phi(r) 
               a,b    a,b      a       b

so that nac_A becomes

               AOs   I,J  /           *    R_A - r
   nac_A = Z_A sum  P     | dr  phi(r)   -----------  phi(r)
               a,b   a,b  /        a     |R_A - r|^3     b

"""
from DFTB.BasisSets import AtomicBasisSet
from DFTB import XYZ, AtomicData
from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_integration

import numpy as np

def non_adiabatic_coupling_integral_rho(atomlist, rhoAB):
    """
    compute the non-adiabatic coupling vector for a transition density

       a,b       /         R_A - r
    nac    = Z_A | dr    ----------- rho (r) 
       A         /       |R_A - r|^3    ab

    by numerical intergration

    Parameters
    ==========
    atomlist            :  list of tuples (Z,[x,y,z]) with atomic numbers
                           and positions
    rhoAB               :  callable, rhoAB(x,y,z) computes the electronic transition
                           density

    Returns
    =======
    nac                 :  numpy array of shape (3,Nat), non-adiabatic coupling vector
    """
    # bring data into a form understood by the module MolecularIntegrals
    Nat = len(atomlist)
    atomic_numbers = np.zeros(Nat, dtype=int)
    atomic_coordinates = np.zeros((3,Nat))
    for i in range(0, Nat):
        Z,pos = atomlist[i]
        atomic_numbers[i] = Z
        atomic_coordinates[:,i] = pos
    # Now we compte the integrals numerically on a multicenter grid.
    # The mesh size is controlled by the following parameters
    radial_grid_factor = 3      # controls size of radial grid 
    lebedev_order = 23          # controls size of angular grid

    # compute coupling integrals for each atom
    nac = np.zeros((3,Nat))
    for i in range(0, Nat):
        Zat,pos = atomlist[i]
        # Solve integrals for each cartesian component
        for xyz in [0,1,2]:  # 

            # print progress
            xyz2str = ["X","Y","Z"]
            print "  %d of %d    atom %s-%d  %s component" % (3*i+xyz+1, 3*Nat,
                                                                 AtomicData.atom_names[Zat-1], i, xyz2str[xyz])

            # define integrand
            def Inac_xyz_integrand(x,y,z):
                # |R_A - r|
                dist = np.sqrt((pos[0]-x)**2 + (pos[1]-y)**2 + (pos[2]-z)**2)
                r = [x,y,z]
                #      R_A - r
                # Za -----------
                #    |R_A - r|^3
                
                # Since only valence electrons are considered, we maybe need to use the 
                # charges of the core (atomic number - nr. cor electrons) instead
                # of the full nuclear charges? But I am not sure.
                Zcore = AtomicData.valence_electrons[AtomicData.atom_names[Zat-1]]
                #Inac = Zat * (pos[xyz]-r[xyz])/dist**3 * rhoAB(x,y,z)
                Inac = Zcore * (pos[xyz]-r[xyz])/dist**3 * rhoAB(x,y,z)
                
                return Inac

            # perform numerical integration
            nac[xyz,i] = multicenter_integration(Inac_xyz_integrand, atomic_coordinates, atomic_numbers,
                                                   radial_grid_factor=radial_grid_factor,
                                                   lebedev_order=lebedev_order)

    return nac
    

def non_adiabatic_coupling_vector(atomlist, basis, Ptrans, en_diff):
    """
    compute non-adiabatic coupling vector from transition density matrix P in
    the AO basis

    Parameters
    ==========
    atomlist      : molecular geometry, list of Nat tuples (Z,[x,y,z])
    basis         : instance of AtomicBasisSet with M basis functions
    Ptrans        : numpy array of shape (Mao,Mao) with transition density matrix
                    Ptras_IJ
    en_diff       : energy difference between the two electronic states, en(J) - en(I)

    Returns
    =======
    nac           : numpy array with shape (3,Nat), non-adiabatic coupling vector
                    nac[:,i] is the vector belonging to atom i
    """
    
    # define function for evaluating transition density on a grid
    def transition_density_func(x,y,z):
        
        # evaluate atomic orbitals on the grid (x,y,z)
        nbfs = len(basis.bfs)
        bfs_xyz = []
        for a in range(0, nbfs):
            # evaluate orbital `a` on the grid
            bfA_xyz = basis.bfs[a].amp(x,y,z)
            bfs_xyz.append(bfA_xyz)
            
        # evaluate transition density on the grid
        rho_trans = 0.0*x
        # rho = sum_{a,b} P_ab orb_a(r) * orb_b(r)
        for a in range(0, nbfs):
            for b in range(a, nbfs):
                # Here we assume that the orbitals are reals
                rho_ab = bfs_xyz[a] * bfs_xyz[b]
                rho_trans += Ptrans[a,b] * rho_ab
                if (b > a):
                    rho_trans += Ptrans[b,a] * rho_ab

        return rho_trans

    # non-adiabatic coupling vector from transition density
    nac = non_adiabatic_coupling_integral_rho(atomlist, transition_density_func)
    # divide by energy difference (E_J-E_I) in denomintaor
    nac /= en_diff
        
    return nac


if __name__ == "__main__":
    import sys
    import os.path
    
    usage = """
       Usage: %s    <.xyz file>   <transition density matrix>  <.nac output file>

    compute the non-adiabatic coupling vector 

                                < S0 | ( d/dR V_{nuc-elec} ) | S1 >
            < S0 | d/dR S1 > = ------------------------------------
                                         E(1) - E(0)

    from the transition density matrix rho_01(r) in the AO basis by numerical integration.
    The energy difference in the denominator is not included (you have to divide the result
    by the excitation energy yourself). 

    Input files:
      - file with molecular geometry in XYZ format
      - file with transition density matrix in AO basis, matrix of shape 
        (num.orbs x num.orbs)
    Output file:
      - file with coupling vector as a matrix of shape (num.atoms,3)

    """ % os.path.basename(sys.argv[0])
    
    if len(sys.argv) < 4:
        print usage
        exit(-1)
    args = sys.argv[1:]

    # input files
    xyz_file     = args[0]
    density_file = args[1]
    # output file
    nac_file     = args[2]
    
    
    # load molecular geometry 
    atomlist = XYZ.read_xyz(xyz_file)[0]
    # and transition density
    Ptrans = np.loadtxt(density_file)
    # create basis of numerical atomic orbitals
    basis = AtomicBasisSet(atomlist)
    # compute NAC vector
    # Set energy difference to 1.0, you have to change this to the actual
    # energy gap between S1 and S0
    en_diff = 1.0
    nac = non_adiabatic_coupling_vector(atomlist, basis, Ptrans, en_diff)

    # and save it to a file
    fh_nac = open(nac_file, "w")
    # write header for NAC vector
    print>>fh_nac, "# non-adiabatic coupling vector ( in bohr^-1 ) "
    print>>fh_nac, "#    X             Y             Z             "
    
    # print non-adiabatic coupling vector as (Nat,3) matrix
    Nat = len(atomlist)
    np.savetxt(fh_nac, nac.transpose(), fmt="%+e")

    print "Non-adiabatic coupling vector written to %s" % nac_file
