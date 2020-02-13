#!/usr/bin/env python
"""
generate the DFT orbitals for first, second and third row atoms in a confining potential
and save them to a file. Also calculate energies of partially charged atoms to estimate 
curvature of energy E(dq) as a function of surplus charge dq, from which the hardness chi 
and the Hubbard parameter U can be derived.

Note:
=====
If we start with an initial density of zero everywhere the eigen energies
after the first iteration can be considerably lower than the final values. Therefore
the energy range has to extend to much lower energies than the expected lowest
orbital energy.
"""
from PseudoAtomDFT import PseudoAtomDFT, occupation_numbers
from DFTB.AtomicData import atom_names
from numpy import linspace, array, inf
import os.path

script_dir = os.path.dirname(os.path.realpath(__file__))
orbdir = os.path.join(script_dir, "free_pseudo_atoms/")

Npts = 3000 # number of radial grid points
rmin = 0.0
rmax = 15.0

# maximal deviation at matching point for Numerov method
numerov_conv = 1.0e-6
# threshold for convergence of SCF calculation
en_conv = 1.0e-5

######################### NEUTRAL ATOMS ######################################

# same code is repeated for all atoms, this functions bundles the repetitions
def run_atomdft(Z=1,Nelec=1, energy_range=[-2,-1,-0.001]):
    # uses global variables orbdir, numerov_conv, en_conv, rmin, rmax, Npts
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import pseudo_atoms_list
        at = pseudo_atoms_list[Z-1]
        atomdft.initialDensityGuess((at.r, at.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    except KeyError:
        # use density of other atom with lower Z as initial guess
        for Z_o in range(Z-1,0-1,-1):
            if pseudo_atoms_list.has_key(Z_o-1):
                at = pseudo_atoms_list[Z_o-1]
                atomdft.initialDensityGuess((at.r, at.radial_density))
                break
        atomdft.initialDensityGuess()        
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "%s.py" % atom_names[Z-1], "w")
    atomdft.saveSolution(fh)
    fh.close()

def hydrogen():
    dE = 0.1
    energy_range = linspace(-2.0, -0.001, 100)

    run_atomdft(Z=1,Nelec=1, energy_range=energy_range)
    """
    Z = 1
    Nelec = 1

    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import h
        atomdft.initialDensityGuess((h.r, h.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "h.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    """

def helium():
    dE = 0.1
    eguess_1s = -0.56989514
    energy_range = array([eguess_1s - dE, eguess_1s + dE])

    run_atomdft(Z=2,Nelec=2, energy_range=energy_range)
    """
    Z = 2
    Nelec = 2
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import he
        atomdft.initialDensityGuess((he.r, he.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "he.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    """

def lithium():
    energy_range = list(linspace(-2.5, -1.5, 5)) \
        + list(linspace(-0.3, -0.001, 5))

    run_atomdft(Z=3,Nelec=3, energy_range=energy_range)
    """
    Z = 3
    Nelec = 3
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import li
        atomdft.initialDensityGuess((li.r, li.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "li.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    """

def beryllium():
    energy_range = list(linspace(-5.0, -2.0, 5)) \
        + list(linspace(-0.3, -0.1, 5))

    run_atomdft(Z=4,Nelec=4, energy_range=energy_range)
    """
    Z = 4
    Nelec = 4
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import be
        atomdft.initialDensityGuess((be.r, be.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "be.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    """
    
def boron():
    energy_range = list(linspace(-7.0, -6.0, 5)) \
        + list(linspace(-0.5, -0.2, 5)) \
        + list(linspace(-0.2, -0.01, 5))

    run_atomdft(Z=5,Nelec=5, energy_range=energy_range)
    """
    Z = 5
    Nelec = 5
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import b
        atomdft.initialDensityGuess((b.r, b.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "b.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    """

def carbon():
#    energy_range = array( list(linspace(-11, -5.0, 20)) \
#                         +list(linspace(-5.0, -0.3, 20)) \
#                         +list(linspace(-0.3,-0.001, 20)) \
#                         +list(linspace(-0.001, 3.0, 20)))
    energy_range = array( list(linspace(-11.0, -5.0, 20)) \
                         +list(linspace(-5.0, -1.0, 20))  \
                         +list(linspace(-1.0, 1.0, 200)) )

    Z = 6
    Nelec = 6
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
#    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, grid_spacing="exponential",\
#                r0=100.0) 
#    # It seems that a free carbon atom does not have a bound 3d orbital.
#    # If one adds a shallow confinement radius all orbitals become bound
#    atomdft.setRadialGrid(rmin, rmax+8.0, 2*Npts)
#    # add polarization functions => 3d orbitals
#    atomdft.setValenceOrbitals(["2s", "2p", "3d"], format="spectroscopic")
    #
    try:
        from free_pseudo_atoms import c
        atomdft.initialDensityGuess((c.r, c.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "c.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def nitrogen():
#    energy_range = list(linspace(-16.0, -12.0, 5)) \
#        + list(linspace(-0.9, -0.5, 5)) \
#        + list(linspace(-0.4, -0.1, 5))
    energy_range = list(linspace(-16.0, -10.0, 10)) \
        + list(linspace(-0.9, -0.5, 10)) \
        + list(linspace(-0.49, -0.05, 10)) \
        + list(linspace(-0.05, 1.5, 200))

    Z = 7
    Nelec = 7
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4)
    atomdft.setRadialGrid(rmin, rmax, Npts)

#    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4, grid_spacing="exponential",\
#                r0=100.0) 
#    # It seems that a free nitrogren atom does not have a bound 3d orbital.
#    # If one adds a shallow confinement radius all orbitals become bound
#    atomdft.setRadialGrid(rmin, rmax+8.0, 2*Npts)
#    # add polarization functions => 3d orbitals
#    atomdft.setValenceOrbitals(["2s", "2p", "3d"], format="spectroscopic")
#    #
    try:
        from free_pseudo_atoms import n
        atomdft.initialDensityGuess((n.r, n.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "n.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def oxygen():
    energy_range = list(linspace(-30.0, -15.0, 5)) \
        + list(linspace(-2.5, -0.5, 10)) \
        + list(linspace(-0.5, 0.001, 10))

    Z = 8
    Nelec = 8
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, maxiter=1000)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import o
        atomdft.initialDensityGuess((o.r, o.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "o.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    
def fluorine():
    energy_range = list(linspace(-25.0, -23.0, 5)) \
                  +list(linspace(-1.5, -0.8, 5)) \
                  +list(linspace(-0.8, -0.001, 5))

    Z = 9
    Nelec = 9
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import f
        atomdft.initialDensityGuess((f.r, f.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "f.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    
def neon():
    energy_range = list(linspace(-31.0, -29.0, 5)) \
                  +list(linspace(-1.5, -1.0, 5)) \
                  +list(linspace(-1.0, -0.001, 5))

    Z = 10
    Nelec = 10
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import ne
        atomdft.initialDensityGuess((ne.r, ne.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "ne.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

# third row atoms
rmax_3rd_row = rmax + 5.0 # increase radial grid by 3 bohr for 3rd row atoms

def natrium():
    energy_range = list(linspace(-45.0, -15.00, 50)) \
                  +list(linspace(-15, -3.0, 50)) \
                  +list(linspace(-3.0, -1.5, 100)) \
                  +list(linspace(-1.5, -0.0001, 200))

    Z = 11
    Nelec = 11
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax_3rd_row, Npts)
    try:
        from free_pseudo_atoms import na
        atomdft.initialDensityGuess((na.r, na.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "na.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

# Add magnesium as it is contained in bacteriochlorophyll
def magnesium():
    energy_range = list(linspace(-50.0, -40.00, 5)) \
                  +list(linspace(-4, -2.0, 5)) \
                  +list(linspace(-2.0, -1.0, 5)) \
                  +list(linspace(-1.0, -0.0001, 20))

    Z = 12
    Nelec = 12
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax_3rd_row, Npts)
    try:
        from free_pseudo_atoms import mg
        atomdft.initialDensityGuess((mg.r, mg.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "mg.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def aluminum():
    energy_range = list(linspace(-60.0, -40.00, 100)) \
                  +list(linspace(-40, -10.0, 50)) \
                  +list(linspace(-10.0, -1.0, 200)) \
                  +list(linspace(-1.0, -0.0001, 200))

    Z = 13
    Nelec = 13
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax_3rd_row, Npts)
    try:
        from free_pseudo_atoms import al
        atomdft.initialDensityGuess((al.r, al.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "al.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def silicon():
    energy_range = list(linspace(-70.0, -40.00, 100)) \
                  +list(linspace(-40, -10.0, 50)) \
                  +list(linspace(-10.0, -0.00001, 100)) \
                  +list(linspace(-0.001, 0.3, 1000))  # d-orbitals have positive energy

    Z = 14
    Nelec = 14
    # A large confinement radius hopefully does not change the orbital energies but facilitates the calculation.
    # So far I couldn't find the d-orbitals with the shooting method if r0=inf (no confinement).
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4, grid_spacing="exponential", r0=15.0)
    atomdft.setValenceOrbitals(["3s","3p","3d"], format="spectroscopic")
    atomdft.setRadialGrid(rmin, rmax_3rd_row+10.0, Npts+3000)
    try:
        from free_pseudo_atoms import si
        atomdft.initialDensityGuess((si.r, si.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "si.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def phosphorus():
    energy_range = list(linspace(-80.0, -40.00, 100)) \
                  +list(linspace(-40, -10.0, 50)) \
                  +list(linspace(-10.0, -1.0, 200)) \
                  +list(linspace(-1.0, -0.0001, 200))

    Z = 15
    Nelec = 15
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax_3rd_row-5.0, Npts+500)
    try:
        from free_pseudo_atoms import p
        atomdft.initialDensityGuess((p.r, p.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "p.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def sulfur():
    energy_range = list(linspace(-100.0, -40.00, 100)) \
                  +list(linspace(-40, -10.0, 100)) \
                  +list(linspace(-10.0, -1.0, 100)) \
                  +list(linspace(-1.0, -0.000001, 100)) \
                  +list(linspace(0.0001, 0.2, 100)) # unoccupied 3d orbital has positive energy

    Z = 16
    Nelec = 16
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4, grid_spacing="exponential")
    # add unoccupied d orbitals to minimal basis
    atomdft.setValenceOrbitals(["3s","3p","3d"], format="spectroscopic")
    atomdft.setRadialGrid(rmin, rmax_3rd_row, Npts)
    try:
        from free_pseudo_atoms import s
        atomdft.initialDensityGuess((s.r, s.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "s.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def chlorine():
    energy_range = list(linspace(-110.0, -40.00, 100)) \
                  +list(linspace(-40, -13.0, 50)) \
                  +list(linspace(-13.0, -1.0, 200)) \
                  +list(linspace(-1.0, -0.0001, 200))

    Z = 17
    Nelec = 17
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax_3rd_row, Npts)
    try:
        from free_pseudo_atoms import cl
        atomdft.initialDensityGuess((cl.r, cl.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "cl.py", "w")
    atomdft.saveSolution(fh)
    fh.close()


def argon():
    energy_range = list(linspace(-120.0, -40.00, 100)) \
                  +list(linspace(-40, -13.0, 50)) \
                  +list(linspace(-13.0, -1.0, 200)) \
                  +list(linspace(-1.0, -0.0001, 200))

    Z = 18
    Nelec = 18
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=5.0e-4, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax_3rd_row, Npts)
    try:
        from free_pseudo_atoms import ar
        atomdft.initialDensityGuess((ar.r, ar.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "ar.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

# fourth row atoms
rmax_4th_row = rmax + 10.0 # increase radial grid by 10 bohr for 4th row atoms
Npts_4th_row = Npts + 1000 

def potassium():
    energy_range = list(linspace(-135.0, -40.00, 100)) \
                  +list(linspace(-40, -13.0, 50)) \
                  +list(linspace(-13.0, -1.0, 200)) \
                  +list(linspace(-1.0, -0.0001, 200))

    Z = 19
    Nelec = 19
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, grid_spacing="exponential")
    # in potassium 4s is filled before 3d: [Ar] 4s^1
    occupation = [occ for occ in occupation_numbers(Nelec)]
    assert occupation[-1] == (1, 3-1, 2) # 1e in 3d
    occupation[-1] = (1,4-1,0) # 1e in 4s
    atomdft.setOccupation(occupation)
    atomdft.setRadialGrid(rmin, rmax_4th_row, Npts_4th_row)
    try:
        from free_pseudo_atoms import k
        atomdft.initialDensityGuess((k.r, k.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "k.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def scandium():
    energy_range = list(linspace(-400.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 21
    Nelec = 21
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, grid_spacing="exponential") #, damping=0.6)    
    # scandium has the electron configuration: [Ar] 3d^1 4s^2
    occupation = [occ for occ in occupation_numbers(Nelec)]
    assert occupation[-1] == (3,3-1,2) # 3e in 3d
    occupation[-1] = (1,3-1,2)     # put 1e in 3d
    occupation.append( (2,4-1,0) ) # and 2e in 4s
    atomdft.setOccupation(occupation)

    atomdft.setValenceOrbitals(["3d", "4s"], format="spectroscopic")
    atomdft.setRadialGrid(0.0000004, 16.0, 6000)
    try:
        from free_pseudo_atoms import sc
        atomdft.initialDensityGuess((sc.r, sc.radial_density))
        #from free_pseudo_atoms import ti
        # load density of titanium and scale it
        #atomdft.initialDensityGuess((ti.r, ti.radial_density * float(Nelec)/float(ti.Nelec)))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "sc.py", "w")
    atomdft.saveSolution(fh)
    fh.close()
    
def titanium():
    energy_range = list(linspace(-400.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 22
    Nelec = 22
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, grid_spacing="exponential") #, damping=0.6)    
    # titanium has the electron configuration: [Ar] 3d^2 4s^2
    occupation = [occ for occ in occupation_numbers(Nelec)]
    assert occupation[-1] == (4,3-1,2) # 4e in 3d
    occupation[-1] = (2,3-1,2)     # put 2e in 3d
    occupation.append( (2,4-1,0) ) # and 2e in 4s
    atomdft.setOccupation(occupation)

    atomdft.setRadialGrid(0.0000004, 16.0, 6000)
    try:
        from free_pseudo_atoms import ti
        atomdft.initialDensityGuess((ti.r, ti.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "ti.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

    
def iron():
    energy_range = list(linspace(-400.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 26
    Nelec = 26
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=1.0e-4, grid_spacing="exponential", damping=0.6)    
    # iron has the electron configuration: [Ar] 3d^(6) 4s^2
    occupation = [occ for occ in occupation_numbers(Nelec)]
    assert occupation[-1] == (8,3-1,2) # 8e in 3d
    occupation[-1] = (6,3-1,2)     # put 6e in 3d
    occupation.append( (2,4-1,0) ) # and 2e in 4s
    atomdft.setOccupation(occupation)

    atomdft.setRadialGrid(0.0000004, 16.0, 6000)
    try:
        from free_pseudo_atoms import fe
        atomdft.initialDensityGuess((fe.r, fe.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "fe.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def copper():
    energy_range = list(linspace(-400.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 29
    Nelec = 29
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=1.0e-4, grid_spacing="exponential", damping=0.6)    
    atomdft.setRadialGrid(0.0000004, 20.0, 6000)
    try:
        from free_pseudo_atoms import cu
        atomdft.initialDensityGuess((cu.r, cu.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "cu.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def zinc():
    energy_range = list(linspace(-400.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 30
    Nelec = 30
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=1.0e-4, grid_spacing="exponential", damping=0.6)    
    atomdft.setRadialGrid(0.0000004, 16.0, 6000)
    try:
        from free_pseudo_atoms import zn
        atomdft.initialDensityGuess((zn.r, zn.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "zn.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def bromine():
    energy_range = list(linspace(-1500.0, -400.0, 200)) \
        + list(linspace(-400.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 35
    Nelec = 35
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=1.0e-4, grid_spacing="exponential", damping=0.6)
    # 3d-shell is closed (occupied by 10 electrons)
    atomdft.setValenceOrbitals(["4s", "4p"], format="spectroscopic")    
    atomdft.setRadialGrid(0.0000004, 16.0, 6000)
    try:
        from free_pseudo_atoms import br
        atomdft.initialDensityGuess((br.r, br.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "br.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def ruthenium():
    energy_range = list(linspace(-1500.0, -800.0, 100)) \
        + list(linspace(-800.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace( -30.0, -10.0, 100)) \
        + list(linspace( -10.0, -4.0, 300)) \
        + list(linspace(  -4.0, +0.04, 300)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 44
    Nelec = 44
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=1.0e-4, grid_spacing="exponential", damping=0.6)
    # The electron configuration of ruthenium is  [Kr] 4d^(7) 5s^(1)
    occupation = [occ for occ in occupation_numbers(Nelec)]
    # `occupation` is a list of tuples
    #   (nr.electrons, quantum number n + 1, quantum number l)
    # that describe the occupied shells.
    # In the default occupation, ruthenium would have the electronic configuration [Kr] 4d^(8).
    assert occupation[-1] == (8,4-1,2)
    # The default occupation is not correct, so we need to adjust it by
    # removing one electron from 4d shell,
    occupation[-1] = (7,4-1,2)
    # and adding a 5s shell with 1 electron.
    occupation.append( (1,5-1,0) )
    # The new occupation is now [Kr] 4d^(7) 5s^(1)
    
    atomdft.setOccupation(occupation)
    atomdft.setValenceOrbitals(["4d", "5s"], format="spectroscopic")
    atomdft.setRadialGrid(0.0000000001, 18.0, 20000)
    try:
        """
        # If no initial density is available for Ru, we can take the density from silver 
        # and scale it, so that it integrates to the number of electrons in ruthenium. 
        from confined_pseudo_atoms import ag
        atomdft.initialDensityGuess((ag.r, float(Nelec)/float(ag.Nelec) * ag.radial_density))
        """
        # If an initial density is available for Ru, we use that one
        from confined_pseudo_atoms import ru
        atomdft.initialDensityGuess((ru.r, float(Nelec)/float(ru.Nelec) * ru.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "ru.py", "w")
    atomdft.saveSolution(fh)
    fh.close()    
    
def ruthenium_2plus():
    energy_range = list(linspace(-1500.0, -800.0, 100)) \
        + list(linspace(-800.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 44
    Nelec = 42     # Ru^(2+)
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv, grid_spacing="exponential", damping=0.6)
    # The electron configuration of Ru^(2+) is  [Kr] 4d^(6)
    occupation = [occ for occ in occupation_numbers(Nelec)]
    # `occupation` is a list of tuples
    #   (nr.electrons, quantum number n + 1, quantum number l)
    # that describe the occupied shells.
    # In the default occupation, ruthenium would have the electronic configuration [Kr] 4d^(6),
    # which is correct.
    assert occupation[-1] == (6,4-1,2)
    
    atomdft.setOccupation(occupation)
    atomdft.setValenceOrbitals(["4d", "5s"], format="spectroscopic")
    atomdft.setRadialGrid(0.000000001, 14.0, 20000)
    try:
        # If no initial density is available for Ru^(2+), we can take the density from silver 
        # and scale it, so that it integrates to the number of electrons in ruthenium 2+. 
        from confined_pseudo_atoms import ag
        atomdft.initialDensityGuess((ag.r, float(Nelec)/float(ag.Nelec) * ag.radial_density))
        """
        # If an initial density is available for Ru^(2+), we use that one
        from confined_pseudo_atoms import ru
        atomdft.initialDensityGuess((ru.r, ru.radial_density))
        """
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    # Although we have calculated Ru^(2+), the pseudoorbitals still have to be saved under the name 'ru.py'.
    # Only one oxidation state can be used for any atom at the same time.
    fh = open(orbdir + "ru.py", "w")
    atomdft.saveSolution(fh)
    fh.close()    
    
def silver():
    energy_range = list(linspace(-1500.0, -800.0, 100)) \
        + list(linspace(-800.0, -200.0, 200)) \
        + list(linspace(-200.0, -30.0, 200)) \
        + list(linspace(-30.0, -5.0, 400)) \
        + list(linspace(-5.0, -0.0000001, 400)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 47
    Nelec = 47
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=1.0e-4, grid_spacing="exponential", damping=0.6)
    # In silver the 5s is filled before 4f: [Kr] 4d^(10) 5s^1
    occupation = [occ for occ in occupation_numbers(Nelec)]
    assert occupation[-1] == (1,4-1,3)
    occupation[-1] = (1,5-1,0)
    atomdft.setOccupation(occupation)
    # I would like to include unoccupied f orbitals in minimal basis but 
    # sofar no Slater rules exist for f orbitals. So include 5p
    atomdft.setValenceOrbitals(["4d", "5s","5p"], format="spectroscopic")    
    atomdft.setRadialGrid(0.000000001, 14.0, 20000)
    try:
        from free_pseudo_atoms import ag
        atomdft.initialDensityGuess((ag.r, ag.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "ag.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

    
def iodine():
    energy_range = list(linspace(-2000.0, -1400.0, 100)) \
        + list(linspace(-1400.0, -150.0, 100)) \
        + list(linspace(-150.0, -30.0, 100)) \
        + list(linspace(-30.0, -10.0, 100)) \
        + list(linspace(-10.0, 3.0, 100)) \
    # equidistant grid, use larger grid for heavy atoms
              
    Z = 53
    Nelec = 53
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=1.0e1, damping=0.5) #, grid_spacing="exponential", damping=0.5)    
    atomdft.setRadialGrid(0.000000001, 30.0, 10000)
    try:
        from free_pseudo_atoms import ag
        atomdft.initialDensityGuess((ag.r, ag.radial_density * 53.0/47.0))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    fh = open(orbdir + "i.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

    
## first row
#hydrogen()
#helium()
## second row
#lithium()
#beryllium()
#boron()
#carbon()
#nitrogen()
#oxygen()
#fluorine()
#neon()
## third row
#natrium()
#magnesium()
#aluminum()
#silicon()
#phosphorus()
#sulfur()
#chlorine()
#argon()
## fourth row
#potassium()

## heavy atoms, maybe the non-relativistic SE is not valid anymore?
# transition metals
scandium()
#titanium()
#iron()
#copper()
#zinc()
#bromine()
#ruthenium()
#ruthenium_2plus()
#silver()
#iodine()

################### PARTIALLY CHARGED ATOMS ######################

def hydrogen_charged():
    energy_range = linspace(-3.0, 3.0, 100)

    Z = 1
    Nelec = 1
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)#, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax, Npts)
    try:
        from free_pseudo_atoms import h
        atomdft.initialDensityGuess((h.r, h.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    print atomdft.Etot
    # removing partial charge from H does not have any effect on the energy
    # as there is no electron-electron interaction in as partially positively charged 1e atom.
    dq = linspace(0.0,0.1,30) # add negative charge => Strange jump at 0.0 because e-e energy
                              # is add discontinuously?
    atomdft.charge_fluctuations(dq)
    fh = open(orbdir + "h.py", "w")
    atomdft.saveSolution(fh)
    fh.close()


def helium_charged():
    energy_range = linspace(-10.0, 10.0, 100)

    Z = 2
    Nelec = 2
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv)#, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax+3.0, Npts+100)
    try:
        from free_pseudo_atoms import he
        atomdft.initialDensityGuess((he.r, he.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    print atomdft.Etot
    dq = linspace(0.0, 0.25, 7)
    atomdft.charge_fluctuations(dq)
    fh = open(orbdir + "he.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def carbon_charged():
    energy_range = linspace(-11, 1.0, 400)

    Z = 6
    Nelec = 6
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=2.0e-4, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax+10.0, Npts+1000)
    try:
        from free_pseudo_atoms import c
        atomdft.initialDensityGuess((c.r, c.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    dq = linspace(-0.1, 0.1, 10)
    atomdft.charge_fluctuations(dq)
    fh = open(orbdir + "c.py", "w")
    atomdft.saveSolution(fh)
    fh.close()

def nitrogen_charged():
    energy_range = linspace(-13, 1.0, 400)

    Z = 7
    Nelec = 7
    atomdft = PseudoAtomDFT(Z,Nelec,numerov_conv,en_conv=2.0e-4)#, grid_spacing="exponential")
    atomdft.setRadialGrid(rmin, rmax+10.0, Npts+1000)
    try:
        from free_pseudo_atoms import n
        atomdft.initialDensityGuess((n.r, n.radial_density))
    except ImportError:
        atomdft.initialDensityGuess()
    atomdft.setEnergyRange(energy_range)
    atomdft.solveKohnSham()
    dq = linspace(-0.1, 0.1, 20)
    atomdft.charge_fluctuations(dq)
    fh = open(orbdir + "n.py", "w")
    atomdft.saveSolution(fh)
    fh.close()


#hydrogen_charged()
#helium_charged()
#carbon_charged()
#nitrogen_charged()
