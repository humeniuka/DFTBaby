#!/usr/bin/env python
"""
The minimal energy conical intersection (MECI) is optimized using Robb's and Bearpark's algorithm.

  M. Bearpark, M. Robb, B. Schlegel,
  "A direct method for the location of the lowest energy point on a potential surface crossing",
  Chem. Phys. Lett. 223 (1994) 269-274.

At the conical intersection SCF will fail. To avoid this, the energy shift method by Harabuchi et. al.
is used:

  Y. Harabuchi, K. Saita, S. Maeda,
  "Exploring radiative and nonradiative decay paths in indole, isoindole, quinoline, and isoquinoline"
  Photochem. Photobiol. Sci., 2018, DOI: 10.1039/C7PP00385D
"""

from DFTB import AtomicData, XYZ
from DFTB.PES import PotentialEnergySurfaces
from DFTB.Optimize.InternalCoords import InternalValenceCoords

import numpy as np
import numpy.linalg as la

class MECI:
    def __init__(self, atomlist0, pes, state1=0, state2=1,
                 coord_system="cartesian",
                 c1=0.9, c2=1.0, epsilon=0.007):
        """
        minimal energy conical intersection (MECI)
        
        Parameters
        ----------
        atomlist0      : initial geometry
        pes            : instance of DFTB.PES.PotentialEnergySurfaces

        Optional
        --------
        state1, state2 : indices of lower and upper electronic states 
                         between which the MECI should be found 
                         (0 - ground state, 1 - first excited state)
        coord_system   : The step along the gradient is taken either directly
                         in 'cartesian' coordinates or in 'internal' redundant
                         coordinates.
        epsilon        : To avoid hitting the CI seam exactly, which would
                         result in convergence problems in the SCF cycle, 
                         an energy gap of `epsilon` is maintained between
                         the upper and the lower states.
        """
        self.atomlist0 = atomlist0
        self.pes = pes
        self.state1 = state1
        self.state2 = state2
        assert self.state1 < self.state2
        self.coord_system = coord_system
        if self.coord_system == "internal":
            self.ic = InternalValenceCoords(self.atomlist0, verbose=self.pes.tddftb.dftb2.verbose)
        # parameters of MECI search algorithm
        self.c1 = c1
        self.c2 = c2
        self.epsilon = epsilon
    def runTDDFTB(self, x):
        """
        Parameters
        ----------
        x       :  cartesian coordinates

        Returns
        -------
        e1,e2   :  total energies of 1st and 2nd electronic state (in Hartree)
        g1,g2   :  gradients of total energies (in a.u.)
        nac     :  non-adiabatic coupling vector between 1st and 2nd state
        """
        # compute energies and gradient of lower state
        energies1, grad1 = self.pes.getEnergiesAndGradient(x, self.state1)
        e1 = energies1[self.state1]
        g1 = grad1
        # compute energies and gradient of lower state
        energies2, grad2 = self.pes.getEnergiesAndGradient(x, self.state2)
        e2 = energies2[self.state2]
        g2 = grad2
        #
        if self.state1 == 0:
            # An approximate NAC vector can be computed for transitions
            # to the ground state.
            nac = self.pes.tddftb.NonAdiabaticCouplingVector(self.state2-1)
            nac = nac.flatten()
        else:
            # NAC vectors between excited states are not available. The optimization
            # algorithm seems to work with any random vector that is not parallel
            # to the gradient. 
            nac = np.ones(len(x))
            
        return e1,e2, g1,g2, nac
        
    def getGradient(self, x):
        """
        The problem with the Bearpark algorithm is that there is no objective function,
        only a gradient. Following in the direction of this vector leads to the MECI.
        """
        # compute electronic structure
        e1,e2,g1,g2,nac = self.runTDDFTB(x)
        # normalized gradient difference vector
        x1 = g2-g1
        x1 /= la.norm(x1)
        # make non-adiabatic coupling vector orthogonal to gradient difference vector
        x2 = nac - np.dot(x1,nac) * x1
        # and normalize it
        x2 /= la.norm(x2)
        # check that x1 and x2 are really orthogonal
        assert abs(np.dot(x1,x2)) < 1.0e-10, "x1 and x2 not orthogonal!"

        # P is the projector onto the 3*N-2 dimensional orthogonal complement to
        # the plance x1,x2
        P = np.eye(len(x1)) - np.outer(x1,x1) - np.outer(x2,x2)

        f = 2*(e2-e1-self.epsilon)*x1
        
        # project gradient dE2/dq onto seam space
        g = np.dot(P, g2)

        #grad = self.c2 * ( self.c1*g + (1.0-self.c1)*f )
        grad = g + f
        
        return e1,e2, grad

    def adjust_shift(self, en_gap):
        """ adjust energy shift depending on the gap between the two states"""
        en_gap_eV = en_gap * AtomicData.hartree_to_eV
        if en_gap_eV > 0.4:
            self.c1 = 0.5
            self.c2 = 2.0
        if en_gap_eV <= 0.4:
            self.c1 = 0.7
            self.c2 = 1.0
        if en_gap_eV <= 0.2:
            self.c1 = 0.9
            
    def opt(self, max_iter=5000, step_size=1.0, gtol=0.005):
        """
        optimize MECI by following the gradient downhill

        Optional
        --------
        max_iter  :  maximum number of steps
        step_size :  The geometry is updated by making a step 
                       x -> x - step_size * grad
        gtol      :  tolerance for the gradient norm
        """
        print "optimize MECI by following the gradient"
        print "  lower state :  %d" % state1
        print "  upper state :  %d" % state2
        print "Intermediate geometries are written to 'meci_path.xyz'"
        print "and a table with energies is written to 'meci_energies.dat'"

        # initial geometry
        x = XYZ.atomlist2vector(self.atomlist0)
        # overwrite geometries from previous run
        mode = "w"
        en_fh = open("meci_energies.dat", "w")
        print>>en_fh, "# optimization of MECI"
        print>>en_fh, "# STEP     ENERGY(1)/Hartree    ENERGY(2)/Hartree    ENERGY(3)-EPS/Hartree"
        
        for i in range(0, max_iter):
            e1, e2, grad = self.getGradient(x)
            # energy gap
            en_gap = e2-e1
            self.adjust_shift(en_gap)

            # save intermediate steps and energies
            atomlist = XYZ.vector2atomlist(x, self.atomlist0)
            XYZ.write_xyz("meci_path.xyz", [atomlist],
                          title="ENERGY= %e  GAP= %e" % (e2, en_gap), mode=mode)
            print>>en_fh, " %4.1d      %+15.10f      %+15.10f      %+15.10f" % (i,e1,e2,e2-self.epsilon)
            en_fh.flush()
            
            # append to trajectory file
            mode = "a"
            #
            
            gnorm = la.norm(grad)
            print " %4.1d    e2= %e  e2-e1= %e   |grad|= %e  (tolerance= %e)" % (i, e2, en_gap, gnorm, gtol)
            if gnorm < gtol:
                break

            if self.coord_system == "cartesian":
                # descend along gradient directly in cartesian coordinates
                x -= step_size * grad
            elif self.coord_system == "internal":
                # use internal redundant coordinates
                # 1) transform cartesian to internal coordinates x -> q
                q = self.ic.cartesian2internal(x)
                # 2) transform cartesian gradient to internal coordinates
                # dE/dx -> dE/dq
                grad_intern = self.ic.transform_gradient(x, grad)
                # 3) take step along gradient in internal coordinates
                q = q - step_size * grad_intern
                # 4) deduce new cartesian coordinates from new internal coordinates q
                x = self.ic.internal2cartesian(q)

        else:
            print "exceeded maximum number of steps"
            
        en_fh.close()
        
            
if __name__ == "__main__":
    import sys
    from os.path import basename
    
    usage = """
  Usage: %s  molecule.xyz  state1  state2 

    find the minimal energy conincal intersection between two
    electronic states

  Arguments:

    molecule.xyz        -   initial geometry in XYZ format
    state1              -   integer, lower electronic state (0 - ground state)
    state2              -   integer, upper electronic state (1 - 1st excited state)
    
  Output Files:

    meci_path.xyz      -  intermediate geometries
    meci_energies.dat  -  energies of lower and upper states 
                          for intermediate geometries

  Type --help to see options for DFTB.

    """ % basename(sys.argv[0])
    
    args = sys.argv[1:]

    if len(args) < 3:
        print usage
        exit(-1)
    
    xyz_file = args[0]         # path to xyz-file
    state1 = int(args[1])      # index of lower electronic state
    state2 = int(args[2])      # index of upper electronic state

    # load initial geometry, we take the last geometry, so it is
    # easier to restart a previous MECI calculation.
    atomlist0 = XYZ.read_xyz(xyz_file)[-1]
    # read the charge of the molecule from the comment line in the xyz-file
    kwds = XYZ.extract_keywords_xyz(xyz_file)
    # initialize the TD-DFTB calculator
    pes = PotentialEnergySurfaces(atomlist0, Nst=state2+1, **kwds)

    meci = MECI(atomlist0, pes,
                #coord_system='internal',
                coord_system='cartesian', 
                state1=state1, state2=state2)
    
    meci.opt(step_size=0.1, gtol=0.001)

