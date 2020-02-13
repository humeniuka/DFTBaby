#!/usr/bin/env python
"""
fit point charges at atomic centers to an electrostatic potential (ESP) provided in a cube file 
using the CHELP(G) algorithm.

Examples
--------
For the exciton model we need the partial charges for 
the ground state and the transition charges for an excited state.
The electrostatic potential for the ground state and excited state densities
(but as far as I know NOT the transition densities) can be extracted from 
the Gaussian 09 checkpoint file. 
The input file below computes the lowest 2 excited states of water with DFT. 

------ water.gjf ------
%Nproc=3
%Chk=water.chk
#P PBEPBE/6-31+G* TD=(Nstates=2,Singlets,Root=1) Density=(Transition=1) Pop=CHELPG  IOp(9/40=5)

water

0,1
  H       -0.29943        0.91863        0.00000
  O       -0.01314       -0.01859        0.00000
  H        0.96591        0.02391        0.00000



------ end of file ------

Run the DFT calcuation:

   g09 < water.gjf > water.out
   formchk water.chk 

Extract the electrostatic potential for the ground state and fit partial charges:

   cubegen 0 potential=SCF water.fchk esp_S0.cube -4 h
   chelpg.py esp_S0.cube partial_charges_S0.dat

This should give the partial charges -0.794 for O and +0.397 for both H-atoms.

Or extract the ESP for 1st excited state and fit partial charges:

   # replace CI Rho(Root) by CI in the formatted checkpoint file
   # If you get the error 'FCGetN is confused' you need to check the number
   # of white spaces so that the data type flag in the formatted checkpoint 
   # file is at the correct position.
   sed -i "s/CI Rho(1) Density/CI Density       /g" water.fchk
   cubegen 0 potential=CI water.fchk esp_S1.cube -4 h
   chelpg.py esp_S1.cube partial_charges_S1.dat
   
This should give the partial charges +0.292 for O and -0.146 for both H-atoms.

Transition densities between ground and excited states can be extracted from formatted checkpoint
files using the 'Multiwfn' program. The IOp(9/40=5) should be added to the route section
so that all excitation and deexcitation coefficients larger than 10^(-5) are written
to the log-file. 

Approximate transition charges can also be obtained from a TD-DFTB calculation 
(using the option --save_transition_charges).

References
----------

 CHELP:  Chirlian,L. Francl, M.  "Atomic Charges Derived from Electrostatic Potentials: A Detailed Study"
         J. Comput. Chem., Vol. 8, No. 6, 894-905 (1987)

 CHELPG: Breneman, C. Wiberg, K. "Determining Atom-Centered Monopoles from Molecular Electrostatic Potentials. The Need for High Sampling Density in Formamide Conformational Analysis"
         J. Comput. Chem., Vol. 11, No. 3, 361-373 (1990)
"""
import numpy as np
from numpy import linalg as la 
from DFTB import AtomicData, XYZ
from DFTB.Analyse import Cube

class CHELPG:
    def __init__(self):
        pass
    def setFitPoints(self, r, V):
        """
        Parameters
        ----------
        r: numpy array with shape (3,Npts), fit points
        V: numpy array with shape (Npts), 
           V[i] is the electrostatic potential at the fit point (x,y,z) = r[:,i]
        """
        self.r = r
        self.V = V
        self.m = len(V)  # number of fit points
    def setAtomicCenters(self, atomlist, charge=0):
        """
        set centers where the monopoles are located

        Parameters
        ----------
        atomlist: list of tuples (Z,[x,y,z]) with atomic number Z and center [x,y,z] (in bohr)

        Optional
        --------
        charge: total charge
        """
        self.n = len(atomlist)  # number of centers
        self.atomlist = atomlist
        self.charge = charge
        # R[:,a] is the 3d position of the a-th center
        self.R = np.zeros((3,self.n))
        # vdw radius for each center
        self.vdw_radii = np.zeros(self.n)
        for a,(Za,posa) in enumerate(self.atomlist):
            self.R[:,a] = np.array(posa)
            self.vdw_radii[a] = AtomicData.vdw_radii[ AtomicData.atom_names[Za-1] ]
    def fitESP(self, verbose=1):
        """
        fit monopoles to electrostatic potential. 

        Returns
        -------
        q: numpy array with shape (n), optimized point charges for each atom


        Minimize the deviation

              sum | V^esp(r) - V^monopole(r) |^2 

        subject to the constrain that the total charge is equal to charge.

        The stationary points of the Lagrangian 

              L(q1,...,qn) = 1/2 sum_i=1^m ( V[i] - sum_j=1^nat q[j]/|r[:,i] - R[j]| )^2  + lambda * (charge - sum_j=1^nat q_j)

        are found by solving the following system of linear equations

                        j=1,...,n

           (  __m                          -1  )  ( q1      )         (  __m           )
           (  \            1               -1  )  ( .       )         ( \       V[i]   )
           (  /      -------------         -1  )  ( qi      )     =   ( /     -------- )
           (  --a=1   d[a,i] d[a,j]         .. )  ( .       )         ( --a=1  d[a,i]  )  i=1,...,n
           (                               -1  )  ( qn      )         (                )
           (  1  1  1 ..................1   0  )  ( lambda  )         (      -charge   )  


                            A                   *      x           =          b

        The first n elements of the solution vector x are the point charges
        """
        # _d[i,j] is the distance between the i-th fit point and the j-th center, |r_i - R_j|
        print "compute distances..."
        _d = np.zeros((self.m, self.n))
        for xyz in range(0, 3):
            for j in range(0, self.n):
                _d[:,j] += (self.r[xyz,:] - self.R[xyz,j])**2
        _d = np.sqrt(_d)

        print "eliminate those fit points that are either"
        print "  - closer than the van der Waals radius to any atom"
        print "  - or farther away than Rmax=2.8 Angstrom from all atoms"
        Rmax = 2.8 / AtomicData.bohr_to_angs
        d = []  # list of fit points outside vdW radius and inside volume with max distance Rmax from each atom
        V = []  # electrostatic potential at the selected fit points
        print "vdW radii: %s" % self.vdw_radii
        print "Rmax= %s bohr" % Rmax
        for i in range(0, self.m): 
            if (_d[i,:] >= self.vdw_radii).all() and (_d[i,:] < Rmax).any():
                d.append( _d[i,:] )
                V.append( self.V[i] )
        m = len(d)
        print "%d of %d fit points left" % (m, self.m)
        d = np.array(d)
        V = np.array(V)

        print "build linear equation system..."
        # A matrix
        A = np.zeros((self.n+1,self.n+1))
        # right hand side of A.x = b
        b = np.zeros(self.n+1)
        #
        for i in range(0, self.n):
            A[i,i] = np.sum( 1.0/(d[:,i]*d[:,i]) )
            for j in range(i+1, self.n):
                A[i,j] = np.sum( 1.0/(d[:,i]*d[:,j]) )
                A[j,i] = A[i,j]
            b[i] = np.sum(V/d[:,i])
        #
        b[-1] = -self.charge
        # last row of A
        A[-1,:self.n] =  1.0
        # last column of A
        A[:self.n,-1] = -1.0
        A[-1,-1] = 0.0

        print "solve A.x = b..."
        x = la.solve(A,b)
        q = x[:-1]

        # root mean square error
        rmsd = np.sqrt( np.sum( (V - np.sum(q/d, axis=1))**2 ) / float(m) )
        # dipole of charge distribution
        dip = sum( [q[i] * self.R[:,i] for i in range(0, self.n)] )
        
        if verbose > 0:
            print ""
            print "root mean square deviation V(electrostatic) - V(monopoles): %e" % rmsd
            print "sum of charges: %s" % sum(q)
            print "dipole vector of charge distribution: %+7.5f  %+7.5f  %+7.5f" % tuple(dip.tolist())
            print ""
            
            print "  #       Element      Monopole charge"
            print "========================================="
            for i,((Zi, posi),qi) in enumerate(zip(self.atomlist, q)):
                print "  %4.1d     %4.2s           %+7.5f" % (i, AtomicData.atom_names[Zi-1].capitalize(), qi)

            print ""
            
        return q

if __name__ == "__main__":
    import os.path
    import sys
    from optparse import OptionParser
    
    usage = "Usage: python %s <.cube file> <.chg output>\n" % os.path.basename(sys.argv[0])
    usage += "  fit monopoles to electrostatic potential (ESP) in a cube file\n"
    usage += "  with the CHELPG method.\n"
    usage += "  The partial charges are written as a column in the output .dat file.\n"
    usage += "  --help shows all options\n"
    
    parser = OptionParser(usage)
    parser.add_option("--fit_centers", dest="fit_centers", type=str, default="", help="Path to an xyz-file with atomic positions that should be used as fit centers instead of the atomic centers in the cube-files. [default: %default]")
    parser.add_option("--format", dest="format", type=str, default="charges", help="The charges can be saved to the file either as a single column ('charges') or together with the atomic positions in the xyz-format ('xyz')")
    
    (opts, args) = parser.parse_args()
    if len(args) < 2:
        print usage
        exit(-1)
        
    esp_cube_file = args[0]
    chg_file = args[1]
    
    # load cube file with electrostatic potential (esp)
    print "load ESP from cube file..."
    atomlist, origin, axes, data = Cube.readCube(esp_cube_file)
    if opts.fit_centers != "":
        print "loading fit centers from '%s'" % opts.fit_centers
        atomlist = XYZ.read_xyz(opts.fit_centers)[0]
    # extract positions of points and values of ESP
    points, epot = Cube.get_points_and_values(origin, axes, data)

    # fit point charges
    chelpg = CHELPG()
    chelpg.setFitPoints(points, epot)
    chelpg.setAtomicCenters(atomlist, charge=0)
    partial_charges = chelpg.fitESP()

    # save point charges
    if opts.format == 'charges':
        np.savetxt(chg_file, partial_charges)
    else:
        XYZ.write_charges(chg_file, atomlist, partial_charges)
    print "partial charges saved to '%s'" % chg_file
