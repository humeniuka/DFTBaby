#!/usr/bin/env python
"""
Optimize a geometry using TD-DFTB, compute the Hessian
by numerical differentiation of the analytical gradients,
analyze the vibrations and sample a set of initial conditions
from the Wigner functions.
"""
from DFTB.Dynamics import HarmonicApproximation
from DFTB.Molden import MoldenExporterSectioned
from DFTB.PES import PotentialEnergySurfaces
from DFTB import XYZ, AtomicData
from DFTB import Thermochemistry
from DFTB.Timer import GlobalTimer as T
from DFTB.Optimize.InternalCoords import NotConvergedError, InternalValenceCoords
from DFTB.Optimize.Optimize import minimize

from DFTB import optparse
# We need to know the arguments of the __init__ functions for
# the following classes and expose them to the user as command line options.
from DFTB.DFTB2 import DFTB2
from DFTB.ImplicitSolvent import SolventCavity
from DFTB.LR_TDDFTB import LR_TDDFTB

from scipy import optimize

import numpy as np
import numpy.linalg as la

class GeometryOptimization:
    def __init__(self,
                 state=0, calc_hessian=0, coord_system="cartesian",
                 grad_tol=1.0e-5, func_tol=1.0e-8, max_steps=100000, method='CG',
                 explicit_bonds="[]",
                 freeze="[]", relaxed_scan="()"):
        """
        Parameters
        ==========
        Geometry Optimization.state: index of electronic state to be optimized (0 - ground state, 1 - first excited state).
        Geometry Optimization.calc_hessian: Should the hessian matrix be computed (1) or not (0)? If yes, the Hessian matrix is saved to the file 'hessian.dat' and the vibrational modes and frequencies are saved to the file 'vib.molden' that can be visualized with the Molden program.
        Geometry Optimization.coord_system: The optimization can be performed either directly in cartesian coordinate ('cartesian') or in redundant internal coordinates ('internal'). Cartesian coordinates are more reliable but make it difficult to converge to the minimum in a floppy molecule. Internal coordinates don't work for disconnected fragments and require a starting geometry with the correct atom connectivity.
        Geometry Optimization.grad_tol: The optimization is finished only if the norm of the gradient is smaller than this threshold.
        Geometry Optimization.func_tol: The optimization is finished only if the energy does not change more than this threshold.
        Geometry Optimization.max_steps: Maximum number of optimization steps.
        Geometry Optimization.method: Choose the optimization algorithm. 'Newton', 'Steepest Descent' and 'BFGS' have their own implementations, 'CG' requests scipy's conjugate gradient method.
        Geometry Optimization.explicit_bonds: Inserts artificial bonds between pairs of atoms. The bonds are specified as a list of tuples (I,J) of atom indices (starting at 1). This allows to join disconnected fragments.
        Geometry Optimization.freeze: Freeze internal coordinates. The internal coordinates that should be kept at their current value during the optimization are specified as a list of tuples of atom indices (starting at 1). Each tuple may contain 2, 3 or 4 atom indices, (I,J) - bond between atoms I and J, (I,J,K) - valence angle I-J-K, (I,J,K,L) - dihedral angle between bonds I-J, J-K and K-L. For example "[(1,2), (4,5,6)]" freezes the bond between atoms 1 and 2 and the angle 4-5-6. The atom indices do not necessarily have to correspond to a 'physical' bond, angle or dihedral. So, for instance, you can also freeze the distance between two atoms that are not connected. 
        Geometry Optimization.relaxed_scan: Perform a relaxed scan along an internal coordinate. The coordinate is incremented from its initial value by `nsteps` steps of size `incr`. In each step the value of the scan coordinate is kept constant while all other degrees of freedom are relaxed. The internal coordinate is specified by 2, 3 or 4 atom indicies as explained for the option `freeze` followed by the number of steps and the increment, which is in Angstrom for bond lengths and degrees for angles. The format is "(I,J, nsteps, incr)" for scanning a bond length, "(I,J,K, nsteps, incr)" for scanning a valence angle and "(I,J,K,L, nsteps, incr)" for scanning a torsion.  For example "(1,2, 5, 0.1)" will scan the bond between atom 1 and 2 in 5 steps of 0.1 Angstrom and "(1,2,3, 9, 10.0)" will scan the angle 1-2-3 in 9 steps of 10.0 degrees. Dihedral angles are limited to the range [0,180], so if the angle is close to 180 degrees a negative increment should be used, otherwise the scan will stop at 180 degrees.
        """
        self.state = state
        self.calc_hessian = calc_hessian
        assert coord_system in ["cartesian", "internal"]
        self.coord_system = coord_system
        self.grad_tol = grad_tol
        self.func_tol = func_tol
        self.maxiter = max_steps
        self.method = method
        # parameters for relaxed scan
        if relaxed_scan != ():
            assert coord_system == "internal", "A relaxed scan requires 'coord_system=internal'!"
            
            if len(relaxed_scan) == 4:
                # scan bond length
                I,J, nsteps, incr = relaxed_scan
                # convert increment from Angstrom to bohr
                incr /= AtomicData.bohr_to_angs
                IJKL = (I,J)
            elif len(relaxed_scan) == 5:
                # scan valence angle
                I,J,K, nsteps, incr = relaxed_scan
                # convert angle from degrees to radians
                incr *= np.pi/180.0
                IJKL = (I,J,K)
            elif len(relaxed_scan) == 6:
                # scan dihedral angle
                I,J,K,L, nsteps, incr = relaxed_scan
                # convert angle from degrees to radians
                incr *= np.pi/180.0
                IJKL = (I,J,K,L)
            else:
                raise ValueError("Format of relaxed scan '%s' not understood!" % relaxed_scan)

            # The scan coordinate has to be frozen in each scan step.
            freeze.append(IJKL)
            
            self.relaxed_scan_nsteps = int(nsteps)
            self.relaxed_scan_incr   = incr
            # shift indices by -1 so that the first index starts at 0
            self.relaxed_scan_IJKL   = tuple(map(lambda I: int(I)-1, IJKL))

            self.optimization_type = "relaxed_scan"
        else:
            self.optimization_type = "minimize"

        # freezing of internal coordinates
        self.freeze = []
        for IJKL in freeze:
            # Indices on the command line start at 1, but internally
            # indices starting at 0 are used.
            IJKL = tuple([I-1 for I in IJKL])
            self.freeze.append(IJKL)
            assert coord_system == "internal", "Freezing of internal coordinates require 'coord_system=internal'!"

            
    def setGeometry(self, atomlist, geom_kwds={}):
        self.geom_kwds = geom_kwds
        self.atomlist = atomlist
    def setOutput(self, xyz_opt="opt.xyz", xyz_scan="scan.xyz", dat_scan="scan.dat"):
        """files where geometries and energy tables created during the optimization and scan are written to"""
        self.xyz_opt = xyz_opt
        self.xyz_scan = xyz_scan
        self.dat_scan = dat_scan
    def getGeometry(self):
        """current geometry"""
        return self.atomlist
    def getEnergy(self):
        """current energy"""
        return self.enI
    def initialize(self):
        """
        This function should be called when the geometry is known (after calling setGeometry(...)).
        """
        # initialize the TD-DFTB calculator
        self.pes = PotentialEnergySurfaces(self.atomlist, Nst=max(self.state+1,2), **self.geom_kwds)
        # initialize internal coordinate system if needed
        if self.coord_system == "internal":
            self.IC = InternalValenceCoords(self.atomlist, freeze=self.freeze, verbose=self.pes.tddftb.dftb2.verbose)
        
    def minimize(self):
        I = self.state

        # convert geometry to a vector
        x0 = XYZ.atomlist2vector(self.atomlist)
    
        # This member variable holds the last energy of the state
        # of interest.
        self.enI = 0.0
        # last available energies of all electronic states that were
        # calculated
        self.energies = None
        
        # FIND ENERGY MINIMUM
        # f is the objective function that should be minimized
        # it returns (f(x), f'(x))
        def f_cart(x):
            #
            if I == 0 and type(self.pes.tddftb.XmY) != type(None):
                # Only ground state is needed. However, at the start
                # a single TD-DFT calculation is performed to initialize
                # all variables (e.g. X-Y), so that the program does not
                # complain about non-existing variables.
                enI, gradI = self.pes.getEnergyAndGradient_S0(x)
                energies = np.array([enI])
            else:
                energies, gradI = self.pes.getEnergiesAndGradient(x, I)
                enI = energies[I]
            self.enI = enI
            self.energies = energies
            print "E = %2.7f     |grad| = %2.7f" % (enI, la.norm(gradI))
            #
            # also save geometries from line searches
            save_xyz(x)

            return enI, gradI

        print "Intermediate geometries will be written to %s" % self.xyz_opt
        # This is a callback function that is executed for each optimization step.
        # It appends the current geometry to an xyz-file.
        def save_xyz(x, mode="a"):
            self.atomlist = XYZ.vector2atomlist(x, self.atomlist)
            XYZ.write_xyz(self.xyz_opt, [self.atomlist], \
                          title="charge=%s energy= %s" % (self.geom_kwds.get("charge",0), self.enI),\
                          mode=mode)
            return x
            
        Nat = len(self.atomlist)

        if self.coord_system == "cartesian":
            print "optimization is performed directly in cartesian coordinates"
            q0 = x0
            objective_func = f_cart
            save_geometry = save_xyz
            max_steplen = None
        elif self.coord_system == "internal":
            print "optimization is performed in redundant internal coordinates"
            # transform cartesian to internal coordinates, x0 ~ q0
            q0 = self.IC.cartesian2internal(x0)
        
            # define functions that wrap the cartesian<->internal transformations
            def objective_func(q):
                # transform back from internal to cartesian coordinates
                x = self.IC.internal2cartesian(q)
                self.IC.cartesian2internal(x)
                # compute energy and gradient in cartesian coordinates
                en, grad_cart = f_cart(x)
                # transform gradient to internal coordinates
                grad = self.IC.transform_gradient(x, grad_cart)
            
                return en, grad

            def save_geometry(q, **kwds):
                # transform back from internal to cartesian coordinates
                x = self.IC.internal2cartesian(q)
                # save cartesian coordinates
                save_xyz(x, **kwds)
                return x

            def max_steplen(q0,v):
                """
                find a step size `a` such that the internal->cartesian
                transformation converges for the point q = q0+a*v
                """
                a = 1.0
                for i in range(0, 7):
                    #print "trying step length a= %e   |a*v|= %e" % (a,la.norm(a*v))
                    q = q0 + a * v
                    try:
                        x = self.IC.internal2cartesian(q)
                    except NotConvergedError as e:
                        # reduce step size by factor of 1/2
                        a /= 2.0
                        continue
                    #print "final step length  a= %e   |a*v|= %e" % (a,la.norm(a*v))
                    break
                else:
                    raise RuntimeError("Could not find a step size for which the transformation from internal to cartesian coordinates would work for q=q0+a*v! Last step size a= %e  |v|= %e  |a*v|= %e" % (a, la.norm(v), la.norm(a*v)) )
                #print "max step length  a = %s" % a
                return a

        else:
            raise ValueError("Unknown coordinate system '%s'!" % self.coord_system)
        # save initial energy and geometry
        objective_func(q0)
        save_geometry(q0, mode="w")

        options = {'gtol': self.grad_tol, 'maxiter': self.maxiter, 'gtol': self.grad_tol, 'norm': 2}
        if self.method == 'CG':
            # The "BFGS" method is probably better than "CG", but the line search in BFGS is expensive.
            res = optimize.minimize(objective_func, q0, method="CG", jac=True, callback=save_geometry, options=options)
            #res = optimize.minimize(objective_func, q0, method="BFGS", jac=True, callback=save_geometry, options=options)

        elif self.method in ['Steepest Descent', 'Newton', 'BFGS']:
            # My own implementation of optimization algorithms
            res = minimize(objective_func, q0,
                           method=self.method,
                           #line_search_method="largest",
                           callback=save_geometry,
                           max_steplen=max_steplen,
                           maxiter=self.maxiter,
                           gtol=self.grad_tol,
                           ftol=self.func_tol)
        else:
            raise ValueError("Unknown optimization algorithm '%s'!" % self.method)
            
        # save optimized geometry
        qopt = res.x
        Eopt = res.fun
        xopt = save_geometry(qopt)
        print "Optimized geometry written to %s" % self.xyz_opt
            
        if self.calc_hessian == 1:
            # COMPUTE HESSIAN AND VIBRATIONAL MODES
            # The hessian is calculated by numerical differentiation of the 
            # analytical cartesian gradients
            def grad(x):
                en, grad_cart = f_cart(x)
                return grad_cart
            print "Computing Hessian"
            hess = HarmonicApproximation.numerical_hessian_G(grad, xopt)
            np.savetxt("hessian.dat", hess)
            masses = AtomicData.atomlist2masses(atomlist)
            vib_freq, vib_modes = HarmonicApproximation.vibrational_analysis(xopt, hess, masses, \
                                                                             zero_threshold=1.0e-9, is_molecule=True)
            # compute thermodynamic quantities and write summary
            thermo = Thermochemistry.Thermochemistry(atomlist, Eopt, vib_freq, self.pes.tddftb.dftb2.getSymmetryGroup())
            thermo.calculate()
        
            # write vibrational modes to molden file
            molden = MoldenExporterSectioned(self.pes.tddftb.dftb2)
            atomlist_opt = XYZ.vector2atomlist(xopt, atomlist)
            molden.addVibrations(atomlist_opt, vib_freq.real, vib_modes.transpose())
            molden.export("vib.molden")

        ## It's better to use the script initial_conditions.py for sampling from the Wigner
        ## distribution
        """
        # SAMPLE INITIAL CONDITIONS FROM WIGNER DISTRIBUTION
        qs,ps = HarmonicApproximation.initial_conditions_wigner(xopt, hess, masses, Nsample=200)
        HarmonicApproximation.save_initial_conditions(atomlist, qs, ps, ".", "dynamics")
        """

    def relaxed_scan(self, IJKL, nsteps, incr):
        """
        perform a relaxed scan along the internal coordinate IJKL. The coordinate is 
        incremented from its initial value by `nsteps` steps of size `incr`. In each
        step the value of the scan coordinate is kept constant while all other degrees 
        of freedom are relaxed. 

        Parameters
        ----------
        IJKL    :  tuple of 2, 3 or 4 atom indices (starting at 0)
                   (I,J)     -   bond between atoms I and J
                   (I,J,K)   -   valence angle I-J-K
                   (I,J,K,L) -   dihedral angle between the bonds I-J, J-K and K-L
        nsteps  :  number of steps
        incr    :  increment in each step, in bohr for bond lengths, in radians
                   for angles
        """
        # length of tuple IJKL determines type of internal coordinate
        typ = len(IJKL)
        coord_type = {2: "bond", 3: "angle", 4: "dihedral"}
        conv_facs = {2: AtomicData.bohr_to_angs, 3: 180.0/np.pi, 4: 180.0/np.pi}
        units = {2: "Angs", 3: "degs", 4: "degs"}
        #
        assert self.coord_system == "internal"
        # freeze scan coordinate at its current value
        self.IC.freeze(IJKL)
        print "  ============ "
        print "  RELAXED SCAN "
        print "  ============ "
        print "  The internal coordinate defined by the atom indices"
        print "    IJKL = %s  " % map(lambda I: I+1, IJKL)
        print "  is scanned in %d steps of size %8.5f %s." % (nsteps, incr*conv_facs[typ], units[typ])
        
        def save_step():
            """function is called after each minimization"""
            scan_coord = self.IC.coordinate_value(xi, IJKL)
            
            print "current value of scan coordinate : %s" % scan_coord

            if i == 0:
                mode = "w"
            else:
                mode = "a"
            # save relaxed geometry of step i
            XYZ.write_xyz(self.xyz_scan, [atomlist],
                          title="charge=%s energy=%s" % (self.geom_kwds.get("charge",0), self.enI), mode=mode)

            # save table with energies along scan
            fh = open(self.dat_scan, mode)
            if i == 0:
                # write header
                print>>fh, "# Relaxed scan along %s defined by atoms %s" % (coord_type[typ], map(lambda I: I+1, IJKL))
                print>>fh, "# state of interest: %d" % self.state
                print>>fh, "# "
                print>>fh, "#  Scan coordinate     Energies "
                print>>fh, "#    %s              Hartree " % units[typ]
                
            print>>fh, "  %8.5f     " % scan_coord,
            for en in self.energies:
                print>>fh, "   %e " % en,
            print>>fh, ""
            fh.close()
            

        for i in range(0, nsteps):
            print "Step %d of relaxed scan" % i
            # relax all other coordinates
            self.minimize()
            # optimized geometry of i-th step
            atomlist = self.getGeometry()
            xi = XYZ.atomlist2vector(atomlist)
            # save geometry
            save_step()
            # take a step of size `incr` along the scan coordinate
            xip1 = self.IC.internal_step(xi, IJKL, incr)
            # update geometry
            atomlist = XYZ.vector2atomlist(xip1, atomlist)
            self.setGeometry(atomlist, geom_kwds=self.geom_kwds)

        print "Scan geometries were written to %s" % self.xyz_scan
        print "Table with scan energies was written to %s" % self.dat_scan
        
    def optimize(self):
        """
        run minimization of energy or relaxed scan
        """
        if self.optimization_type == "minimize":
            self.minimize()
        elif self.optimization_type == "relaxed_scan":
            self.relaxed_scan(self.relaxed_scan_IJKL,
                              self.relaxed_scan_nsteps,
                              self.relaxed_scan_incr)
        else:
            raise ValueError("BUG? optimization_type = %s" % self.optimization_type)
            
if __name__ == "__main__":
    import sys
    import os.path

    usage = """

   Usage: %s  molecule.xyz

      optimize the geometry and compute the hessian matrix

      Intermediate geometries are written to 'molecule_opt.xyz'. The last
      geometry is the optimized one.

      Type --help to see all options.
      To reduce the amount of output add the option --verbose=0.
      
      Examples:
          GeometryOptimization.py molecule.xyz --state=0 --calc_hessian=1
      optimizes the molecule on the ground state and computes the Hessian.
          GeometryOptimization.py molecule.xyz --state=1 --coord_system='internal'
      finds the minimum on the 1st excited state using redundant internal coordinates.
    """ % os.path.basename(sys.argv[0])

    # This wrapper makes the optional parameters of the python functions visible
    # as optional command line argument.
    parser = optparse.OptionParserFuncWrapper(
        [
            # options for geometry optimization
            GeometryOptimization.__init__,
            # electronic structure options
            DFTB2.__init__, 
            DFTB2.runSCC,
            SolventCavity.__init__,
            LR_TDDFTB.getEnergies
        ],
        usage, section_headers=["DFTBaby", "GeometryOptimization"])
    # extract optional parameters from command line
    (options,args) = parser.parse_args(GeometryOptimization.__init__)

    if len(args) < 1:
        print usage
        exit(-1)
    
    GOpt = GeometryOptimization(**options)
    #
    xyz_file = args[0]
    # Read the geometry from the xyz-file
    atomlist = XYZ.read_xyz(xyz_file)[0]
    # read the charge of the molecule from the comment line in the xyz-file
    geom_kwds = XYZ.extract_keywords_xyz(xyz_file)

    # set initial geometry
    GOpt.setGeometry(atomlist, geom_kwds=geom_kwds)
    # initialize calculator and internal coordinates
    GOpt.initialize()

    # output files for optimizations and scans
    xyz_opt = xyz_file.replace(".xyz", "_opt.xyz")
    xyz_scan = xyz_file.replace(".xyz", "_scan.xyz")
    dat_scan = xyz_file.replace(".xyz", "_scan.dat")
    GOpt.setOutput(xyz_opt=xyz_opt,
                   xyz_scan=xyz_scan,
                   dat_scan=dat_scan)
    
    # run optimization
    GOpt.optimize()
    
    # optimized geometry, we don't need it here actually, the final
    # geometry in <molecule>_opt.xyz is the optimized one
    atomlist_opt = GOpt.getGeometry()

    print "final energy: %s Hartree" % GOpt.getEnergy()
    # timing
    print T
    
    print "FINISHED"
    
