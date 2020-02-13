#!/usr/bin/env python
"""
Optimize a geometry using TD-DFTB, compute the Hessian
by numerical differentiation of the analytical gradients,
analyze the vibrations and sample a set of initial conditions
from the Wigner functions (saved to /tmp/dynamics_####.in)
"""
from DFTB.Dynamics import HarmonicApproximation
from DFTB.Molden import MoldenExporterSectioned
from DFTB.PES import PotentialEnergySurfaces
from DFTB import XYZ, utils, AtomicData
from DFTB import Thermochemistry
from DFTB.Timer import GlobalTimer as T

import numpy as np
import numpy.linalg as la

if __name__ == "__main__":
    ################################################################
    # The following code is only an example to illustrate
    # how the DFTB.PES class can be used. It performs an optimization
    # on an excited state.
    # You should call this program as 
    #      DFTB/optimize.py 
    # from the DFTB-### folder to ensure that the paths are correct
    ################################################################
    import sys
    from scipy import optimize

    if len(sys.argv) < 3:
        print "Usage: %s <xyz-file> I [H]" % sys.argv[0]
        print "  optimize geometry on the I-th electronic state"
        print "  and compute the hessian if the optional 3rd command line argument is H"
        print "  Type --help to see all options."
        print "  To reduce the amount of output add the option --verbose=0."
        print "  Examples:"
        print "      optimize.py molecule.xyz 0 H"
        print "  optimizes the molecule on the ground state and computes the Hessian."
        print "      optimize.py molecule.xyz 1"
        print "  finds the minimum on the 1st excited state."
        exit(-1)
    
    #
    xyz_file = sys.argv[1]    # path to xyz-file
    I = int(sys.argv[2])      # index of electronic state
    # Should the Hessian be calculated as well?
    calc_hessian = False
    if len(sys.argv) > 3:
        # optional 3rd argument
        if sys.argv[3].upper() == "H":
            calc_hessian = True
    # Read the geometry from the xyz-file
    atomlist = XYZ.read_xyz(xyz_file)[0]
    # read the charge of the molecule from the comment line in the xyz-file
    kwds = XYZ.extract_keywords_xyz(xyz_file)
    # initialize the TD-DFTB calculator
    pes = PotentialEnergySurfaces(atomlist, Nst=max(I+1,2), **kwds)

    # convert geometry to a vector
    x0 = XYZ.atomlist2vector(atomlist)

    # FIND ENERGY MINIMUM
    # f is the objective function that should be minimized
    # it returns (f(x), f'(x))
    def f(x):
        #
        if I == 0 and type(pes.tddftb.XmY) != type(None):
            # only ground state is needed. However, at the start
            # a single TD-DFT calculation is performed to initialize
            # all variables (e.g. X-Y), so that the program does not
            # complain about non-existing variables.
            enI, gradI = pes.getEnergyAndGradient_S0(x)
        else:
            energies, gradI = pes.getEnergiesAndGradient(x, I)
            enI = energies[I]
        print "E = %2.7f     |grad| = %2.7f" % (enI, la.norm(gradI))
        #
        # also save geometries from line searches
        save_xyz(x, info="energy=%s" % enI)

        return enI, gradI

    xyz_opt = xyz_file.replace(".xyz", "_opt.xyz")
    print "Intermediate geometries will be written to %s" % xyz_opt
    # This is a callback function that is executed by numpy for each optimization step.
    # It appends the current geometry to an xyz-file.
    def save_xyz(x, mode="a", info=""):
        atomlist_opt = XYZ.vector2atomlist(x, atomlist)
        XYZ.write_xyz(xyz_opt, [atomlist_opt], \
                      title="charge=%s %s" % (kwds.get("charge",0), info),\
                      mode=mode)
    # save initial geometry
    save_xyz(x0, mode="w")
        
    Nat = len(atomlist)
    options = {'gtol': 1.0e-7, 'norm': 2}
    # The "BFGS" method is probably better than "CG", but the line search in BFGS is expensive.
    res = optimize.minimize(f, x0, method="CG", jac=True, callback=save_xyz, options=options)
    #res = optimize.minimize(f, x0, method="BFGS", jac=True, callback=save_xyz, options=options)
    # save optimized geometry
    xopt = res.x
    Eopt = res.fun
    save_xyz(xopt, info="energy=%s" % Eopt)
    print "Optimized geometry written to %s" % xyz_opt

    if calc_hessian == True:
        # COMPUTE HESSIAN AND VIBRATIONAL MODES
        # The hessian is calculated by numerical differentiation of the 
        # analytical gradients
        def grad(x):
            if I == 0:
                enI, gradI = pes.getEnergyAndGradient_S0(x)
            else:
                energies, gradI = pes.getEnergiesAndGradient(x, I)
            return gradI
        print "Computing Hessian"
        hess = HarmonicApproximation.numerical_hessian_G(grad, xopt)
        np.savetxt("hessian.dat", hess)
        masses = AtomicData.atomlist2masses(atomlist)
        vib_freq, vib_modes = HarmonicApproximation.vibrational_analysis(xopt, hess, masses, \
                                                                         zero_threshold=1.0e-9, is_molecule=True)
        # compute thermodynamic quantities and write summary
        thermo = Thermochemistry.Thermochemistry(atomlist, Eopt, vib_freq, pes.tddftb.dftb2.getSymmetryGroup())
        thermo.calculate()
        
        # write vibrational modes to molden file
        molden = MoldenExporterSectioned(pes.tddftb.dftb2)
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

    # timing
    print T
    
