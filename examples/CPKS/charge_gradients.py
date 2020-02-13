#!/usr/bin/env python
#
# run this script as
#
#    ./charge_gradients.py water.xyz
#

def calculate_charge_derivative():
    """compute analytical gradients of Mulliken charges"""
    #
    # This long prologue serves for reading all options
    # from the command line or the dftbaby.cfg file.
    # The command-line option
    #
    #    --cpks_solver='direct' or 'iterative'
    #
    # determines whether the full CPKS matrix should be constructed ('direct')
    # or whether the system of linear equations should be solved iteratively using the Krylov
    # subspace method.
    #
    from DFTB.XYZ import read_xyz, extract_keywords_xyz
    from DFTB.DFTB2 import DFTB2
    from DFTB import utils
    from DFTB.LR_TDDFTB import LR_TDDFTB
    from DFTB.ExcGradients import Gradients
    
    import sys
    import os

    usage  = "Usage: %s <xyz-file>\n" % sys.argv[0]
    usage += "   --help option will give more information\n"

    parser = utils.OptionParserFuncWrapper([\
       DFTB2.__init__, DFTB2.runSCC, \
       LR_TDDFTB.getEnergies, LR_TDDFTB.saveAbsorptionSpectrum, LR_TDDFTB.analyseParticleHole, \
       Gradients.getGradients], \
                usage)

    (options,args) = parser.parse_args(DFTB2.__init__)

    if len(args) < 1:
        print usage
        exit(-1)

    xyz_file = args[0]
    atomlist = read_xyz(xyz_file)[0]
    # number of atoms
    Nat = len(atomlist)
    kwds = extract_keywords_xyz(xyz_file)

    tddftb = LR_TDDFTB(atomlist, **options)

    (options,args) = parser.parse_args(tddftb.getEnergies)
    (scf_options,args) = parser.parse_args(tddftb.dftb2.runSCC)
    options.update(scf_options)
    tddftb.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
    tddftb.getEnergies(**options)

    #
    # The important part starts here.
    #
    grad = Gradients(tddftb)
    # A CPKS calculation has to be preceeded by a gradient calculation. 
    # The option `save_intermediates_CPKS=1` tells the program to save intermediate
    # variables that are needed later during the CPKS calculation.
    grad.gradient(I=0, save_intermediates_CPKS=1)

    # This runs a CPKS calculation and computes the gradients of the charges
    dQdp = grad.getChargeGradients()
    # The gradient of the charge on atom B w/r/t the position of atom A is
    #   d(Q_B)/dR_A = dQdp[3*A:3*(A+1), B]

    A = 0       # first atom
    B = Nat-1   # last atom
    print "Gradient of Mulliken charge on atom A=%d w/r/t nucleus B=%d  d(Q_B)/dR_A = %s" \
        % (A, B, dQdp[3*A:3*(A+1), B])

    
if __name__ == "__main__":
    calculate_charge_derivative()    

    
