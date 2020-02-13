#!/usr/bin/env python
"""
evaluate force field and exciton model at a fixed geometry
"""

import numpy.linalg as la

from DFTB import XYZ, AtomicData

from DFTB.ForceField.PeriodicForceField import PeriodicForceField, read_force_field, read_transition_charges, save_exciton_spectrum, plot_exciton_spectrum

if __name__ == "__main__":
    import numpy as np
    from DFTB import XYZ, utils
    import sys
    from optparse import OptionParser

    usage  = "Usage: python %s <.ff file>\n" % sys.argv[0]
    usage += "  constructs force field and computes total energy and\n"
    usage += "  exciton spectrum if transition charges are available.\n"
    usage += "  see --help for all options\n"

    parser = OptionParser(usage)
    parser.add_option("--transition_charges", dest="transition_charges", type=str, default="", help="Path to .chromo file with excitation energies, transition charges and magnetic transition dipoles for building a Frenckel exciton model. [default: %default]")
    parser.add_option("--spectrum_file", dest="spectrum_file", type=str, default="exciton_spectrum.dat", help="Save excitonic absorption and circular dichroism spectra to this file [default: %default]")
    parser.add_option("--state", dest="state", type=int, default=0, help="Excitonic state for which the gradient should be calculated (0 for ground state). [default: %default]")
    parser.add_option("--plot_spectra", dest="plot_spectra", type=str, default="", help="Choose units ('nm', 'Hartree', 'cm-1' or 'eV') for plotting absorption and circular dichroism spectra. [default: %default]")
    
    (opts,args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)

    ff_file = args[0]  #"h2.ff" #"ethene.ff" #"pyrene_crystal_expanded.ff" #
    # read force field definition
    atomlist, atomtypes, partial_charges, lattice_vectors = read_force_field(ff_file)
    # read transition charge for exciton model (if available)
    if opts.transition_charges != "":
        chromophores = list(read_transition_charges(opts.transition_charges))
    else:
        chromophores = []
    pff = PeriodicForceField(atomlist, atomtypes, partial_charges, lattice_vectors, chromophores)
    coords = XYZ.atomlist2vector(atomlist)

    # evaluate force field once
    energy, grad = pff.getEnergyAndGradient(coords, state=opts.state)
    print "Total energy: %s" % energy
    print "|gradient|  : %s" % la.norm(grad)  
    # compute exciton spectrum
    en, T, M = pff.getTransitionDipoles(verbose=1)
    save_exciton_spectrum(opts.spectrum_file, en, T, M)
    #
    if opts.plot_spectra != "":
        plot_exciton_spectrum(en, T, M, units=opts.plot_spectra)
