#!/usr/bin/env python
"""
Optimize the geometry of a molecular crystal using a periodic force field.
"""
import numpy as np
import numpy.linalg as la
import os.path

from DFTB import XYZ
from DFTB.ForceField.PeriodicForceField import read_force_field, read_transition_charges, PeriodicForceField, enlarged_unitcell, save_exciton_spectrum

if __name__ == "__main__":
    import sys
    from scipy import optimize
    import optparse

    usage  = "Usage: %s <.ff force field file>\n" % os.path.basename(sys.argv[0])
    usage += "  optimize geometry using periodic force field\n"
    usage += "  Force field files can be found in the folder DFTB/ForceField/examples/\n"
    usage += "  type --help to see all options.\n"

    parser = optparse.OptionParser(usage)
    parser.add_option("--nr_unit_cells", dest="nr_unit_cells", type=int, default=0, help="In the output the N neighbouring unit cells are added to the geometry to give an impression of a larger block of the crystal. [default: %default]")
    parser.add_option("--save_every", dest="save_every", type=int, default=1, help="Save every N-th optimization step. [default: %default]")
    parser.add_option("--transition_charges", dest="transition_charges", type=str, default="", help="Path to .chromo file with excitation energies, transition charges and magnetic transition dipoles for building a Frenckel exciton model. [default: %default]")
    parser.add_option("--spectrum_file", dest="spectrum_file", type=str, default="exciton_spectrum.dat", help="For the optimized geometry excitonic absorption and circular dichroism spectra are saved to this file [default: %default]")
    parser.add_option("--state", dest="state", type=int, default=0, help="Excitonic state that should be optimized. If no transition charges were provided, only the ground state (0) can be optimized. [default: %default]")
    
    (opts, args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)

    
    ff_file = args[0]
    # read force field definition    
    atomlist, atomtypes, partial_charges, lattice_vectors = read_force_field(ff_file)
    # read transition charge for exciton model (if available)
    if opts.transition_charges != "":
        chromophores = list(read_transition_charges(opts.transition_charges))
    else:
        chromophores = []

    force_field = PeriodicForceField(atomlist, atomtypes, partial_charges, lattice_vectors, chromophores)

    # convert geometry to a vector
    x0 = XYZ.atomlist2vector(atomlist)

    # FIND ENERGY MINIMUM
    # f is the objective function that should be minimized
    # it returns (f(x), f'(x))
    def f(x):
        # also save geometries from line searches
        #save_xyz(x)
        #
        energy, grad = force_field.getEnergyAndGradient(x, state=opts.state)
        print "E = %2.7f     |grad| = %2.7f" % (energy, la.norm(grad))
        #
        return energy, grad

    #xyz_opt = os.path.join("/tmp/", os.path.basename(ff_file).replace(".ff", "_opt.xyz"))
    xyz_opt = ff_file.replace(".ff", "_opt.xyz")
    print "Intermediate geometries will be written to %s" % xyz_opt
    # This is a callback function that is executed by numpy for each optimization step.
    # It appends the current geometry to an xyz-file.
    step=0  # count optimization steps
    def save_xyz(x, mode="a"):
        global step
        if step % opts.save_every != 0:
            #print "skip step %d" % step
            step += 1
            return
        #print "write step %d" % step
        step += 1
        atomlist_opt = XYZ.vector2atomlist(x, atomlist)
        atomlist_opt, dummy = enlarged_unitcell(atomlist_opt, lattice_vectors,
                                                                   nmax=opts.nr_unit_cells)
        XYZ.write_xyz(xyz_opt, [atomlist_opt], \
                      mode=mode)
        
    # save initial geometry
    save_xyz(x0, mode="w")
        
    options = {'gtol': 1.0e-7, 'norm': 2}
    # The "BFGS" method is probably better than "CG", but the line search in BFGS is expensive.
    res = optimize.minimize(f, x0, method="CG", jac=True, callback=save_xyz, options=options)
    # save optimized geometry
    xopt = res.x
    step = 0  # reset optimization counter, so that the last step is written independently of `opts.save_every`
    save_xyz(xopt)
    print "Optimized geometry written to %s" % xyz_opt

    if opts.transition_charges != "":
        # compute exciton spectrum for optimized geometry
        en, T, M = force_field.getTransitionDipoles(verbose=1)
        save_exciton_spectrum(opts.spectrum_file, en, T, M)
    
