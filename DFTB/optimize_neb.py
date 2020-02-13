#!/usr/bin/env python
"""
Optimize a reaction path using the nudged elastic band algorithm
"""
from DFTB.PES import PotentialEnergySurfaces
from DFTB import NEB
from DFTB import XYZ, utils, AtomicData
from DFTB.Timer import GlobalTimer as T

import numpy as np
import numpy.linalg as la

class PESforNEB(PotentialEnergySurfaces):
    name = ""
    def getEnergiesAndGradient(self, x, I):
        #
        if I == 0 and type(self.tddftb.XmY) != type(None):
            # only ground state is needed. However, at the start
            # a single TD-DFT calculation is performed to initialize
            # all variables (e.g. X-Y), so that the program does not
            # complain about non-existing variables.
            enI, gradI = super(PESforNEB, self).getEnergyAndGradient_S0(x)
            energies = np.array([enI])
        else:
            energies, gradI = super(PESforNEB, self).getEnergiesAndGradient(x, I)
            enI = energies[I]
        #
        return energies, gradI

    def setName(self, name):
        # The name is added to the file names of the intermediate output
        self.name = name
    # callback function for saving intermediate reaction paths
    def plot(self, images=[], energies=[], istep=0):
        geometries = [XYZ.vector2atomlist(im, self.atomlist) for im in images]
        xyz_out = "neb_%s_%d.xyz" % (self.name, istep)
        XYZ.write_xyz(xyz_out, geometries)
        print "wrote path for iteration %d to %s" % (istep, xyz_out)
        if energies != []:
            # first column: index of image
            # second column: energy of image
            data = np.vstack((np.arange(0, len(images)), energies)).transpose()
            np.savetxt("path_energies_%s_%d.dat" % (self.name, istep), data)

if __name__ == "__main__":
    import sys
    from scipy import optimize
    from os.path import basename

    if len(sys.argv) < 2:
        print "Usage: %s <xyz-file with educt, intermediates and product> I" % basename(sys.argv[0])
        print "  optimize reaction path on the I-th electronic state"
        print "  type --help to see all options"
        print "  to reduce the amount of output add the option --verbose=0"
        exit(-1)
    
    #
    xyz_file = sys.argv[1]    # path to xyz-file
    I = int(sys.argv[2])      # index of electronic state
    # Read the geometry from the xyz-file
    atomlists = XYZ.read_xyz(xyz_file)
    atomlist = atomlists[0]
    # read the charge of the molecule from the comment line in the xyz-file
    kwds = XYZ.extract_keywords_xyz(xyz_file)
    # initialize the TD-DFTB calculator
    pes = PESforNEB(atomlist, Nst=I+2, **kwds)
    pes.setName(basename(xyz_file.replace(".xyz", "")))

    # convert geometry to a vector
    x0 = XYZ.atomlist2vector(atomlist)

    neb = NEB.NEB(force_constant=1.0, mass=1.0)
    neb.setEnergyCalculator(pes)
    images = [XYZ.atomlist2vector(atomlist) for atomlist in atomlists]
    neb.setImages(images, states=[I for im in images])
    neb.addImagesLinearly(2)
    neb.plot() # save initial path
    neb.findMEP(tolerance=0.02, nsteps=1000, dt=0.1, friction=0.2, optimize_endpoints=True)

    me = neb.splineMEProfile()
    mep = neb.splineMEPath()
    

    # timing
    print T
    
