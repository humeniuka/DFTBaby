#!/usr/bin/env python
"""
compute electronic forces (without repulsive forces) along a trajectory of geometries using DFTB.
"""
from DFTB import XYZ
from DFTB.PES import PotentialEnergySurfaces

import optparse
import os.path

if __name__ == "__main__":
    import sys
    import os
    usage ="Usage: %s <geometries file> <energies file> <forces file>\n" % os.path.basename(sys.argv[0])
    usage+="  Compute electronic energies and forces using DFTB for all geometries in <geometries file>"
    usage+="  and write them to <energies file> and <forces file>, respectively.\n"
    usage+="  --help shows all options.\n"
    
    parser = optparse.OptionParser(usage)
    parser.add_option("--charge", dest="charge", help="Set total charge. [default: %default]", default=0, type=int)

    (opts,args) = parser.parse_args()

    # delete --charge option from command line, since it would interfere with DFTBaby, which
    # does not understand this option. This is really ugly, but..
    for i in range(0, len(sys.argv)):
        if "--charge=" in sys.argv[i]:
            del sys.argv[i]
            
    if len(args) < 3:
        print usage
        exit(-1)
        
    geom_file = args[0]
    energy_file = args[1]
    force_file = args[2]

    print "Compute forces with DFTB"
    print "========================"
    print ""
    fh_en = open(energy_file, "w")
    print>>fh_en, "# ELECTRONIC ENERGY / HARTREE"
    
    # first geometry
    atomlist = XYZ.read_xyz(geom_file)[0]
    # read charge from title line in .xyz file
    kwds = XYZ.extract_keywords_xyz(geom_file)
    charge = kwds.get("charge", opts.charge)
    
    pes = PotentialEnergySurfaces(atomlist, charge=charge)
    # dftbaby needs one excited states calculation to set all variables
    x = XYZ.atomlist2vector(atomlist)
    pes.getEnergies(x)
    
    for i,atomlist in enumerate(XYZ.read_xyz(geom_file)):
        # compute electronic ground state forces with DFTB
        x = XYZ.atomlist2vector(atomlist)
        en = pes.getEnergy_S0(x)
        # total ground state energy including repulsive potential
        en_tot = en[0]
        print "Structure %d   enTot= %s Hartree" % (i, en_tot)
        # electronic energy without repulsive potential
        en_elec = pes.tddftb.dftb2.getEnergies()[2]
        gradVrep, gradE0, gradExc = pes.grads.gradient(I=0)
        # exclude repulsive potential from forces
        grad = gradE0 + gradExc
        forces = XYZ.vector2atomlist(-grad, atomlist)
        if i == 0:
            mode = 'w'
        else:
            mode = 'a'
        print>>fh_en, "%2.7f" % en_elec
        
        XYZ.write_xyz(force_file, [forces],
                      title="forces with DFTB, total_energy=%2.7f  elec_energy=%2.7f" % (en_tot, en_elec),
                      units="hartree/bohr", mode=mode)

    fh_en.close()
