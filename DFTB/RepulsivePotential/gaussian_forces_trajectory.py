#!/usr/bin/env python
"""
compute forces along a trajectory of geometries using Gaussian.
"""
from DFTB import XYZ
import Gaussian
import getpass
from optparse import OptionParser
import os.path

if __name__ == "__main__":
    import sys
    import os
    usage = "Usage: %s <geometries file> <energies file> <forces file>\n\n" % os.path.basename(sys.argv[0])
    usage += "Compute total energies and forces for all geometries in <geometries file>"
    usage += "and write them to <energies file> and to <forces file>, respectively.\n"
    usage += "Before using this program, load the Gaussian environment: 'module load g09'"
    if len(sys.argv) < 4:
        print usage
        exit(-1)
    parser = OptionParser(usage=usage)
    parser.add_option("--method", dest="method", help="quantum chemical method [default: %default]", default="LC-PBEPBE")
    parser.add_option("--basis_set", dest="basis_set", help="name of basis set [default: %default]", default="6-311+G*")
    parser.add_option("--charge", dest="charge", type="int", help="total charge [default: %default]", default=0)
    parser.add_option("--multiplicity", dest="multiplicity", type="int", help="multiplicity [default: %default]", default=1)

    (options,args) = parser.parse_args()
    geom_file = args[0]
    energy_file = args[1]
    force_file = args[2]
    tmp_dir = "./g09_temp/"
    os.system("mkdir -p %s" % tmp_dir)
    mode = 'w'
    print "Compute forces with GAUSSIAN"
    print "============================"
    print "Method:     %s" % options.method
    print "Basis set:  %s" % options.basis_set
    print ""
    
    fh_en = open(energy_file, "w")
    print>>fh_en, "# TOTAL ENERGY / HARTREE     %s/%s" % (options.method, options.basis_set)

    
    kwds = XYZ.extract_keywords_xyz(geom_file)
    if kwds.has_key("charge"):
        # charge keyword in file overrides command line option --charge
        options.charge = int(kwds["charge"])
    
    for i,atomlist in enumerate(XYZ.read_xyz(geom_file)):
        tmp_com = tmp_dir + "gaussian.com"
        tmp_out = tmp_dir + "gaussian.log"
        tmp_chk = tmp_dir + "gaussian.chk"
        if i == 0:
            guess=""
        else:
            # read starting MOs from checkpoint file
            guess="Guess=Read"
        Gaussian.write_input(tmp_com, atomlist, \
                        route="# %s/%s %s Force Integral=(Grid=Ultrafine) Geom=NoCrowd" % (options.method, options.basis_set, guess), \
                        chk_file=tmp_chk, \
                        title="Gradients for fitting repulsive potential", \
                        charge=options.charge, multiplicity=options.multiplicity)
        # compute forces with gaussian
        os.system("g09 < %s > %s" % (tmp_com, tmp_out))

        en_tot = Gaussian.read_total_energy(tmp_out)
        print "Structure %d   enTot= %s Hartree" % (i, en_tot)
        forces_dic = Gaussian.read_forces(tmp_out)
        forces = []
        for center in forces_dic.keys():
            forces.append(forces_dic[center])
        print>>fh_en, "%2.7f" % en_tot
        XYZ.write_xyz(force_file, [forces],
                      title="forces with %s/%s, total_energy=%2.7f" % (options.method, options.basis_set, en_tot),
                      units="hartree/bohr",
                      mode=mode)
        mode = 'a'

    fh_en.close()
    
