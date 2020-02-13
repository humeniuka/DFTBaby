#!/usr/bin/env python
"""
convert an xyz-file to the .gen-format expected by Seifert's DFTB program
Only works for molecules containing H,C,N and O atoms
"""
from DFTB import XYZ
from DFTB import AtomicData

def write_gen(gen_out, atomlist):
    # map Z to DFTB atomtypes, H->4, C->3, N->2, O->1
    atomtypes = {8:1, 7:2, 6:3, 1:4}
    Nat = len(atomlist)
    txt = "  %d C\n" % Nat
    txt +=" O N C H\n"
    for i,(Zi,posi) in enumerate(atomlist):
#        x,y,z = posi
        x,y,z = map(lambda v: AtomicData.bohr_to_angs*v, posi)
        txt += " %d  %d  %20.7f %20.7f %20.7f\n" % (i+1,atomtypes[Zi],x,y,z)
    txt += "%d\n" % (3*Nat)
    for i in range(0, Nat):
        txt += "%d 1\n" % (i+1)
        txt += "%d 2\n" % (i+1)
        txt += "%d 3\n" % (i+1)
    
    fh = open(gen_out, "w")
    print>>fh, txt
    fh.close()

if __name__ == "__main__":
    import sys
    usage = "Usage: python %s <xyz input> <gen output>\n" % sys.argv[0]
    usage += "  converts and xyz-file to the gen-format expected by Seifert's DFTB program\n"

    if len(sys.argv) < 3:
        print usage
        exit(-1)
    xyz_in = sys.argv[1]
    gen_out = sys.argv[2]

    write_gen(gen_out, XYZ.read_xyz(xyz_in)[0])


    
