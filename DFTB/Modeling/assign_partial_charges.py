import numpy as np
import sys

if len(sys.argv) < 4:
    print "Usage: %s <solvent box .xyz> <partial charges of solv. molecule> <output file>" % sys.argv[0]
    print "  Assigns partial charges to each solvent molecule in solvent box."
    exit(-1)

xyz_in = sys.argv[1]
charges_in = sys.argv[2]
charges_out = sys.argv[3]

# partial charges for each atom in a solvent molecule
charges = np.loadtxt(charges_in)

# solvent box
fh_in = open(xyz_in, "r")
lines = fh_in.readlines()
fh_in.close()

fh_out = open(charges_out, "w")

solvent = lines[2:]

# replace atomic labels in xyz-file with partial charges 
iat=0
nsolv = len(charges) # number of atoms in single solvent molecule
header = lines[:2]
for hl in header:
    print>>fh_out, hl.strip()

for l in solvent:
    iat = (iat + 1) % len(charges)
    parts = l.strip().split()
    parts[0] = "%s" % charges[iat]
    lrepl = " ".join(parts)
    print>>fh_out, lrepl

fh_out.close()
    
