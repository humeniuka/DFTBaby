#!/usr/bin/env python
"""
plot radial wave functions of pseudo orbitals. You need to check that
the radial grid is large enough so that the tails of all orbitals decay to 
zero.
"""

from numpy import *
from matplotlib.pyplot import *
import os.path
from DFTB import utils
from DFTB import AtomicData

def plot_pseudoatom(atom, pngfile, elmfile=None):
    """
    plot all radial functions associated with a pseudo atom
    and save the plot to a png-file

    Parameters:
    ===========
    atom: atom module, as written by generate_pseudoatoms.py
    pngfile: file name
    elmfile (optional): path to a hotbit .elm file with orbitals,
       corresponding orbitals from the .elm file will be plotted 
       as dotted lines for comparison
    """
    if elmfile != None:
        from hotbit_format import parseHotbitParameters
        hbdata = parseHotbitParameters(elmfile)
    cla()
    title("Z = %s, Nelec = %s" % (atom.Z, atom.Nelec))
    xlabel("r / bohr", fontsize=15)
    print "Atom %s" % AtomicData.atom_names[atom.Z-1]
    print "========"
    for i,orb_name in enumerate(atom.orbital_names):
        print "Orbital %s:  %s  %s" % (orb_name, \
                         ("%2.5f Hartree" % atom.energies[i]).rjust(15), \
                         ("%2.5f eV" % (atom.energies[i]*AtomicData.hartree_to_eV)).rjust(20))
        plot(atom.r, atom.radial_wavefunctions[i], lw=2.5, label="%s en=%s" % (orb_name, atom.energies[i]))
        # compare with hotbit
        if elmfile != None and hbdata.has_key("radial_wavefunction_%s" % orb_name):
            ru = hbdata["radial_wavefunction_%s" % orb_name]
            en = hbdata["orbital_energy_%s" % orb_name]
            plot(ru[:,0], ru[:,1], ls="-.", label="%s en=%s (unconfined)" % (orb_name, en))
    legend()
    show()
#    savefig(pngfile)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: python %s <pseudo atom module>.py [<.elm file for comparison>]" % sys.argv[0]
        print "Plots the atomic orbitals of a pseudo atom."
        exit(-1)
    modfile = sys.argv[1]
    if len(sys.argv) > 2:
        elmfile = sys.argv[2]
    else:
        elmfile = None
    mod = utils.dotdic()
    execfile(modfile, mod)
    plot_pseudoatom(mod, "/tmp/%s.png" % AtomicData.atom_names[mod.Z-1], elmfile)
