"""
This script is a plugin for VMD that facilitates the visualization of cube files 
that are created by DFTBaby: molecular orbitals and transition densities.
VMD has to be compiled with support for Tk, Tcl and Python.
"""

import sys
try:
    import VMD
    import Molecule
    from graphics import color
except ImportError as e:
    print e
    print "This script is a plugin for VMD (Visual Molecular Dynamics)."
    print "It should be imported from the python console of VMD or be called as"
    print "  vmd -python -e %s -args <pattern for cube files>" % sys.argv[0]
    exit(-1)

import os
from os.path import join
import Tkinter, tkFileDialog

class CubeBrowser:
    def __init__(self):
        self.top = Tkinter.Tk()
        self.top.title("Cube Browser")
        self.cube_files = []
    def loadCubes(self, cube_files):
        self.cube_files += cube_files
        self.cube_files.sort()
        # SCROLLBAR
        self.scrollbar = Tkinter.Scrollbar(self.top)
        self.scrollbar.pack(side=Tkinter.RIGHT, fill=Tkinter.Y)
        # LIST OF CUBE FILES
        self.Lb = Tkinter.Listbox(self.top, selectmode=Tkinter.SINGLE)
        # populate listbox with the files
        for i,f in enumerate(cube_files):
            # try to extract the type of the volumetric data, orbital or transition density
            # from the name of the file
            
            # <molecule name>_<type>_<index>.cube
            self.Lb.insert(i,basename(f))
        self.Lb.config(width=0)
        self.Lb.pack()
        self.Lb.bind("<<ListboxSelect>>", self.select)
        # attach listbox to scrollbar
        self.Lb.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.Lb.yview)

        # show the first cubefile in the list
        self.mol = show_cube(cube_files[0])
    def select(self, event):
        """
        This function is called is the users double clicks on a cube file in the
        listbox.
        """
        # index of the selected file
        print self.Lb.curselection()
        isel = int(self.Lb.curselection()[0])
        print "Selected %d-th cube file %s" % (isel, self.cube_files[isel])
        # delete old molecule...
        self.mol.delete()
        # and load new one
        self.mol = show_cube(self.cube_files[isel])

import glob
from os.path import basename

def show_cube(cube_file, iso_value=0.002):
    """
    loads cube file and displays the geometry and the volumetric data in the VMD main window
    """
    mol = Molecule.Molecule()
    mol.load(cube_file)
    mol.rename(cube_file.replace(".cube", ""))
    mol.clearReps()
    # 
    lines = Molecule.MoleculeRep(style="Lines", material="Opaque")
    mol.addRep(lines)
    # isosurface data for positive density
    iso_plus = Molecule.MoleculeRep(style="IsoSurface", material="Transparent")
    iso_plus.changeStyle("IsoSurface %f 0 0 0" % iso_value)
    iso_plus.changeColor("ColorID 0")
    mol.addRep(iso_plus)
    # negative density
    iso_minus = Molecule.MoleculeRep(style="IsoSurface", material="Transparent")
    iso_minus.changeStyle("IsoSurface %f 0 0 0" % (-iso_value))
    iso_minus.changeColor("ColorID 1")
    mol.addRep(iso_minus)

    return mol

if __name__ == "__main__":
    import sys
    browser = CubeBrowser()
    # 
    cube_files = []
    for pattern in sys.argv[1:]:
        cube_files += glob.glob(pattern)
    #
    if len(cube_files) == 0:
        print "No cube files found for patterns '%s'" % "".join(sys.argv[1:])
        exit(-1)

    browser.loadCubes(cube_files)
