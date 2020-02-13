"""
extract data blocks from a hotbit parameter .par file or a hotbit element .elm file
"""
import numpy
import re
import DFTB.AtomicData as AtomicData

class _Block(object):
    """base class, all derived blocks have to define their own __init__, process and close methods."""
    def __init__(self, args=[]):
        """called at the beginning of a block"""
        pass
    def name(self):
        """the block is saved to the data dictionary under this name"""
        return "Block Name"
    def process(self, line):
        pass
    def close(self):
        """called at the end of a block"""
        pass
    def getData(self):
        """returns a representation of the block data"""
        pass

class _IgnoreBlock(_Block):
    pass

def parseHotbitParameters(filename, convert=True):
    """
    parseHotbitParameters(filename)  -  extract some data blocks from a 
    hotbit parameter (.par) files or hotbit element (.elm) files

    The following blocks can be processed:
    - *_*_table
    - repulsion

    Returns:
    ========
    a dictonary of data blocks
    """
    Data = {}
    """A dictonary containing all the blocks read."""
    fh = open(filename, 'r')
    block = _IgnoreBlock()
    for line in fh.readlines():
        if (line[0] == "#"):
            """comment"""
            continue
        elif "=" in line:
            """A label followed by = starts a new data block. Close the previous block and 
            open a new one. The name of the block is used to find the correct derived class 
            (by appending _Block). E.g. if the block "repulsion" is encountered an instance 
            of class repulsion_Block is created, its "__init__"-function is called at the 
            beginning of the block with the tokens on the first line  as a list of arguments, 
            the method "process" is called on each line within the block and the method "close" 
            when another block is encountered or at the end of the file."""
            block.close()
            """close previous block"""
            key,val = line.strip().split("=")
            if not re.match("^[\w-]+$", key.strip()):
                key = "invalid"
            args = key.split("_")
            if ".par" in filename:
                block_type = args[-1]
                args = args[:-1]
            elif ".elm" in filename:
                block_type = args[0]
                args = args[1:] + [val]
            try:
                block = eval("_" + block_type + "_Block")(args)
                """call the constructor of a class that bears the same name as the current block"""
                Data[block.name()] = block
            except NameError:
                block = _IgnoreBlock()
        else:
            block.process(line.strip())
            """The process-method of the current block is used."""
    block.close()
    fh.close()
    PyData = {}
    if convert == True:
        """convert blocks into format suitable for PyQuante"""
        for (block_name, block) in Data.iteritems():
            PyData[block_name] = block.getData()
        return PyData
    return Data

########## blocks that can occur in .par files ##########

class _table_Block(_Block):
    def __init__(self, args=[]):
        self.atom1 = args[0]
        self.atom2 = args[1]
        self.data = []
    def name(self):
        return "slako_integrals_%s_%s" % (self.atom1, self.atom2)
    def process(self, line):
        if line != "":
            self.data.append( numpy.fromstring(line, sep=" ") )
    def close(self):
        # there are 1 column for distances, 10 columns for H integrals and 10 columns for S integrals
        self.data = numpy.reshape(numpy.hstack(self.data), (len(self.data), 21))
    def getData(self):
        return self.data

class _repulsion_Block(_Block):
    def __init__(self, args=[]):
        self.data = []
    def name(self):
        return "repulsive_potential"
    def process(self, line):
        if line != "":
            self.data.append( numpy.fromstring(line, sep=" ") )
    def close(self):
        self.data = numpy.reshape(numpy.hstack(self.data), (len(self.data), 2))
    def getData(self):
        return self.data

####### blocks that can occur in .elm files #############

class _u_Block(_repulsion_Block): 
    """radial wave functions"""
    def __init__(self, args=[]):
        self.orbname = args[0].strip()
        self.data = []
    def name(self):
        return "radial_wavefunction_%s" % self.orbname

class _epsilon_Block(_Block):
    """orbital energy"""
    def __init__(self, args=[]):
        self.orbname = args[0].strip()
        self.energy = float(args[1])
    def name(self):
        return "orbital_energy_%s" % self.orbname
    def getData(self):
        return self.energy

class _U_Block(_Block):
    """Hubbard U parameter"""
    def __init__(self, args=[]):
        self.U = float(args[0].strip())
    def name(self):
        return "hubbard_U"
    def getData(self):
        return self.U

class _symbol_Block(_Block):
    """atom name"""
    def __init__(self, args=[]):
        self.atom_name = args[0].strip()
    def name(self):
        return "atom_name"
    def getData(self):
        return self.atom_name.lower()

def load_hotbit_pseudoatom(elmfile):
    """
    load description of pseudo atom from a hotbit element (.elm) file

    Parameters:
    ===========
    elmfile: path to element description

    Returns:
    ========
    atom module with data stored as members
    """
    Data = parseHotbitParameters(elmfile)
    class Atom: pass
    atom = Atom()
    atom.Z = AtomicData.atomic_number(Data["atom_name"])
    atom.energies = []
    atom.valence_orbitals = []
    atom.nshell = []
    atom.angular_momenta = []
    atom.orbital_occupation = []
    atom.r = None
    atom.radial_wavefunctions = []
    spectr2l = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5, "i": 6}
    for k in Data.keys():
        if "orbital_energy" in k:
            atom.energies.append(float(Data[k]))
            dummy,spectr = k.rsplit("_",1)
            n = int(spectr[0])
            l = int(spectr2l[spectr[1]])
            occnr = AtomicData.valence_occ_numbers[Data["atom_name"]]
            if occnr.has_key(spectr):
                atom.orbital_occupation.append(occnr[spectr])
            else:
                atom.orbital_occupation.append(0)
            atom.valence_orbitals.append(len(atom.valence_orbitals))
            atom.nshell.append(n)
            atom.angular_momenta.append(l)
            if Data.has_key("radial_wavefunction_%s" % spectr):
                atom.r = Data["radial_wavefunction_%s" % spectr][:,0]
                atom.radial_wavefunctions.append(Data["radial_wavefunction_%s" % spectr][:,1])
            atom.hubbard_U = Data["hubbard_U"]
    return atom

if __name__ == "__main__":
    import sys
    parfile = sys.argv[1]
    Data = parseHotbitParameters(parfile)
    print Data.keys()
    if ".elm" in parfile:
        print load_hotbit_pseudoatom(parfile).__dict__
