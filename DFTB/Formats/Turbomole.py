import re
import numpy as np
from numpy import sqrt
from DFTB.AtomicData import atom_names, atomic_number

class _Block(object):
    """base class, all derived blocks have to define their own __init__, process and close methods."""
    def __init__(self, args=[]):
        """called at the beginning of a block"""
        pass
    def process(self, line):
        pass
    def close(self):
        """called at the end of a block"""
        pass
    def getPyQ(self):
        """returns a representation of the block that can be used with PyQuante"""
        return(None)
    def tm_repr(self):
        """convert the block back to a representation that is understood by turbomole"""
        return ""

class _IgnoreBlock(_Block):
    pass

_coef_pattern = re.compile(r"""[-\d]            # digit or minus sign before colon
                              \.               # colon
                              \d{13,14}        # 13 or 14 digits after colon
                              D[+-]\d{2}       # exponent
                                      """, re.VERBOSE)
_eigval_pattern = re.compile(r"""\s*(\d+)       # number of the eigenvector
                                \s+
                                eigenvalue\s*=\s*
                                ([+-]*\d*.\d+D[+-]\d{2})  # eigenenergy with format e.g. 0.1493787507771279D+00
                                      """, re.VERBOSE)

class _eigenpairs_Block(_Block):
    """reads eigenvector file"""
    def __init__(self, args=[]):
        self.nstates = 0
        self.energies = []
        self.vecs = []
        self.nstates = 0
        self.eigval_pattern = _eigval_pattern
    def process(self, line):
        result = self.eigval_pattern.match(line)
        if (result == None):
            coefs = _coef_pattern.findall(line)
            for cstr in coefs:
                c = float(cstr.replace('D', 'E'))
                self.vecs[self.nstates-1].append(c)
        else:
            self.nstates += 1
            en = float(result.group(2).replace('D', 'E'))
            self.vecs.append([])
            self.energies.append(en)
    def close(self):
        pass
    def getPyQ(self):
        en = np.array(self.energies)
        eigvecs = np.array(self.vecs).transpose()
        return( (en, eigvecs) )
    def __str__(self):
        text = "eigenpairs-Block"
        for st in xrange(self.nstates):
            text += "Energy = %s\n" % self.energies[st]
            text += "Coefs = %s\n" % (self.vecs[st][0:1] + ["..."] + self.vecs[st][-2:])
        return(text)

class _uhf_Block(_Block):
    pass

def order_excs(nocc, nvirt, maxvirts=10e10):
    """
    order_excs(nocc, nvirt)  -  create ordered list of excitations
    which allows to map the components of a CIS eigenvector to the corresponding
    single excitations (spin,o,v) where the occupied orbital o is replaced by the virtual
    v both with the same spin (0 := spin-up, 1 := spin-down).

    nocc    tuple (nalpha, nbeta) with number of occupied orbitals in reference 
            for spin-up and spin-down
    nvirt   tuple (nvirta, nvirtb) with number of virtual orbitals for spin-up and spin-down 
    
    Option   Default   Description
    --------------------------------
    maxvirt   10e10    Only excitations up to at most the maxvirt-th spin orbital are included
    
    return  list of lists with (spin,o,v) tuples, excs = [[(spin1,o1,v2)],[(spin2,o2,v2)],...]
            If cis is an eigenvector in the basis of single excited determinants,
            then cis[i] is the component in the direction of the excs[i] excitation
    """
    excs_a = []
    """excitations with spin-up"""
    for o in xrange(nocc[0]):
        for v in xrange(nocc[0], min(maxvirts, nocc[0]+nvirt[0])):
            excs_a.append( [(0, o, v)] )
    excs_b = []
    """excitations with spin-down"""
    for o in xrange(nocc[1]):
        for v in xrange(nocc[1], min(maxvirts, nocc[1]+nvirt[1])):
            excs_b.append( [(1, o, v)] )
    excitations = excs_a + excs_b
    return(excitations)

_moeigval_pattern = re.compile(r"""\s*(\d+)               # number of the eigenvector
                                \s+
                                \S+\s+                     # symmetry label
                                eigenvalue\s*=\s*
                                ([+-]*\d*.\d+D[+-]\d{2})  # eigenenergy with format e.g. 0.1493787507771279D+00
                                \s+nsaos\s*=\s*
                                (\d+)         # number of basis function
                                      """, re.VERBOSE)

class _scfmo_Block(_eigenpairs_Block):
    """read molecular orbital coefficients in a Turbomole mo-file"""
    def __init__(self, args=[]):
        super(_scfmo_Block, self).__init__()
        """the mo blocks have the same structre as sing_a files except
        that the number of basis functions and the symmetry label is 
        specified on the first line."""
        self.eigval_pattern = _moeigval_pattern

class _uhfmo_alpha_Block(_scfmo_Block):
    """spin-up orbitals from unresticted HF"""
    pass

class _uhfmo_beta_Block(_scfmo_Block):
    """spin-down orbitals from unresticted HF"""
    pass

class _basis_Block(_Block):
    def __init__(self, args=[]):
        self.basis_set = ""
        """name of basis set"""
        self.basis_data = {}
        self.prims = []
        self.found_basis = False
        self.atno = None
        """atomic number of the element to which the basis belongs"""
        self.ncgbf = 0  
        """number of primitive gaussians in the current contraction"""
        self.sym = ""   
        """angular momentum labels 'S', 'P', 'D', or 'F'"""
    def process(self, line):
        """an element name in the first column indicates the start of a basis set"""
        if (line[0].isalpha()):
            self.found_basis = True
            self.ncgbf = 0
            (atom, self.basis_set) = line.strip().split(' ', 1)
            self.atno = atomic_number(atom)
            self.basis_data[self.atno] = []
        elif (line[0] == ' '):
            if (self.found_basis):
                (a, b) = line.strip().split()
                if (self.ncgbf == 0):
                    """read number of contractions and angular momentum, e.g.  5 s"""
                    self.ncgbf = int(a)
                    self.sym = b.strip().upper()    
                    self.prims = []
                else:
                    """read self.ncgbf exponent and coefficient pairs"""
                    exp, coef = float(a), float(b)
                    self.prims.append( (exp, coef) )
                    self.ncgbf -= 1
                    if (self.ncgbf == 0):
                        """no more coefficients to read, save contraction"""
                        self.basis_data[self.atno].append( (self.sym, self.prims) )
                        
                    """decrease the number of coefficients which are still to be read"""
    def getPyQ(self):
        """In Turbomole basis file the shells are not alway sorted in the order s,p,d,f"""
        symorder = {'S': 0, 'P': 1, 'D': 2, 'F': 3}
        """order of orbitals, first s then p, d, f"""
        for atno,contr in self.basis_data.iteritems():
            contr.sort(cmp=lambda x,y: symorder[x[0]]-symorder[y[0]])
        return(self.basis_data)
    def __str__(self):
        text = "Basis Set %s\n" % self.basis_set
        text += str(self.basis_data)
        return(text)

class _coord_Block(_Block):
    def __init__(self, args=[]):
        self.atomlist = []
    def process(self, line):
        parts = line.strip().split()
        x,y,z = float(parts[0]), float(parts[1]), float(parts[2])
        atno = atomic_number(parts[3])
        self.atomlist.append( (atno, (x,y,z)) )
    def getPyQ(self):
        return(self.atomlist)
    def __str__(self):
        return(str(self.atomlist))

class _point_charges_Block(_Block):
    def __init__(self, args=[]):
        self.point_charges = []
    def process(self, line):
        parts = line.strip().split()
        x,y,z, charge = map(float, parts)
        self.point_charges.append( (charge, (x,y,z)) )
    def getPyQ(self):
        return(self.point_charges)
    def __str__(self):
        return(str(self.point_charges))
        

class _energy_Block(_Block):
    def __init__(self, args=[]):
        self.iter_energies = []
    def process(self, line):
        p = line.strip().split()
        it,scf,scfkin,scfpot = int(p[0]), float(p[1]), float(p[2]), float(p[3])
        self.iter_energies.append( (scf, scfkin, scfpot) )
    def getPyQ(self):
        try:
            scf,scfkin,scfpot = self.iter_energies[-1]
            """return converged total energy"""
        except:
            scf = None
        return(scf)

class _closed_Block(_Block):
    """find number of closed shells"""
    def __init__(self, args=[]):
        self.nocc = 0
    def process(self, line):
        p = line.strip().split()
        (symlabel, states, occ) = p[0], p[1].split('-'), p[2]
        if (len(states) == 2):
            """a range is given, e.g.  1-4"""
            nstates = int(states[1]) - int(states[0]) + 1
            self.nocc += nstates
        else:
            """only one state"""
            self.nocc += 1
    def getPyQ(self):
        return(self.nocc)

class _alpha_Block(_closed_Block):
    """find number of spin-up shells"""
    pass

class _beta_Block(_closed_Block):
    """find number of spin-down shells"""
    pass


class _grad_Block(_Block):
    def __init__(self, args=[]):
        self.atomlist = []
        self.coordvector = []
        self.gradient = []
    def process(self, line):
        if "cycle" in line:
            """rest, we only want the converged gradient"""
            self.atomlist = []
            self.coordvector = []
            self.gradient = []
            return
        parts = line.strip().split()
        if len(parts) == 4:
            x,y,z = float(parts[0].replace("D","E")), float(parts[1].replace("D","E")), float(parts[2].replace("D","E"))
            atno = atomic_number(parts[3])
            self.atomlist.append( (atno, (x,y,z)) )
            self.coordvector += [x,y,z]
        elif len(parts) == 3:
            gradx,grady,gradz = float(parts[0].replace("D","E")), float(parts[1].replace("D","E")), float(parts[2].replace("D","E"))
            self.gradient += [gradx,grady,gradz]
    def getPyQ(self):
        return(self.atomlist, np.array(self.coordvector), np.array(self.gradient))
    def __str__(self):
        return(str(self.atomlist))

class _couplingvector_Block(_Block):
    def __init__(self, args=[]):
        self.coupling_vector = []
    def process(self, line):
        parts = line.strip().split()
        self.coupling_vector += [float(parts[0].replace("D","E")), float(parts[1].replace("D","E")), float(parts[2].replace("D","E"))]
    def getPyQ(self):
        return np.array(self.coupling_vector)
        
class _actual_Block(_Block):
    def __init__(self, args=[]):
        self.failed_step = args[1]
    def getPyQ(self):
        return self.failed_step

class _HashBlock(_Block):
    """Each line in the block consists of a key-value-pair separated by whitespaces"""
    def __init__(self, args=[]):
        self.hash = {}
    def process(self, line):
        (key, value) = line.split(None,1)
        self.hash[key.strip()] = value.strip()
    def close(self):
        """called at the end of a block"""
        pass
    def getPyQ(self):
        """returns a representation of the block that can be used with PyQuante"""
        return(self.hash)
    def tm_repr(self):
        block_type = self.__class__.__name__[1:].replace("_Block","")
        txt = "$%s\n" % block_type
        for k,v in self.hash.iteritems():
            txt += "\t%s %s\n" % (k,v)
        return txt

class _dft_Block(_HashBlock):
    pass

"""
In Turbomole pure d-, f- and g-orbitals are used, while PyQuante uses cartesian basis functions.
Ordering in PyQuante of cart. d-orbitals: (2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)
Ordering in PyQuante of cart. f-orbitals: (3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
           (0,3,0),(0,2,1),(0,1,2), (0,0,3)
Ordering in Turbomole of pure d-orbitals: d0, d1a, d1b, d2a, d2b
Ordering in Turbomole of pure f-orbitals: f0, f1a, f1b, f2a, f2b, f3a, f3b

Also there are some knacks with the normalization, that need be taken care of.
The matrix _TM_U transforms pure basis functions into properly 
normalized cartesian basis functions in the usual order.
"""

_TM_U = {
    'D': np.array([
            [-0.5    , 0 , 0 , 0          , sqrt(3.0)/2.0],
            [0      , 0 , 0 , 1          , 0            ],
            [0      , 1 , 0 , 0          , 0            ],
            [-0.5   , 0 , 0 , 0          ,-sqrt(3.0)/2.0],
            [0      , 0 , 1 , 0          , 0            ],
            [1      , 0 , 0 , 0          , 0            ]]
               ),
    'F': np.array([
            [        0           ,   -sqrt(15.0/40.0)  ,        0          ,       0      ,        0        , sqrt(15.0/24.0)    ,       0              ],
            [        0           ,          0          , -sqrt(3.0/40.0)   ,       0      ,        0        ,         0          , -3.0*sqrt(3.0/24.0)  ],
            [-3.0*sqrt(3.0/60.0) ,          0          ,        0          ,       0      ,   sqrt(3.0)/2.0 ,         0          ,       0              ],
            [        0           , -sqrt(3.0/40.0)     ,        0          ,       0      ,        0        , -3.0*sqrt(3.0/24.0),       0              ],
            [        0           ,          0          ,        0          ,       1      ,        0        ,         0          ,       0              ],
            [        0           , 4.0*sqrt(3.0/40.0)  ,        0          ,       0      ,        0        ,         0          ,       0              ],
            [        0           ,          0          , -sqrt(15.0/40.0)  ,       0      ,        0        ,         0          , sqrt(15.0/24.0)      ],
            [ -3.0*sqrt(3.0/60.0),          0          ,        0          ,       0      ,  -sqrt(3.0)/2.0 ,         0          ,       0              ],
            [        0           ,          0          , 4.0*sqrt(3.0/40.0),       0      ,        0        ,         0          ,       0              ],
            [ 2.0*sqrt(15.0/60.0),          0          ,        0          ,       0      ,        0        ,         0          ,       0              ]
            ])
    }


_ang2sym = ['S', 'P', 'D', 'F']
"""convert angular momentum numbers 0,1,2,3 to spectroscopic symbols"""
_orbdim_pure = [1, 3, 5, 7]
"""number of basis functions in pure S, P, D and F shells"""
_orbdim_cart = [1, 3, 6, 10]
"""number of basis functions in cartesian S,P,D and F shells"""

def count_bfs(bs):
    """count_bfs(bs)  - count number of basis functions in pure and cartesian basis
    The basis set bs only contains cartesian functions, len(bs.bfs) of them.
    However count_bfs deduces from the angular momentum how many functions there should
    be if these cartesian functions are replaced by spherical ones.
    
    bs  -  BasisSet object
    
    return   tuple (npure, ncart) with the # of pure and the # of cartesian functions
    """
    nbfs = len(bs.bfs)
    np = 0
    nc = 0
    while(nc < nbfs):
        bf = bs.bfs[nc]
        L = bf.ang_mom
        np += _orbdim_pure[L]
        nc += _orbdim_cart[L]
    assert(nc == nbfs)
    return( (np, nc) )

def canonical_mos(bs, TMorbs):
    """
    canonical_mos(bs, TMorbs) - convert Turbomole SAOs to canonical molecular orbitals
                                and reorder mo coefficients for d- and f-orbitals
    
    bs    - BasisSet object
    TMorbs  - molecular orbitals in Turbomole format

    return -  transformed, canonical orbitals in cartesian basis

    The molecular orbital coefficients generated by Turbomole are not canonical.
    Also the d- and f-orbitals are ordered differently in Turbomole. 
    To obtain canonical properly ordered coefficients, the mos belonging to d-functions 
    have to be multiplied by the 6x5 matrix TM_U["D"] 
    (assuming that 5d/7f cartesian basis functions are used in Turbomole), 
    the mos belonging to f-functions have to be multiplied by the 10x7 matrix TM_U["F"].
    """
    (Np, Nc) = count_bfs(bs)
    orbs_can = np.zeros( (Nc, Np) )
    """holds transformed orbitals"""
    if (TMorbs.shape[1] == Nc) and (Np != Nc):
        raise Exception("Turbomole MOs apparently are in cartesian basis (6d/10f), please use pure basis functions (5d/7f) instead.")
    assert( Np == TMorbs.shape[1] )
    """make sure mos are in pure basis"""
    np = 0
    """count basis functions in pure basis"""
    nc = 0
    """count basis functions in cartesian basis"""
    while(nc < Nc):
        bf = bs.bfs[nc]
        L = bf.ang_mom
        if L in [2, 3]:
            dimp = _orbdim_pure[L]
            dimc = _orbdim_cart[L]
            orbs_can[nc:nc+dimc,:] = np.dot(_TM_U[_ang2sym[L]], TMorbs[np:np+dimp,:])
            """transform portion of coefficient matrix belonging to this shell"""
            np += dimp
            nc += dimc
            """skip other basis functions belonging to the same shell"""
        elif L > 3:
            raise("transformation for G orbitals not implemented yet")
        else:
            orbs_can[nc,:] = TMorbs[np,:]
            np += 1
            nc += 1
            """basis functions of S and P orbitals don't require transformation"""
    assert(np == Np)
    assert(nc == Nc)
    Nspurious = Nc-Np
    """Nspur is the difference in the number of basis functions if the cartesian
    basis is used as opposed to the cartesian basis which always results in less
    basis functions. If Nspur is not zero the number of molecular orbitals does
    not equal the number of basis functions any more and orbs is not a square matrix."""
#    print "Nspurious = %s !!!" % Nspurious
    return(orbs_can)

def parseTurbomole(filename, convert=True):
    """
    parseTurbomole(filename)  -  extract some data blocks from a Turbomole file
    The following files can be processed:
    - coord
    - basis
    - mos
    - sing_a, ucis_a
    - energy
    return   -   a dictonary of data blocks
    -Data["coord"] contains a structure that allows to construct a Molecule object
    in PyQuante,   mol = Molecule('', atomlist=Data["coords"])
    -Data["basis"] contains a structure that can be used to initialize a BasisSet
    object in PyQuante,    bs = BasisSet(atoms, basis_data=Data["basis"])
    -Data["scfmo"] contains a tuple with the orbital energies and the molecular
    orbital coefficients,  (orbe, orbs) = Data["scfmo"]
    (orbe is an array, orbs is an array where orbs[:,i] contains the i-th eigenvector)
    use canonical_mos(bs, orbs) to order and normalize the orbitals properly.
    -Data["eigenpairs"] contains a tuple with the excitation energies and the 
    CIS coefficients, to be able to assign the vectors to particle-hole excitations
    (spin,o,v) the number of occupied and virtual orbitals and the number of
    spin-up and spin-down electrons need to be known.
    (Ecis,CIS) = Data["eigenpairs"], CIS[:,i] contains the eigenvector in the basis of single
    excitations, CIS[:nalpha,i] holds the spin-up excitations and CIS[nalpha:,i] the spin-down
    excitations. 
    -Data["energy"] contains the converged total ground state energy from an SCF calculation
    """
    Data = {}
    """A dictonary containing all the blocks read."""
    fh = open(filename, 'r')
    block = _IgnoreBlock()
    for line in fh.readlines():
        if (line[0] == "#"):
            """comment"""
            continue
        elif (line[0] == "$"):
            """A dollar sign starts a block of information. Close the previous block and 
            open a new one. The name of the block is used to find the correct derived class 
            (by appending _Block). E.g. if the block "spectrum" is encountered an instance 
            of class spectrum_Block is created, its "__init__"-function is called at the 
            beginning of the block with the tokens on the first line  as a list of arguments, 
            the method "process" is called on each line within the block and the method "close" 
            when another block is encountered or at the end of the file."""
            args = line[1:].split()
            block_type = args[0]
            try:
                block.close()
                """close previous block"""
                block = eval("_" + block_type + "_Block")(args[1:])
                """call the constructor of a class that bears the same name as the current block"""
                Data[block_type] = block
            except:
                block = _IgnoreBlock()
        else:
            block.process(line)
            """The process-method of the current block is used."""
    block.close()
    fh.close()
    PyData = {}
    if convert == True:
        """convert blocks into format suitable for PyQuante"""
        for (block_type, block) in Data.iteritems():
            PyData[block_type] = block.getPyQ()
        return(PyData)
    else:
        return Data

def importTurbomole(tm_dir):
    """importTurbomole(tm_dir) -  import coordinates, basis, energies, mo and cis coefficients from Turbomole
    tm_dir     -   path to directory with output files from converged ground state calculation
    and excited state calculation with dft and cis
    return  tuple (mol,bs,en,nocc,nvirt,orbe,orbs,Ecis,UCIS,excitations)
    mol - Molecule object
    bs  - BasisSet object
    en  - ground state energy
    nocc - tuple (nocca, noccb) with occupied spin-up and occupied spin-down orbitals
    nivrt - tuple (nvirta, nvirtb) with virtual orbitals for spin-up and spin-down
    orbe  - tuple (orbea, orbeb) orbital energies for spin-up and spin-down
    orbs  - tuple (orbsa, orbsb) spin-up and spin-down mo coefficients
    Ecis  - excited state energies
    UCIS  - cis coefficients
    excitations  -  list of excitations to map between components in a cis vector UCIS[:,st] 
    and particle-hole pairs (spin, o, v)

    tm_dir has to contain the following files:
    for neutral molecule:   coord, basis, control, mos, sing_a, energy
    for ionized molecule:   coord, basis, control, alpha, beta, unrs_a, energy
    """

    from PyQuante.Molecule import Molecule
    from PyQuante.Basis.basis import BasisSet
    # read coord file
    Data = parseTurbomole(tm_dir + "/coord")
    mol = Molecule('', atomlist = Data["coord"])
    # read basis
    Data = parseTurbomole(tm_dir + "/basis")
    bs = BasisSet(mol, basis_data = Data["basis"])
    # check for unresticted calculation and determine functional
    Data = parseTurbomole(tm_dir + "/control")
    if Data.has_key("dft"):
        functional = Data["dft"]["functional"].lower()
    else:
        raise Exception("Can only import data from DFT calculations.")
    if (not Data.has_key("uhf")):
        """rhf data. In an rhf calculations molecular orbitals
        and CIS coefficients for alpha and beta electrons are the same.
        In order to treat uhf and rhf data with the same functions
        duplicate orbitals and cis coefficients."""
        # find occupation
        (nbfsp, nbfsc) = count_bfs(bs)
        nocca = Data["closed"]
        noccb = nocca
        nvirta = nbfsp-nocca
        nvirtb = nbfsp-noccb
        # read mo coefficients and energies from restricted SCF calculation and duplicate
        Data = parseTurbomole(tm_dir + "/mos")
        orbe,orbs = Data["scfmo"]
        orbea, orbeb = orbe, orbe
        orbsa, orbsb = orbs, orbs
        """orbsa[:,i] holds mo coefficients of i-th spin-up orbital"""
        # read CIS coefficients
        Data = parseTurbomole(tm_dir + "/sing_a")
        Ecis,CISclosed = Data["eigenpairs"]
        #
#        print "got closed shell"
        # check format of CIS coefficients
        if (CISclosed.shape[0] == nocca*nvirta):
            CIS = np.vstack( (CISclosed, CISclosed) )/sqrt(2.0)
            """CIS coefficients are the same for spin-up and spin-down excitations, so we
            duplicate them and adjust the normalization:
                 CIS = vstack( (CISclosed, CISclosed) )/sqrt(2.0)"""
        elif (CISclosed.shape[0] == 2*nocca*nvirta):
            d = nocca*nvirta
            X = CISclosed[:d,:]
            Y = CISclosed[d:,:]
            XpY = (X+Y)/2.0
            CIS = np.vstack( (XpY, XpY) )/sqrt(2.0)
#            print "Closed shell with hybrid functional."
        else:
            # got something strange
            raise Exception("Format of CIS coefficients from closed shell calculation not understood. Wrong dimension %s" % str(CISclosed.shape))
    else:
        """uhf data"""
        # find occupation
        (nbfsp, nbfsc) = count_bfs(bs)
        nocca = Data["alpha"]#27
        noccb = Data["beta"] #26
        nvirta = nbfsp-nocca#135-nocca
        nvirtb = nbfsp-noccb#135-noccb
        # read mo coefficients and energies from unrestricted HF calculation
        Data = parseTurbomole(tm_dir + "/alpha")
        orbea,orbsa = Data["uhfmo_alpha"]
        Data = parseTurbomole(tm_dir + "/beta")
        orbeb,orbsb = Data["uhfmo_beta"]
        # read CIS coefficients
        Data = parseTurbomole(tm_dir + "/unrs_a")
        Ecis,CISopen = Data["eigenpairs"]
        #
#        print "got open shell"
        # check format of CIS coefficients
        if (CISopen.shape[0] == nocca*nvirta + noccb*nvirtb):
            """got CIS coefficients from calculation with non-hybrid functional, coefficients 
            for spin-up and spin-down are different"""
            CIS = CISopen
        elif (CISopen.shape[0] == 2*(nocca*nvirta + noccb*nvirtb)):
            """got CIS coefficients from calculation with hybrid functional such as pbe0"""
            d = nocca*nvirta + noccb*nvirtb
            X = CISopen[:d,:]
            Y = CISopen[d:,:]
            CIS = (X+Y)/2.0
#            print "Open shell with hybrid functional"
        else:
            # got something strange
            raise Exception("Format of CIS coefficients from open shell calculation not understood. Wrong dimension %s" % str(CISopen.shape))
    if functional in ["pbe", "slater-dirac-exchange", "s-vwn", "vwn", "s-vwn_gaussian", "pwlda", "becke-exchange", "b-lyp", "b-vwn", "lyp", "b-p", "pbe"]:
        """sqrt of eigenvalues in sing_a/unrs_a file are the excitation energies."""
        Ecis = sqrt(Ecis) 
    elif functional in ["bh-lyp", "b3-lyp", "pbe0"]:
        """eigenenergies equal eigenvalues"""
    else:
        raise Exception("Functional %s not supported." % functional)
    """The excitation energy of the ground state is 0 but needs to be inserted
    for consistency."""
    Ecis = np.insert(Ecis, 0, 0.0)
    """transform from Turbomole's internal mo format to canonical orbitals"""
    orbsa = canonical_mos(bs,orbsa)
    orbsb = canonical_mos(bs,orbsb)
    """check if mos make sense"""
#    check_mos(bs.bfs, orbsa)
#    check_mos(bs.bfs, orbsb)
    """CIS calculation is done in the basis spanned by single excitations,
       Insert the ground state."""
    nexs,nst = CIS.shape
    assert(nexs == nocca*nvirta + noccb*nvirtb)
    UCIS = np.zeros( (nexs+1, nst+1) )
    grst = np.zeros(nexs+1)
    grst[0] = 1.0
    UCIS[:,0] = grst
    """ground state = (1,0,0,0,...)"""
    for st in xrange(1,nst+1):
        vec = CIS[:,st-1]
        UCIS[0,st] = 0.0
        """projection on ground state is 0"""
        UCIS[1:,st] = vec[:]
    groundstate = [ [] ]
    excitations = groundstate + order_excs((nocca,noccb), (nvirta, nvirtb))
    Data = parseTurbomole(tm_dir + "/energy")
    en = Data["energy"]

    return( (mol, bs, en, (nocca, noccb), (nvirta, nvirtb), (orbea, orbeb), (orbsa, orbsb), Ecis, UCIS, excitations) )

def readTMorbitals(tm_dir):
    from PyQuante.Molecule import Molecule
    from PyQuante.Basis.basis import BasisSet
    """read TM molecular orbitals and energies"""
    # read coord file
    Data = parseTurbomole(tm_dir + "/coord")
    mol = Molecule('', atomlist = Data["coord"])
    # read basis
    Data = parseTurbomole(tm_dir + "/basis")
    bs = BasisSet(mol, basis_data = Data["basis"])

    Data = parseTurbomole(tm_dir + "/control")    
    if (not Data.has_key("uhf")):
        Data = parseTurbomole(tm_dir + "/mos")
        orbe,orbs = Data["scfmo"]
        orbea, orbeb = orbe, orbe
        orbsa, orbsb = orbs, orbs
    else:
        Data = parseTurbomole(tm_dir + "/alpha")
        orbea,orbsa = Data["uhfmo_alpha"]
        Data = parseTurbomole(tm_dir + "/beta")
        orbeb,orbsb = Data["uhfmo_beta"]
    """transform from Turbomole's internal mo format to canonical orbitals"""
    orbsa = canonical_mos(bs,orbsa)
    orbsb = canonical_mos(bs,orbsb)
    """check if mos make sense"""
#    check_mos(bs.bfs, orbsa)
#    check_mos(bs.bfs, orbsb)
    return (orbea, orbsa), (orbeb, orbsb)

############ Functions for writing coord and modifying control files##################
def writeTMcoords(mol):
    """translate the geometry of the Molecule mol in the Turbomole format"""
    text = "$coord\n"
    for at in mol.atoms:
        text += "\t%s %s %s %s\n" % \
            (     ("%0.12f" % at.r[0]).ljust(18), \
                  ("%0.12f" % at.r[1]).ljust(18), \
                  ("%0.12f" % at.r[2]).ljust(18), atom_names[at.atno-1].lower().ljust(3))
    text += "$user-defined bonds\n"
    text += "$end\n"
    return text

class Flag(_Block):
    def __init__(self, flagname, value):
        self.flagname = flagname
        self.value = value
    def tm_repr(self):
        return "$%s %s\n" % (self.flagname, self.value)

def set_block(lines, block_type, block):
    """lines should contain the lines from the control file.
    This function replaces the block with name block_type
    by a new block
    If block = None it removes the block.
    Imagine you wish to set the flag $exopt to the value 3. Then read the control
    files using readlines() and act with set_block() on the list of lines
    to include the block into the list:
       flblock = Flag("exopt", 3)
       lines = set_block(lines, "exopt", flblock)
    """
    modified_lines = []
    found = False
    skip = False
    for l in lines:
        l = l.rstrip()
        if l[0] == "$":
            skip = False
            parts = l[1:].split(" ", 1)
            if parts[0] == block_type:
                if block != None:
                    modified_lines += [l + "\n" for l in block.tm_repr().split("\n")[:-1]]
                found = True
                skip = True
                continue
            elif parts[0] == "end":
                continue
        if skip == False:
            modified_lines.append(l + "\n")
    if found == False and block != None:
        """flag was not present, so set it"""
        modified_lines += [l + "\n" for l in block.tm_repr().split("\n")[:-1]]
    modified_lines.append("$end\n")
    return modified_lines

########## reading other output files ########################
def read_spectrum(spectrum_file):
    """parse spectrum file"""
    fh = open(spectrum_file, "r")
    ex_energies = []
    osc_strengths = []
    for line in fh.readlines():
        if line[0] == "#":
            continue
        p = line.strip().split()
        en,osc = float(p[0]), float(p[1])
        ex_energies.append(en)
        osc_strengths.append(osc)
    fh.close()
    return ex_energies, osc_strengths

        
if __name__ == "__main__":
    import sys
    if (len(sys.argv) < 2):
        print "Usage: %s <TM directory>" % sys.argv[0]
        exit(-1)
    tm_dir = sys.argv[1]
    (mol,bs,en,nocc,nvirt,orbe,orbs,Ecis,UCIS,excitations) = importTurbomole(tm_dir)
    print mol
    print bs
    print "Occupation: occupied = %s   virtual = %s" % (nocc, nvirt)
    print "SCF ground state energy = %s" % en
    print "excitation energies = %s" % Ecis
