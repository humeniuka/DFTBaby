#!/usr/bin/env python
"""
Implementation of the DREIDING force field for QM/MM calculations with periodic boundary conditions
in molecular crystals. On top of the force field a Frenckel exciton model can be added. 
"""
import numpy as np
import numpy.linalg as la

from DFTB import XYZ, AtomicData

# class that wraps C-code
class PeriodicForceField:
    def __init__(self, atomlist, atomtypes, partial_charges, lattice_vectors,
                 chromophores, verbose=1, **kwds):
        """
        Builds a forces field with periodic boundary conditions.
        valence terms:
          bond stretching, angle bending, torsion and spectroscopic inversion
        non-bonding terms:
          van der Waals (Lennard Jones potential)
        The Coulomb interaction is not included, so all atoms should be approximately neutral. 
        The parameters are taken from the DREIDING force field. 

        Parameters:
        ===========
        atomlist: geometry which defined the bonding topology, 
          list of tuples (Zat,[x,y,z]) for each atom
        atomtypes: list with indeces of atom types
        partial_charges: list with partial charges 
        lattice_vectors: 3 lists with translation vectors of lattice
          [[ax,ay,az],[bx,by,bz],[cx,cy,cz]]
        chromophores: list of tuples (indeces, excitation_energies, transition_charges)
          for each chromophore (see `read_transition_charges()`)

        Optional
        ========
        connectivity_matrix: integer numpy array of shape (Nat,Nat) with connectivity
          matrix that is used to construct bonds, angles and dihedrals.
          If atoms I and J are connected, connectivity_matrix[I,J] = 1, otherwise 0.
          If absent, atoms are assumed to be connected if their distance is smaller than
          1.3 times the sum of the covalent radii.
        verbose: print additional information about the constructed force field
        """
        import ff
        self.force_field = ff.ForceField(atomlist, atomtypes, partial_charges,
                                         lattice_vectors, chromophores, verbose=verbose, **kwds)
        
    def getEnergyAndGradient(self, coords, state=0):
        """
        Computes the energy and gradient on the atoms in the central unit cell.
        The number and order of atoms should remain the same as used when defining the force field.

        Parameters:
        ===========
        coords: numpy array with shape (3*Nat), coords[3*i:3*(i+1)] are the x,y,z coordinates 
           of the i-th atom in bohr
        state: index of excitonic state, 0 stands for ground state

        Returns:
        ========
        energy: total force field energy in Hartree
        grads: numpy array with shape (3*Nat), grads[3*i:3*(i+1)] is the gradient on the i-th atom in Hartree/bohr
        """
        energy, grad = self.force_field.getEnergyAndGradient(coords, state)
        return energy, grad

    def getInternalCoordinateDefinitions(self):
        """ 
        definitions of all (redundant) internal coordinates

        Returns:
        ========
        coord_types  :  list of characters designating type of internal coordinate
                        ('B' - bond, 'A' - valence angle, 'D' - dihedral angle)
        atom_indices :  list of tuples with atom indices (starting at 0) involved in internal coordinate
                        (I,J) - bond , (I,J,K) - valence angle, (I,J,K,L) 
        """
        coord_types, _atom_indices = self.force_field.getInternalCoordinateDefinitions()
        atom_indices = []
        for k,ktype in enumerate(coord_types):
            # indices of atoms involved in internal coordinate k
            I, J, K, L = _atom_indices[4*k], _atom_indices[4*k+1], _atom_indices[4*k+2], _atom_indices[4*k+3]
            if ktype == 'B':
                bond = (I,J)
                atom_indices.append( bond )
            elif ktype == 'A':
                angle = (I,J,K)
                atom_indices.append( angle )
            elif ktype == 'D':
                dihedral = (I,J,K,L)
                atom_indices.append( dihedral )
            elif ktype == 'I':
                inversion = (I,J,K,L)
                atom_indices.append( inversion )
            else:
                raise ValueError("Unknown type of internal coordinate '%s'!" % ktype)
                
        return coord_types, atom_indices
                
    def getAngleIndices(self):
        """indices of internal coordinates which are angles (valence angles or dihedrals)"""
        coord_types, atom_indices = self.force_field.getInternalCoordinateDefinitions()
        angle_indices = np.where((coord_types == 'A') | (coord_types == 'D') | (coord_types == 'I'))[0]
        return angle_indices

    def getRedundantInternalCoordinates(self, coords, state=0):
        """
        Convert the cartesian coordinates into redundant internal coordinates and compute the 
        Wilson B-matrix.

        Parameters:
        ===========
        coords  : numpy array with shape (3*Nat), coords[3*i:3*(i+1)] are the x,y,z coordinates 
                  of the i-th atom in bohr
        
        Returns:
        ========
        internal  : numpy array with shape (M) with internal coordinates, first all bonds (in bohr), then
                    cosines of angles and torsions.
        bmatrix   : numpy array with shape (M,3*Nat) with Wilson B-matrix,
                    bmatrix[k,3*i:3*(i+1)] is the gradient of the k-th internal coordinate on
                    atom i.
        """
        internal, bmatrix = self.force_field.getRedundantInternalCoordinates(coords, 0)
        return internal, bmatrix
    
    def getExcitonStates(self):
        """
        Fetches the eigen energies and coefficients of the exciton wavefunctions.
        Calls to this functions must be preceeded by a call to .getEnergyAndGradient(...).
        
        Returns:
        ========
        en     : numpy array with shape (Nst), eigen energies of exciton states (in Hartree)
        coeffs : numpy array with shape (Nst,Nst), coefficients of the exciton wavefunction
                 coeffs[:,i] are the eigenvectors belonging to the eigenenergy en[i]
        """
        energies, coefficients = self.force_field.getExcitonStates()
        # check that wavefunctions are orthogonal
        olap = np.dot(coefficients.transpose(), coefficients)
        err = la.norm(olap - np.eye(len(energies)))
        assert err  < 1.0e-10, "exciton wavefunctions not orthogonal, |S - Id|= %e" % err
        return energies, coefficients
    def getTransitionDipoles(self, verbose=0):
        """
        Fetches the electric and magnetic transition dipole moments between the ground state (S0)
        and the excitonic eigenstates states.

        Returns:
        ========
        en: numpy array with shape (Nst), excitation energies in Hartree
        T : numpy array with shape (3,Nst), electric transition dipoles
        M : numpy array with shape (3,Nst), magnetic transition dipoles
        """
        elec_tdip, magn_tdip = self.force_field.getBasicTransitionDipoles()
        en, C = self.force_field.getExcitonStates()
        #  To get the transition dipoles between S0 and the excitonic
        # eigenstates, one has to multiply with the excitation coefficients
        T = np.dot(elec_tdip, C)
        M = np.dot(magn_tdip, C)

        if verbose > 0:
            print elec_tdip
            print magn_tdip
            print C
            print "transition dipoles"
            print T
            print "magnetic transition dipoles"
            print M

        return en, T, M
    
    # The following methods are only defined to make the interface compatible with Gaussian.UFF_handler
    def calc(self, atomlist, **kwds):
        coords = XYZ.atomlist2vector(atomlist)
        self.energy, self.grad = self.getEnergyAndGradient(coords)
    def get_MM_Energy(self):
        return self.energy
    def get_MM_Gradient(self):
        return self.grad

def save_exciton_spectrum(spec_file, en, T, M):
    """
    compute oscillator and rotational strengths and save absorption and
    circular dichroism spectra to file.

    Parameters
    ==========
    spec_file: filename where spectra are to be stored
    en: numpy array with shape (Nst), excitation energies in Hartree
    T : numpy array with shape (3,Nst), electric transition dipoles
    M : numpy array with shape (3,Nst), magnetic transition dipoles
    """
    # oscillator strengths and rotational strength
    nst = len(en)
    f = np.zeros(nst)
    R = np.zeros(nst)
    angleEM = np.zeros(nst)
    for i in range(0, nst):
        f[i] = 2.0/3.0 * en[i] * np.dot(T[:,i].conjugate(), T[:,i])
        # Im{ T_(0->i) * M[i->0] }
        R[i] = np.dot(T[:,i], M[:,i])
        # angle between T and M in radians
        TMnorm = la.norm(T[:,i])*la.norm(M[:,i])
        if TMnorm > 0.0:
            angleEM[i] = np.arccos( np.dot(T[:,i], M[:,i])/TMnorm )

    txt  = ""
    txt += "# Exciton             excitation energy         oscillator    rotatory\n"
    txt += "#  state           au      eV           nm             strength           E-M angle\n"
    txt += "#==================================================================================\n"
    for i in range(0, nst):
        txt += "   %4.1d          %7.5f  %7.5f  %9.2f    %7.5f    %+7.5f      %+7.3f\n" \
            % (i,
               en[i], en[i]*AtomicData.hartree_to_eV, AtomicData.hartree_to_nm/en[i],
               f[i], R[i], angleEM[i]*180.0/np.pi)
    txt += "\n"

    print txt
    fh = open(spec_file, "w")
    print>>fh, txt
    fh.close()
    print "exciton spectrum written to '%s'" % spec_file

def plot_exciton_spectrum(en, T, M, broadening=0.005, units="nm"):
    """
    plot absorption and circular dichroism spectra on top of each other

    Parameters
    ----------
    en: numpy array with shape (Nst), excitation energies in Hartree
    T : numpy array with shape (3,Nst), electric transition dipoles
    M : numpy array with shape (3,Nst), magnetic transition dipoles
    """
    # oscillator strengths and rotational strength
    nst = len(en)
    f = np.zeros(nst)
    R = np.zeros(nst)
    angleEM = np.zeros(nst)
    for i in range(0, nst):
        f[i] = 2.0/3.0 * en[i] * np.dot(T[:,i].conjugate(), T[:,i])
        # Im{ T_(0->i) * M[i->0] }
        R[i] = np.dot(T[:,i], M[:,i])


    import matplotlib.pyplot as plt
    from DFTB.Analyse.absorption_spectrum import broadened_spectrum, convert_energy

    ens, ab_spec = broadened_spectrum(en, f, broadening)
    ens, cd_spec  = broadened_spectrum(en, R, broadening)

    ax_cd = plt.subplot(211)
    ax_ab = plt.subplot(212, sharex = ax_cd)

    # convert energies from Hartree to nm
    ens = convert_energy(ens, "Hartree", units)
    
    fs = 18.0  # fontsize
    ax_cd.set_xlabel("Energy / %s" % units, fontsize=fs)
    # circular dichroism spectrum
    ax_cd.plot(ens, cd_spec)
    ax_cd.set_ylabel("Rotational Strength / arb. units", fontsize=fs)

    # absorption spectrum
    ax_ab.plot(ens, ab_spec)
    ax_cd.set_ylabel("Oscillator Strength / arb. units", fontsize=fs)

    plt.show()
    
def read_force_field(ff_file, units="Angstrom", fragment_id=AtomicData.atomic_number):
    """
    A force field file is structured like an xyz-file.
    1st line: number of atoms Nat
    2nd line: comment (line must not be empty)
    lines 3 to Nat+3: element name, x,y,z and atom type
       e.g.:      C   0.0 0.0 0.0    6
    The 6th column may optionally contain a partial charge
       e.g.:      O   0.0 0.0 0.0    15  -0.2
    The next lines are optional and contain the lattice vectors defining
    the unit cell starting with the symbol TV (translation vector):
       e.g.       Tv  1.0 0.0 0.0 
                  Tv  0.0 1.0 0.0
                  Tv  0.0 0.0 1.0


    Returns:
    ========
    atomlist: list of tuples (Zati,[xi,yi,zi])
    atomtypes: list of atom types for each atom
    partial_charges: list of partial charges for each atom
    lattice_vectors: list of at most 3 vectors defining the unit cell
    """
    fh = open(ff_file)
    line = fh.readline()
    # read number of atoms
    nat = int(line.split()[0])
    # skip comment
    fh.readline()
    #
    atomlist = []
    atomtypes = []
    partial_charges = []
    lattice_vectors = []
    for i in xrange(nat+3):
        # read nat atoms and at most 3 lattice vectors
        line = fh.readline()
        if not line:
            # end of file reached
            break
        words = line.split()
        x,y,z = map(float, words[1:4])
        if units == "Angstrom":
            x,y,z = map(lambda c: c/AtomicData.bohr_to_angs, [x,y,z])
        if words[0] == "Tv":
            lattice_vectors.append( [x,y,z] )
            continue
        atno = fragment_id(words[0])
        atomlist.append((atno, [x,y,z]))
        atomtypes.append(int(words[4]))
        if len(words) > 5:
            # 6th column contains partial charges
            partial_charges.append( float(words[5]) )
        else:
            partial_charges.append( 0.0 )    
    fh.close()
    if lattice_vectors == []:
        print "No lattice vectors found"
        # no lattice vectors were provided
        # HACK: By setting a lattice vectors to 0, we tell
        #       the function 'build_force_field' that we do not want
        #       a periodic calculation in this direction.
        #       If all lattice vectors are 0, only the atoms in the central
        #       unit cell are included (number of replica cells == 0)
        lattice_vectors = [ [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0] ]

    assert len(lattice_vectors) == 3, "Need 3 lattice vectors, got %d!" % len(lattice_vectors)
    return atomlist, atomtypes, partial_charges, lattice_vectors

def write_force_field(ff_file, atomlist, atomtypes, partial_charges, lattice_vectors=[],
                      units="Angstrom", title="", mode="w"):
    """
    save geometry with atom type assignments and partial charges

    Parameters:
    ===========
    ff_file: path to output file
    atomlist: list of tuples (Zati,[xi,yi,zi])
    atomtypes: list of atom types for each atom
    partial_charges: list of partial charges for each atom
    lattice_vectors: list of at most 3 vectors defining the unit cell, 
                     if the list is empty, [], no lattice vectors are saved
    """
    fh = open(ff_file, mode=mode)
    nat = len(atomlist)

    print>>fh, " %d" % nat
    print>>fh, " %s" % title
    for i in range(nat):
        Zi,posi = atomlist[i]
        atname = AtomicData.atom_names[Zi-1].capitalize()
        type_num = atomtypes[i]
        charge = partial_charges[i]
        x,y,z = posi[0], posi[1], posi[2]
        if units == "Angstrom":
            x,y,z = map(lambda c: c*AtomicData.bohr_to_angs, [x,y,z])
        print>>fh, " %4s    %15.10f %15.10f %15.10f    %3.1d  %+5.3f" % (atname.capitalize(), x,y,z, type_num, charge)
    # lattice vectors
    if len(lattice_vectors) == 3:
        for i in range(0, 3):
            print>>fh, "Tv    %10.7f   %10.7f   %10.7f" % tuple(lattice_vectors[i])
    fh.close()
    print "geometry, atomtypes and partial charges written to '%s'" % ff_file

def read_transition_charges(chromo_file, units="Angstrom"):
    """
    A .chromo file contains the excitation energies, atom-centered transition charges 
    and magnetic transition dipoles for different molecular fragments (chromophores) 
    that can be excited. This information is the basis for building an exciton model. 

    The file is structured like an xyz-file.
    1st line: number of atoms Nat
    2nd line: excitation energies (in eV), one for each of the Nst excited states
       e.g.    0.1 0.2 0.3
    line 3 to Nat+3: element name, x,y,z, ...
       The columns 1-4 are not used, they are only for convenience so that the file
       can be visualized in molden.
       The 5th column contains an index (starting from 0) into the atoms in the force field file
       The columns 6 to 6+4*Nst contain the transition charges q and magnetic dipoles [mx,my,mz]
       for the excited states I=1,...,Nst
       
              ELEMENT  x  y  z   INDEX   q(1) mx(1) my(1) mz(1) .... q(Nst) mx(Nst) my(Nst) mz(Nst)

    Returns
    =======
    iterator of tuples (indeces, excitation_energies, transition_charges) 
      for each geometry in the .chromo file
    """
    fh = open(chromo_file)
    while True:
        line = fh.readline()
        if not line:
            # end of file
            raise StopIteration
        # read number of atoms
        Nat = int(line.split()[0])
        # read excitation energies
        line = fh.readline()
        excitation_energies = map(lambda f: float(f) / AtomicData.hartree_to_eV, line.split())
        # number of excited states for this chromophore
        Nst = len(excitation_energies)
        indeces = []
        transition_charges = []
        for i in xrange(Nat):
            line = fh.readline()
            if not line:
                # end of file reached
                raise StopIteration
            words = line.split()
            # ignore columns 1-4
            # column 5 contains index into atoms in force field file
            index = int(words[4])
            indeces.append(index)
            # columns 6 to 6+4*Nst contain transition charges and magnetic transition dipoles
            values = map(float, words[5:])
            tqs = [] # list of [q,mx,my,mz] for each state on atom i
            for j in range(0, Nst):
                if len(values) >= 4*j+1:
                    q = values[4*j]  # transition charge on atom i for state j
                else:
                    q = 0.0
                if len(values) >= 4*j+4:
                    mx,my,mz = values[4*j+1: 4*j+4] # transition dipole on atom i for state j
                else:
                    mx,my,mz = 0.0, 0.0, 0.0
                tqs += [q,mx,my,mz]
            transition_charges.append(tqs)
            
        yield (indeces, excitation_energies, transition_charges)

def write_transition_charges(chromo_file, atomlist,
                             indeces, excitation_energies, transition_charges,
                             index_offset=0,
                             units="Angstrom", mode="w"):
    """
    write geometry, indeces, excitation energies and transition charges and 
    magnetic transition dipolesfor a single chromophore to a .chromo file

    Optional
    --------
    index_offset: this offset is added to the indeces
    mode: 'w' write or 'a' append
    """
    Nat = len(indeces)
    fh = open(chromo_file, mode)
    # number of atoms
    print>>fh, "%d" % Nat
    # excitation energies
    Nst = len(excitation_energies)
    print>>fh, "%s" % " ".join(map(lambda f: "%10.5f" % (f*AtomicData.hartree_to_eV),
                                   excitation_energies))
    #
    for i in range(0, Nat):
        atno,pos = atomlist[i]
        x,y,z = pos[0], pos[1], pos[2]
        if units == "Angstrom":
            x,y,z = map(lambda c: c*AtomicData.bohr_to_angs, [x,y,z])
        try:
            atname = AtomicData.atom_names[atno-1]
        except TypeError:
            # leave it as is
            atname = atno
        
        row = "%4s %15.10f %15.10f %15.10f    %5.1d    " \
              % (atname.capitalize(),x,y,z, indeces[i]+index_offset)
        for j in range(0, Nst):
            q = transition_charges[i][4*j]
            mvec = transition_charges[i][4*j+1:4*j+4]
            row += "%+7.5f  %+7.5f %+7.5f %+7.5f    " % (q, mvec[0], mvec[1], mvec[2])
        print>>fh, row
    
    fh.close()
        
def enlarged_unitcell(atomlist, lattice_vectors, nmax=1):
    """
    add neighbouring n cells to geometry
    """
    if nmax == 0:
        # only unitcell
        return atomlist, lattice_vectors
    a1 = np.array(lattice_vectors[0])
    a2 = np.array(lattice_vectors[1])
    a3 = np.array(lattice_vectors[2])
    # norms of unit vectors
    norm2 = np.zeros(3)
    norm2[0] = la.norm(a1)
    norm2[1] = la.norm(a2)
    norm2[2] = la.norm(a3)
    # If the norm of a lattice vector is 0, the periodicity in that direction is suppressed
    ns = [[0],[0],[0]]
    for dim in range(0, 3):
        if norm2[dim] == 0.0:
            ns[dim] = [0]
        else:
            # replicate unit cell along translation vector nmax times
            ns[dim] = range(-nmax,nmax+1)
    print "replica cells: %s" % map(len, ns)
            
    atomlist_enlarged = []
    for n1 in ns[0]:
        for n2 in ns[1]:
            for n3 in ns[2]:
                for Zi,posi in atomlist:
                    tvec = n1*a1 + n2*a2 + n3*a3
                    posi_translated = np.array(posi) + tvec
                    atomlist_enlarged.append( (Zi, posi_translated) )
    # lattice vectors of enlarged unit cell
    lattice_vectors_enlarged = np.zeros((3,3))
    for dim in range(0, 3):
        fac = len(ns[dim])
        lattice_vectors_enlarged[dim,:] = fac * np.array(lattice_vectors[dim])
        
    return atomlist_enlarged, lattice_vectors_enlarged.tolist()
                    
if __name__ == "__main__":
    import numpy as np
    from DFTB import XYZ, utils
    import sys
    from optparse import OptionParser

    usage = "Usage: python %s <.ff file>\n" % sys.argv[0]
    usage+= "  constructs force field and checks that numerical and analytical gradients agree\n"

    parser = OptionParser(usage)
    parser.add_option("--transition_charges", dest="transition_charges", type=str, default="", help="Path to .chromo file with excitation energies, transition charges and magnetic transition dipoles for building a Frenckel exciton model. [default: %default]")
    parser.add_option("--state", dest="state", type=int, default=0, help="Excitonic state that should be optimized. If no transition charges were provided, only the ground state (0) can be optimized. [default: %default]")
    
    (opts,args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)

    ff_file = args[0]  #"h2.ff" #"ethene.ff" #"pyrene_crystal_expanded.ff" #
    # read force field definition
    atomlist, atomtypes, partial_charges, lattice_vectors = read_force_field(ff_file)
    # read transition charge for exciton model (if available)
    if opts.transition_charges != "":
        chromophores = list(read_transition_charges(opts.transition_charges))
    else:
        chromophores = []
    pff = PeriodicForceField(atomlist, atomtypes, partial_charges, lattice_vectors, chromophores)
    pff.getInternalCoordinateDefinitions()
    coords = XYZ.atomlist2vector(atomlist)
    
    def energy_func(x):
        energy, grad = pff.getEnergyAndGradient(x, state=opts.state)
        return energy

    grad_num = utils.numerical_gradient(energy_func, coords, 1.0e-4)
    energy, grad_analytic = pff.getEnergyAndGradient(coords, state=opts.state)

    Nat = len(atomlist)
    print "numerical gradient"
    for i in range(0, Nat):
        print grad_num[3*i:3*(i+1)];
    print "analytical gradient"
    for i in range(0, Nat):
        print grad_analytic[3*i:3*(i+1)]

    grad_dif = grad_num - grad_analytic
    print "difference"
    for i in range(0, Nat):
        print grad_dif[3*i:3*(i+1)]

    print "error: %s" % np.sum(abs(grad_dif))

    # check the exciton states are orthogonal
    exc_energies, coefs = pff.getExcitonStates()
    print "excitation energies"
    print exc_energies
    print "coefficients of selected exciton wavefunction"
    print coefs[:,opts.state-1]
    print "Check that exciton wavefunctions are orthogonal"
    print "overlap"
    S = np.dot(coefs.transpose(), coefs)
    print S
    err = la.norm(S - np.eye(S.shape[0]))
    print "|S - Id|= %e" % err
    
    
