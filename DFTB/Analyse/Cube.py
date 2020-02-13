"""
 Code for writing Cube files
"""
from numpy import mgrid, sqrt, where
import numpy as np
from scipy import linalg as sla
import sys
import time

from DFTB.BasisSets import AtomicBasisSet
from DFTB import AtomicData
from DFTB.Timer import GlobalTimer as T
from DFTB.GammaApproximation import AuxiliaryBasisSet

def get_bbox(atomlist,**opts):
    """
    get_bbox(atoms,**opts) - find bounding box around molecule
    
    atoms  -  Molecule object

    Options     Default value     Description
    -----------------------------------------
    dbuff          5            extra space around molecule in a.u.

    return     -    limits of bbox ((xmin,xmax),(ymin,ymax),(zmin,zmax)) 
    """
    dbuff = opts.get('dbuff',5)
    big = opts.get('big',10000)
    xmin = ymin = zmin = big
    xmax = ymax = zmax = -big
    for (Zi, posi) in atomlist:
        x,y,z = posi
        xmin = min(xmin,x)
        ymin = min(ymin,y)
        zmin = min(zmin,z)
        xmax = max(xmax,x)
        ymax = max(ymax,y)
        zmax = max(zmax,z)
    xmin -= dbuff
    ymin -= dbuff
    zmin -= dbuff
    xmax += dbuff
    ymax += dbuff
    zmax += dbuff
    return (xmin,xmax),(ymin,ymax),(zmin,zmax)

def function_to_cubefile(atomlist, func, **opts):
    import sys
    (xmin,xmax),(ymin,ymax),(zmin,zmax) = get_bbox(atomlist, **opts)
    dx,dy,dz = xmax-xmin,ymax-ymin,zmax-zmin
    ppb = opts.get("ppb", 2.0) # Points per bohr
    filename = opts.get("filename", '')
    if (filename != ''):
        savefh = sys.stdout
        # print now writes directly to a file
        # all other print statements have to be changed to sys.stderr.write(...) !
        sys.stdout = open(filename, 'w')
    spacing = 1.0/ppb
    nx,ny,nz = int(dx*ppb),int(dy*ppb),int(dz*ppb)
    print "CUBE FILE"
    print "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    print "%5i %11.6f %11.6f %11.6f" %  (len(atomlist),xmin,ymin,zmin)
    print "%5i %11.6f %11.6f %11.6f" %  (nx,spacing,0,0)
    print "%5i %11.6f %11.6f %11.6f" %  (ny,0,spacing,0)
    print "%5i %11.6f %11.6f %11.6f" %  (nz,0,0,spacing)

    # The second record here is the nuclear charge, which differs from the
    #  atomic number when a ppot is used. Since I don't have that info, I'll
    #  just echo the atno
    for (Zi, posi) in atomlist:
        x,y,z = posi
        print "%5i %11.6f %11.6f %11.6f %11.6f" %  (Zi,Zi,x,y,z)
    grid = mgrid[xmin:xmax:spacing,ymin:ymax:spacing,zmin:zmax:spacing]

    ampGrid = func(grid, pow(spacing,3))

    for i in xrange(nx):
        for j in xrange(ny):
            for k in xrange(nz):
                amp = ampGrid[i][j][k]
                if abs(amp) < 1e-12: 
                    amp = 0.0
                print " %11.5e" % amp.real,
                if k % 6 == 5: print "\n ",
            print "\n ",
    """write grid to cube file"""
    # close file
    if (filename != ''):
        sys.stdout.close()
        sys.stdout = savefh



def orbital_amplitude(grid, bfs,mo, threshold=0.0, cache=True):
    """evaluate the amplitude of the orbital with coefficients mo in a
    basis set given by bfs on the grid (x,y,z)"""
    nbf = len(bfs)
    # It takes a long time to evaluate the basis functions on the grid.
    # Therefore bf_grid is stored as a static variable and whenever this function
    # is called the cached bf_grid is used unless the grid has changed
    
    # compare grids
    reuse_grid = False
    if hasattr(orbital_amplitude, "cached_grid"):
        diff = sum([np.sum(abs(grid[xyz] - orbital_amplitude.cached_grid[xyz])) for xyz in [0,1,2]])
        if diff < 1.0e-6:
            reuse_grid = True
    if reuse_grid == True:
        sys.stderr.write("Use cached grid\n")
        # use cached grid
        bf_grid = orbital_amplitude.bf_grid
    else:
        sys.stderr.write("Compute atomic orbital amplitudes on the grid\n")
        bf_grid = []
        for m in range(0, nbf):
            bf_grid.append( bfs[m].amp(*grid) )
        # static variables
        if cache == True:
            orbital_amplitude.cached_grid = grid
            orbital_amplitude.bf_grid = bf_grid
    sys.stderr.write("Compute molecular orbitals on the grid\n")
    xg,yg,zg = grid
    ampGrid = np.zeros(xg.shape, dtype=complex)

    MOmax = abs(mo).max()
    MOthresh = threshold * MOmax
    for ibf in xrange(nbf):
        if abs(mo[ibf]) > MOthresh:
            ampGrid += bf_grid[ibf] * mo[ibf]
    return(ampGrid)

def orbital2grid(atomlist,bfs,mo,**opts):
    """
    orbital2grid(atomlist,bfs,mo,**opts)   -   evaluate orbital amplitude on a grid and write to cube file
    
    atomlist
    bfs  - list of atomic basis functions
    mo  - vector of mo coefficients
    
    Option      Default Value         Description
    ppb              2.0           points per Bohr, resolution of grid
    filename      sys.stdout       if a file name is provided, a file handle is opened
                                   and the cube file is written there
    ampdens      "amplitude"       determines which quantity is printed 
                                   "amplitude" -> print amplitude of orbital
                                   "density" -> print absolute value squared
    threshold        0.0           neglect coefficients with relative weight below
                                   this threshold
    """
    ampdens = opts.get("ampdens", "amplitude")
    if (ampdens == "amplitude"):
        convert = lambda amp: amp
    else:
        convert = lambda amp: pow(abs(amp),2)
    def func(grid, dV):
        ampGrid = orbital_amplitude(grid, bfs, mo, threshold=opts.get("threshold",0.0))
        return convert(ampGrid)
    function_to_cubefile(atomlist, func, **opts)

def electron_density(grid, bfs, P, threshold=0.0):
    """
    Parameters:
    ===========
    grid: tuple (xg,yg,zg) with coordinates on the grid
    """
    rho = 0.0
    nbf = len(bfs)
    # It takes a long time to evaluate the basis functions on the grid.
    # Therefore bf_grid is stored as a static variable and whenever this function
    # is called the cached bf_grid is used unless the grid has changed
    
    # compare grids
    reuse_grid = False
    if hasattr(electron_density, "cached_grid") and grid[0].shape == electron_density.cached_grid[0].shape:
        diff = sum([np.sum(abs(grid[xyz] - electron_density.cached_grid[xyz])) for xyz in [0,1,2]])
        if diff < 1.0e-6:
            reuse_grid = True
    if reuse_grid == True:
#        sys.stderr.write("Use cached grid\n")
        # use cached grid
        bf_grid = electron_density.bf_grid
    else:
#        sys.stderr.write("Compute atomic orbital amplitudes on the grid\n")
        bf_grid = []
        for m in range(0, nbf):
            bf_grid.append( bfs[m].amp(*grid) )
        # static variables
        electron_density.cached_grid = grid
        electron_density.bf_grid = bf_grid
#    sys.stderr.write("Compute electron density on the grid\n")
    
#    for m in range(0, nbf):
#        Cm = 0.0
#        for n in range(0, nbf):
#            Cm += P[m,n] * bf_grid[n]      # NxN scalar*grid multiplications
#        rho += Cm * bf_grid[m].conjugate() # N grid*grid multiplications

    Pmax = abs(P).max()
    Pthresh = threshold * Pmax
    for m in range(0, nbf):
        for n in range(0, nbf):
            if abs(P[m,n]) > Pthresh:
                rho += P[m,n] * bf_grid[m].conjugate() * bf_grid[n] # NxN grid*grid multiplications
    return rho        


def density2grid(atomlist,bfs,P,**opts):
    """
    density2grid(atomlist,bfs,P,**opts)   -   evaluate electron density on a grid and write to cube file
    
    atomlist
    bfs  - list of atomic basis functions
    P - density matrix
    
    Option      Default Value         Description
    ppb              2.0           points per Bohr, resolution of grid
    filename      sys.stdout       if a file name is provided, a file handle is opened
                                   and the cube file is written there
    threshold                      atomic orbitals with coefficients below 
                                   threshold*(largest coefficient) are neglected when 
                                   evaluating MO amplitudes on a cube grid.

    """
    def func(grid, dV):
        ta = time.time()
        densGrid = electron_density(grid, bfs, P, \
                       threshold=opts.get("threshold", 0.0))
        tb = time.time()
        sys.stderr.write("Timing: %20.10f sec\n" % (tb-ta))
        return densGrid
    function_to_cubefile(atomlist, func, **opts)

def auxiliary_density2grid(atomlist, aux_bs, dq, ddip,**opts):
    """
    auxiliary_density2grid(aux_bs, dq, ddip,**opts)   -   evaluate electron density derived from the auxiliary basis on a grid
    
    atomlist
    aux_bfs  - auxiliary basis set
    dq       - partial Mulliken charges
    ddip     - partial Mulliken dipoles
    
    Option      Default Value         Description
    ppb              2.0           points per Bohr, resolution of grid
    filename      sys.stdout       if a file name is provided, a file handle is opened
                                   and the cube file is written there

    """
    def func(grid, dV):
        ta = time.time()
        x,y,z = grid
        densGrid = aux_bs.partial_density(dq, ddip, x,y,z)
        tb = time.time()
        sys.stderr.write("Timing: %20.10f sec\n" % (tb-ta))
        return densGrid
    function_to_cubefile(atomlist, func, **opts)


class CubeExporter:
    def __init__(self, dftb, name=""):            
        self.dftb = dftb
        self.name = name
        if self.dftb != None:
            self.bs = AtomicBasisSet(self.dftb.atomlist)
        self.ppb = 2.0
    def setGrid(self, points_per_bohr=2.0, cube_threshold=0.0):
        """
        points_per_bohr: separation between grid points
        cube_threshold: coefficients below threshold*(max coefficient) are neglected
           when evaluating density matrices or MOs on the grid
        """
        self.ppb = points_per_bohr
        self.cube_threshold=cube_threshold
    @T.timer
    def exportOrbital(self, orbital="HOMO", orbital_cubefile="/tmp/orbital.cube"):
        """
        orbital: name of the orbital whose amplitude should be written to a cube file, "HOMO" or "LUMO" or the integer index of the orbital
        orbital_cubefile: path to the output cube file
        """
        if orbital == None:
            return
        occ_index = where(self.dftb.f > 0.0)[0]
        virt_index = where(self.dftb.f == 0.0)[0]
        if orbital == "HOMO":
            orb = occ_index[-1]
        elif orbital == "LUMO":
            orb = virt_index[0]
        elif orbital == "HOMO-1":
            orb = occ_index[-2]
        elif orbital == "LUMO+1":
            orb = virt_index[+1]
        else:
            orb = int(orbital)
        orbital2grid(self.dftb.atomlist, self.bs.bfs, self.dftb.orbs[:,orb], \
                     filename=orbital_cubefile, ppb=self.ppb,\
                     threshold=self.cube_threshold)
    @T.timer
    def exportDensity(self, density_cubefile="/tmp/density.cube"):
        """
        density_cubefile: path to cube file
        points_per_bohr: separation between grid points
        """
        if density_cubefile == None:
            return
        density2grid(self.dftb.atomlist, self.bs.bfs, self.dftb.P, \
                         filename=density_cubefile, ppb=self.ppb, \
                         threshold=self.cube_threshold)
    def exportPartialDensity(self, density_cubefile="/tmp/partial_density.cube"):
        """
        density_cubefile: path to cube file
        points_per_bohr: separation between grid points
        """
        if density_cubefile == None:
            return
        density2grid(self.dftb.atomlist, self.bs.bfs, self.dftb.P-self.dftb.P0, \
                         filename=density_cubefile, ppb=self.ppb, \
                         threshold=self.cube_threshold)
    def exportAuxiliaryDensity(self, density_cubefile="/tmp/auxiliary_density.cube"):
        if density_cubefile == None:
            return
        bs_aux = AuxiliaryBasisSet(self.dftb.atomlist, self.dftb.hubbard_U)
        auxiliary_density2grid(self.dftb.atomlist, bs_aux, self.dftb.dq, self.dftb.ddip, \
                         filename=density_cubefile, ppb=self.ppb, \
                         threshold=self.cube_threshold)
    def exportCubes(self, cubedir=None, cube_orbitals="[]", points_per_bohr=2.0, cube_threshold=1.0e-3):
        """
        Cube-Files.cubedir: directory where cube files are stored
        Cube-Files.cube_orbitals: list indeces of molecular orbitals which should be saved to a cube file, e.g. '[1,2,3]' (if None, the HOMO-1, HOMO, LUMO and LUMO+1 are saved)
        Cube-Files.points_per_bohr: resolution of grid for cube files
        Cube-Files.cube_threshold: atomic orbitals with coefficients below threshold*(largest coefficient) are neglected when evaluating MO amplitudes on a cube grid.
        """
        if cubedir == None:
            return
        from os.path import join, expandvars, expanduser
        cubedir = expanduser(expandvars(cubedir))
        self.setGrid(points_per_bohr=points_per_bohr, cube_threshold=cube_threshold)
        if cube_orbitals == []:
            print "Write cube files for HOMO-1, HOMO, LUMO, LUMO+1"
            self.exportOrbital("HOMO-1", orbital_cubefile=join(cubedir, "%s_HOMO-1.cube" % self.name))
            self.exportOrbital("HOMO", orbital_cubefile=join(cubedir, "%s_HOMO.cube" % self.name))
            self.exportOrbital("LUMO", orbital_cubefile=join(cubedir, "%s_LUMO.cube" % self.name))
            self.exportOrbital("LUMO+1", orbital_cubefile=join(cubedir, "%s_LUMO+1.cube" % self.name))
        else:
            print "Write cube files for orbitals %s" % cube_orbitals
            for o in cube_orbitals:
                o = int(o)
                self.exportOrbital(o-1, orbital_cubefile=join(cubedir, "%s_orbital_%.4d.cube" % (self.name, o)))
        print "Write cube file for density"
        self.exportDensity(density_cubefile=join(cubedir, "%s_density.cube" % self.name))
        print "Write cube file for partial density"
        self.exportPartialDensity(density_cubefile=join(cubedir, "%s_partial_density.cube" % self.name))
        print "Write cube file for auxiliary density"
        self.exportAuxiliaryDensity(density_cubefile=join(cubedir, "%s_auxiliary_density.cube" % self.name))
        

class CubeExporterEx(CubeExporter):
    """
    for excited state calculations
    """
    def __init__(self, tddftb, name=""):
        self.tddftb = tddftb
        self.dftb = tddftb.dftb2
        self.name = name
        if self.dftb != None:
            self.bs = AtomicBasisSet(self.dftb.atomlist)
        self.ppb = 2.0
        # auxiliary basis of Gaussian flucuation functions F_A(r)
        self.bs_aux = AuxiliaryBasisSet(self.dftb.atomlist, self.dftb.hubbard_U)
    @T.timer
    def exportTransitionDensity(self, state, \
                                tdense_cubefile="/tmp/transition_density.cube", \
                                save_transition_density=0):
        """
        Write the transition density for the excited state to a file.

        state: index of excited state (starting from 0)
        tdense_cubefile: path to cube file
        """
        Ptrans = self.tddftb.TransitionDensityMatrix(state)
        if save_transition_density == 1:
            # transition density P_ab in AO basis
            tdense_file = tdense_cubefile[:-5] + "-AO" + ".mat"
            np.savetxt(tdense_file, Ptrans)
            print "Transition density matrix in AO basis for state %d written to %s" % (state+1, tdense_file)
            # orthogonalize transition density
            # P  ->  S^(1/2). P . S^(1/2)^T
            # If the orthogonalized transition densities for states I and J are interpreted as vectors in the
            # space R^[nocc x nvirt], they form an orthonormal basis:
            #  sum_(a,b) P^I_(a,b) * P^J_(a,b) = delta_IJ
            #
            Ptrans_ortho = np.dot(self.Ssq, np.dot(Ptrans, self.Ssq.transpose()))
            tdense_file_ortho = tdense_cubefile[:-5] + ".mat"
            np.savetxt(tdense_file_ortho, Ptrans_ortho)
            print "Orthogonalized transition density matrix in AO basis for state %d written to %s" % (state+1, tdense_file_ortho)
        density2grid(self.dftb.atomlist, self.bs.bfs, Ptrans, \
                         filename=tdense_cubefile, ppb=self.ppb, \
                         threshold=self.cube_threshold)
    
    def exportDifferenceDensity(self, state, \
                                difdense_cubefile="/tmp/difference_density.cube"):
        """
        Write the density difference between the excited state and the ground state to a file.

        Parameters
        ----------
        state             : index of excited state (starting from 0)
        difdense_cubefile : path to cube file
        """
        P0, PI = self.tddftb.ExcitedDensityMatrix(state)
        # difference density matrix
        Pdif = PI-P0
        
        density2grid(self.dftb.atomlist, self.bs.bfs, Pdif, \
                     filename=difdense_cubefile, ppb=self.ppb, \
                     threshold=self.cube_threshold)
        
    def exportDifferenceDensity_approx(self, state, \
                                difdense_cubefile="/tmp/difference_density.cube"):
        """
        Write the density difference between the excited state and the ground state to a file.
        The approximations of eqn. (F9) and F(10) in reference [1] are used which depend on 
        particle and hole charges.

        Parameters
        ----------
        state             : index of excited state (starting from 0)
        difdense_cubefile : path to cube file
        
        References
        ----------
        [1] A. Humeniuk and R. Mitric,
            "Long-range correction for tight-binding TD-DFT"
            J. Chem. Phys.143, 134120 (2015)
        """
        # particle-hole charges
        particle_charges, hole_charges = self.tddftb.ParticleHoleCharges(state)
        # partial charges on ground state
        dq0 = self.tddftb.dftb2.getPartialCharges()
        # partial charges on excited state
        dqI = particle_charges + hole_charges + dq0
        auxiliary_density2grid(self.dftb.atomlist, self.bs_aux, dqI, 0.0*self.dftb.ddip, \
                               filename=difdense_cubefile, ppb=self.ppb, \
                               threshold=self.cube_threshold)
    def exportCubes(self, cubedir=None, cube_states="[1,2,3]", density_type="transition", points_per_bohr=2.0, cube_threshold=1.0e-3, save_transition_density=0):
        """
        Cube-Files.cubedir: directory where cube files of transition or difference densities are stored
        Cube-Files.cube_states: indeces of excited states for which transition or difference densities should be computed, the first excited state has index 1!
        Cube-Files.density_type: specifies whether transition densities ('transition') or density differences ('difference') between the excited states and the ground state are calculated
        Cube-Files.points_per_bohr: resolution of grid for cube files
        Cube-Files.cube_threshold: atomic orbitals with coefficients below threshold*(largest coefficient) are neglected when evaluating MO amplitudes on a cube grid
        Cube-Files.save_transition_density: if 1, the orthogonalized transition density P = S^(1/2).P.S^(1/2)^T in the AO basis is also saved for each state in cube_states as a matrix in a plain text file with suffix '.mat'
        """
        if cubedir == None:
            return
        print "Export cube files for excited states %s" % cube_states
        from os.path import join, expandvars, expanduser
        cubedir = expanduser(expandvars(cubedir))
        self.setGrid(points_per_bohr=points_per_bohr, cube_threshold=cube_threshold)
        #
        if save_transition_density == 1:
            # compute the matrix square root of the overlap matrix S^(1/2)
            S = self.tddftb.dftb2.getOverlapMatrix()
            self.Ssq = sla.sqrtm(S)
        for st in cube_states:
            if density_type == "transition":
                self.exportTransitionDensity(st-1, \
                                         tdense_cubefile=join(cubedir, "%s_tdense_%.4d%s.cube" % (self.name, st, self.tddftb.Irreps[st-1])), \
                                         save_transition_density=save_transition_density)
            elif density_type == "difference":
                self.exportDifferenceDensity(st-1, \
                                         difdense_cubefile=join(cubedir, "%s_difdense_%.4d%s.cube" % (self.name, st, self.tddftb.Irreps[st-1])))
            else:
                raise ValueError("option 'density_type' can be 'difference' or 'transition' but not %s" % density_type)
###### READING AND WRITING GRID DATA TO CUBE FILES ###########
def readCube(cubefile):
    """
    Parameters:
    ===========
    cubefile: path to cube file

    Returns:
    ========
    atomlist: list of tuples (Zi, posi) with atomic positions
    origin: origin of axes, 3d vector
    axes: 3 vectors spanning a voxel
    data: values for each voxel,   data[i,j,k] is the value of the voxel
       at position  x_ijk = origin + i*axes[0] + j*axes[1] + k*axes[2]
    """
    fh = open(cubefile)
    l = fh.readlines()
    fh.close()
    # first two lines
    #assert l[0].strip() == "CUBE FILE"
    #assert l[1].strip() == "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    # number of atoms and origin
    parts = l[2].split()
    natoms, origin = int(parts[0]), np.array(map(float, parts[1:4]))
    # voxels and axis vectors
    Nvox = [0, 0, 0]
    axes = [None, None, None]
    for xyz in range(0, 3):
        parts = l[3+xyz].split()
        Nvox[xyz], axes[xyz] = int(parts[0]), np.array(map(float, parts[1:]))

    atomlist = []
    for i in range(0, natoms):
        parts = l[6+i].split()
        Zi, posi = int(parts[0]), map(float, parts[2:])
        atomlist.append((Zi, posi))
    # read grid into long list of floats
    data1d = []
    for dataline in l[7+i:]:
        data1d += map(float, dataline.split())
    Nx,Ny,Nz = Nvox
    assert len(data1d) == Nx*Ny*Nz
    data = np.reshape(data1d, (Nx,Ny,Nz))
    return atomlist, origin, axes, data

def get_points_and_values(origin, axes, data):
    """
    compute positions of grid points in a cube file and associated values

    Parameters
    ----------
    origin: 3d vector with origin of cube file
    axes: 3 vectors spanning a voxel
    data: values for each voxel,   data[i,j,k] is the value of the voxel
       at position  x_ijk = origin + i*axes[0] + j*axes[1] + k*axes[2]

    Returns
    -------
    points: (x,y,z) = points[i] are the coordinates of the i-th grid point
    values: values[i] is the value of the cube at the i-th grid point
    """
    Nx,Ny,Nz = data.shape
    c = 0
    points = np.zeros((3,Nx*Ny*Nz))
    values = np.zeros(Nx*Ny*Nz)
    for i in xrange(Nx):
        for j in xrange(Ny):
            for k in xrange(Nz):
                points[:,c] = origin + i*axes[0] + j*axes[1] + k*axes[2]
                values[c] = data[i,j,k]
                c += 1
    return points, values

def writeCube(cubefile, atomlist, origin, axes, data):
    """
    write grid data to cube file

    Parameters:
    ===========
    cubefile: path to the cube file
    atomlist: list of tuples (Zi,posi) with nuclear geometry
    origin: 3d vector with origin of cube file
    axes: 3 vectors spanning a voxel
    data: values for each voxel,   data[i,j,k] is the value of the voxel
       at position  x_ijk = origin + i*axes[0] + j*axes[1] + k*axes[2]
    """
    fh = open(cubefile, "w")
    nx,ny,nz = data.shape
    xmin,ymin,zmin = origin
    print>>fh, "CUBE FILE"
    print>>fh, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    print>>fh, "%5i %11.6f %11.6f %11.6f" %  (len(atomlist),xmin,ymin,zmin)
    print>>fh, "%5i %11.6f %11.6f %11.6f" %  (nx,axes[0][0],axes[0][1],axes[0][2])
    print>>fh, "%5i %11.6f %11.6f %11.6f" %  (ny,axes[1][0],axes[1][1],axes[1][2])
    print>>fh, "%5i %11.6f %11.6f %11.6f" %  (nz,axes[2][0],axes[2][1],axes[2][2])

    for (Zi, posi) in atomlist:
        x,y,z = posi
        print>>fh, "%5i %11.6f %11.6f %11.6f %11.6f" %  (Zi,Zi,x,y,z)

    for i in xrange(nx):
        for j in xrange(ny):
            for k in xrange(nz):
                print>>fh, " %11.5e" % data[i,j,k],
                if k % 6 == 5: print>>fh, "\n ",
            print>>fh, "\n ",

    fh.close()

if __name__ == "__main__":
    from numpy import array, mgrid, ones, diag, zeros

    atomlist = [(6, array([0.1, 1.0, 2.0]))]
    bs = AtomicBasisSet(atomlist)
    dx = 0.40
    grid = mgrid[-6.0:6.0:dx,-6.0:6.0:dx,-6.0:6.0:dx]
    nbfs = len(bs.bfs)
    for i,bfs in enumerate(bs.bfs):
        mo = zeros(nbfs)
        mo[i] = 1.0
        orbital2grid(atomlist, bs.bfs, mo, filename="/tmp/orbital_%d.cube" % i)
