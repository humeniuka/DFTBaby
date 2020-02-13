#!/usr/bin/env python
"""
a PyQt4 widget that can display molecular structures and volumetric data
"""

# First, and before importing any Enthought packages, set the ETS_TOOLKIT
# environment variable to qt4, to tell Traits that we will use Qt.
import os
os.environ['ETS_TOOLKIT'] = 'qt4'

# To be able to use PySide or PyQt4 and not run in conflicts with traits,
# we need to import QtGui and QtCore from pyface.qt
from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
        SceneEditor

import numpy as np
import numpy.linalg as la
from DFTB import XYZ, AtomicData
from DFTB.Analyse import MolecularGraph
from DFTB.Modeling import MolecularCoords as MolCo
from DFTB.Modeling.MinimalEnclosingBox import MoleculeBox
from DFTB.Modeling import BondOrders
from DFTB import Parameters
from DFTB import ImplicitSolvent
from DFTB import NACsApprox
from DFTB.optparse import OptionParserFuncWrapper

class CubeData:
    """
    This class holds volumetric data 
    """
    def __init__(self):
        self.data = []
        self.atomlist = []
        self.grid = []
    def loadFromFile(self, cube_file):
        fh = open(cube_file)
        # first two lines contain comments
        file_type = fh.readline().strip()
        #assert file_type == "CUBE FILE", "%s does not seem to be a cube file!" % cube_file
        comment = fh.readline()
        # read number of atoms and lower corner of cube
        parts = fh.readline().split()
        nat, xmin, ymin, zmin = abs(int(parts[0])), float(parts[1]), float(parts[2]), float(parts[3])
        # axes of a voxel
        parts = fh.readline().split()
        n1, dx1, dy1, dz1 = int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
        assert abs(dy1) == 0.0 and abs(dz1) == 0.0, "x-axis of voxels not parallel to global x-axis"
        parts = fh.readline().split()
        n2, dx2, dy2, dz2 = int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
        assert abs(dx2) == 0.0 and abs(dz2) == 0.0, "y-axis of voxels not parallel to global y-axis"
        parts = fh.readline().split()
        n3, dx3, dy3, dz3 = int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
        assert abs(dx3) == 0.0 and abs(dy3) == 0.0, "z-axis of voxels not parallel to global z-axis"
        # volume element
        self.dV = dx1 * dy2 * dz3
        # read atomic positions
        self.atomlist = []
        for i in range(0, nat):
            parts = fh.readline().split()
            Zi = int(parts[0])
            pos = map(float, parts[-3:])
            self.atomlist.append( (Zi, pos) )
        # read data
        lines = fh.readlines()
        try:
            # In some cube files there is an extra line after the atomic position
            # that contains two integers
            #  e.g.   1 2201
            int1, int2 = map(int, lines[0].split())
            # skipt extra line
            lines = lines[1:]
        except ValueError:
            pass
        data_str = ' '.join(lines)
        fh.close()

        self.data = np.fromstring(data_str, sep=' ')
        assert len(self.data) == n1*n2*n3, "Not enough data points in cube file, got %d but expected %d." % (len(self.data), n1*n2*n3)
        self.data.shape = (n1,n2,n3)

        xmax, ymax, zmax = xmin+dx1*n1, ymin+dy2*n2, zmin+dz3*n3
        self.grid = np.mgrid[xmin:xmax:n1*1j, ymin:ymax:n2*1j, zmin:zmax:n3*1j]
        # 
    def integrateIso(self, iso):
        """
        sum over all values in the cube file which have values >= iso or values <= -iso
        and weigh them with the volume element
        """
        # integral over f >= iso
        Iplus = np.sum(self.data[self.data >= iso]) * self.dV
        # integral over f <= iso
        Iminus = np.sum(self.data[self.data <= -iso]) * self.dV
        return Iplus, Iminus
        

atom_colours_rgb = {"h": (0.7, 0.7, 0.7),
                    "c": (0.5019607843137255, 0.5019607843137255, 0.5019607843137255),
                    "o": (1.0, 0.0, 0.0),
                    "n": (0.0, 0.0, 1.0),
                    "s": (1.0, 1.0, 0.0),
                    "p": (0.0, 1.0, 0.0),
                    }
################################################################################
import time

def disable_rendering(func):
    def wrapper(self, *args, **kwds):
        mlab = self.scene.mlab
        fig = mlab.gcf()
        fig.scene.disable_render = True
        ta = time.time()
        func(self, *args, **kwds)
        tb = time.time()
        #print "Function %s took %s seconds" % (func.__name__, tb-ta)
        fig.scene.disable_render = False
    return wrapper

class MoleculeObject:
    def __init__(self, scene, show_flags):
        self.scene = scene
        self.show_flags = show_flags
        # keep track of what is actually shown
        self.shown_atoms = False
        self.shown_lewis_structure = False
        
        self.primitives = []
    def draw(self, atomlist, animation=False):
        need_recreate = False
        if animation == False or len(self.primitives) == 0:
            need_recreate = True
        if self.show_flags["atoms"] != self.shown_atoms:
            need_recreate = True
        if self.show_flags["Lewis structure"] != self.shown_lewis_structure:
            need_recreate = True
        ##
        if need_recreate == False:
            self._update_visualization(atomlist)
        else:
            self.remove()
            self._create_visualization(atomlist)
    def remove(self):
        for o in self.primitives:
            o.remove()
        self.primitives = []
    @disable_rendering
    def _create_visualization(self, atomlist):
        mlab = self.scene.mlab
        # draw atoms as balls
        vec = XYZ.atomlist2vector(atomlist)
        x, y, z = vec[::3], vec[1::3], vec[2::3]
        Zs = np.array([Z for (Z,pos) in atomlist])
        atom_names = [AtomicData.atom_names[Z-1] for (Z,pos) in atomlist]
        rads = np.array([AtomicData.covalent_radii.get(atname) for atname in atom_names])
        atom_radii = np.sqrt(rads)*1.8

        s = atom_radii
        if self.show_flags["atoms"] == False:
            # make atoms so small that one does not see them
            scale_factor = 0.01
            self.shown_atoms = False
        else:
            scale_factor = 0.45
            self.shown_atoms = True
        atoms = mlab.quiver3d(x,y,z, s,s,s, scalars=Zs.astype(float),
                              mode="sphere", scale_factor=scale_factor, resolution=20,
                              figure=self.scene.mayavi_scene)

        # atoms are coloured by their atomic number
        atoms.glyph.color_mode = "color_by_scalar"
        atoms.glyph.glyph_source.glyph_source.center = [0,0,0]
        self.lut = atoms.module_manager.scalar_lut_manager.lut.table.to_array()
        for atname in set(atom_names):
            Z = AtomicData.atomic_number(atname)
            self.lut[Z,0:3] = ( np.array( atom_colours_rgb.get(atname, (0.0, 0.75, 0.75)) )*255.0 ).astype('uint8')
        atoms.module_manager.scalar_lut_manager.lut.table = self.lut
        atoms.module_manager.scalar_lut_manager.data_range = (0.0, 255.0)
        
        # draw bonds
        C = XYZ.connectivity_matrix(atomlist)
        Nat = len(atomlist)
        connections = []
        bonds = mlab.points3d(x,y,z, scale_factor=0.15, resolution=20,
                              figure=self.scene.mayavi_scene)
        for i in range(0, Nat):
            Zi, posi = atomlist[i]
            for j in range(i+1,Nat):
                if C[i,j] == 1:
                    Zj, posj = atomlist[j]
                    connections.append( (i,j) )
        bonds.mlab_source.dataset.lines = np.array(connections)
        bonds.mlab_source.dataset.modified()
        
        tube = mlab.pipeline.tube(bonds, tube_radius=0.15, #tube_radius=0.05,
                                  figure=self.scene.mayavi_scene)
        tube.filter.radius_factor = 1.0
        surface = mlab.pipeline.surface(tube, color=(0.8,0.8,0.0),
                                        opacity=0.7,
                                        figure=self.scene.mayavi_scene)

        self.atoms = atoms
        self.bonds = bonds
#        self.tube = tube
#        self.surface = surface
        self.atom_radii = s
        self.Zs = Zs.astype(float)
        self.primitives.append(atoms)
        self.primitives.append(bonds)
        self.primitives.append(tube)
        self.primitives.append(surface)
        # Lewis structure
        if self.show_flags["Lewis structure"] == True:
            bondsTuples, bond_orders, lone_pairs, formal_charges = BondOrders.assign_bond_orders(atomlist, C, charge=0)
            # plot DOUBLE bonds
            double_bonds = mlab.points3d(x,y,z, scale_factor=0.15, resolution=20, figure=self.scene.mayavi_scene)
            connections_double = []
            for i,bond in enumerate(bondsTuples):
                if bond_orders[i] == 2:
                    # double bond between atoms I and J
                    connections_double.append( bond )
            double_bonds.mlab_source.dataset.lines = np.array(connections_double)
            double_bonds.mlab_source.dataset.modified()
            
            tube = mlab.pipeline.tube(double_bonds, tube_radius=0.35, figure=self.scene.mayavi_scene)
            tube.filter.radius_factor = 1.0
            surface = mlab.pipeline.surface(tube, color=(0.8,0.8,0.0), figure=self.scene.mayavi_scene)
            #
            self.double_bonds = double_bonds
            self.primitives.append(double_bonds)
            self.primitives.append(tube)
            self.primitives.append(surface)
            #
            # plot TRIPLE bonds
            triple_bonds = mlab.points3d(x,y,z, scale_factor=0.15, resolution=20, figure=self.scene.mayavi_scene)
            connections_triple = []
            for i,bond in enumerate(bondsTuples):
                if bond_orders[i] == 3:
                    # triple bond between atoms I and J
                    connections_triple.append( bond )
            triple_bonds.mlab_source.dataset.lines = np.array(connections_triple)
            triple_bonds.mlab_source.dataset.modified()
            
            tube = mlab.pipeline.tube(triple_bonds, tube_radius=0.35, figure=self.scene.mayavi_scene)
            tube.filter.radius_factor = 1.0
            surface = mlab.pipeline.surface(tube, color=(0.8,0.8,0.0), figure=self.scene.mayavi_scene)
            #
            self.triple_bonds = triple_bonds
            self.primitives.append(triple_bonds)
            self.primitives.append(tube)
            self.primitives.append(surface)

            self.shown_lewis_structure = True
    def _update_visualization(self, atomlist):
        """
        only change the underlying data but do not recreate visualization.
        This is faster and avoids jerky animations.
        """
        vec = XYZ.atomlist2vector(atomlist)
        x, y, z = vec[::3], vec[1::3], vec[2::3]
        s = self.atom_radii
        if self.show_flags["atoms"] == False:
            # make atoms so small that one does not see them
            scale_factor = 0.01
            self.shown_atoms = False
        else:
            scale_factor = 0.3
            self.shown_atoms = True
        # update atom positions
        self.atoms.mlab_source.set(x=x,y=y,z=z,u=s,v=s,w=s, scalars=self.Zs, scale_factor=scale_factor)
        # atoms are coloured by their atomic number
        self.atoms.glyph.color_mode = "color_by_scalar"
        self.atoms.glyph.glyph_source.glyph_source.center = [0,0,0]
        self.atoms.module_manager.scalar_lut_manager.lut.table = self.lut
        self.atoms.module_manager.scalar_lut_manager.data_range = (0.0, 255.0)

        # update bonds
        self.bonds.mlab_source.reset(x=x,y=y,z=z, scale_factor=1.0)

        C = XYZ.connectivity_matrix(atomlist)
        Nat = len(atomlist)
        connections = []
        for i in range(0, Nat):
            Zi, posi = atomlist[i]
            for j in range(i+1,Nat):
                if C[i,j] == 1:
                    Zj, posj = atomlist[j]
                    connections.append( (i,j) )
        self.bonds.mlab_source.dataset.lines = np.array(connections)
        self.bonds.mlab_source.dataset.modified()

class ChargeCloudRadii:
    def __init__(self):
        hubbard_U_byZ = Parameters.get_hubbard_parameters("homegrown")
        self.sigmas = {}  # width of Gaussian charge fluctuations F(r)
        self.radii = {}  # radius where, F(r0) = F(0) * isovalue
        for Z,U in hubbard_U_byZ.iteritems():
            self.sigmas[Z] = 1.0/(np.sqrt(np.pi)*U)
    def setIsoValue(self, iso):
        logIso = np.log(iso)
        for Z,sA in self.sigmas.iteritems():
            r0 = np.sqrt(-2*sA**2 * logIso)
            self.radii[Z] = r0
    def getRadii(self, atomlist):
        r0s = []
        for Z,pos in atomlist:
            r0s.append( self.radii[Z] )
        return r0s
        
class ChargesObject:
    def __init__(self, scene, show_flags, charge_cloud_radii):
        self.scene = scene
        self.show_flags = show_flags
        self.charge_cloud_radii = charge_cloud_radii
        # keep track of what is actually shown
        self.shown_charges = False
        self.shown_frag_charges = False
        self.shown_charge_clouds = False
        self.shown_dipole_moment = False
        self.shown_enclosing_boxes = False
        self.shown_solvent_screening_charges = False
        self.shown_nacs = False
        
        self.nr_fragments = 0
        self.charge_labels = []  # list of objects created by mlab
        self.frag_charge_labels = []
        self.dipole_vectors = []
        self.charge_clouds = []
        self.frag_enclosing_boxes = []
        self.sas_points = []
        self.nac_vectors = []
    def draw(self, atomlist, charges, fragments=[], animation=False):
        #
        if type(charges) != type(None) and len(charges.shape) == 2:
            # sum particle and hole charges 
            charges = np.sum(charges, axis=1)
        #
        
        need_recreate = False
        if animation == False or len(self.charge_labels + self.frag_charge_labels + self.charge_clouds) == 0:
            # draw for the first time
            need_recreate = True
        if self.show_flags["frag. charges"] == True and len(fragments) != self.nr_fragments:
            # number of fragments has changed
            need_recreate = True
        # view differs from what should be shown
        if self.show_flags["charges"] != self.shown_charges:
            need_recreate = True
        if self.show_flags["frag. charges"] != self.shown_frag_charges:
            need_recreate = True
        if self.show_flags["charge clouds"] != self.shown_charge_clouds:
            need_recreate = True
        if self.show_flags["dipole moment"] != self.shown_dipole_moment:
            need_recreate = True
        if self.show_flags["enclosing box"] != self.shown_enclosing_boxes:
            need_recreate = True
        if self.show_flags["screening charges (COSMO)"] != self.shown_solvent_screening_charges:
            need_recreate = True
        if self.show_flags["non-adiab. coupling vectors"] != self.shown_nacs:
            need_recreate = True
        ##
        if need_recreate == False:
            self._update_visualization(atomlist, charges, fragments)
        else:
            self.remove()
            self.nr_fragments = len(fragments)
            self._create_visualization(atomlist, charges, fragments)
    def remove(self):
        for o in (self.charge_labels + self.frag_charge_labels + self.charge_clouds + self.dipole_vectors + self.frag_enclosing_boxes + self.sas_points + self.nac_vectors):
            o.remove()
        self.charge_labels = []
        self.frag_charge_labels = []
        self.charge_clouds = []
        self.dipole_vectors = []
        self.frag_enclosing_boxes = []
        self.sas_points = []
        self.nac_vectors = []
    @disable_rendering
    def _create_visualization(self, atomlist, charges, fragments):
        mlab = self.scene.mlab
        
        if self.show_flags["charges"] == True and type(charges) != type(None):
            self.charge_labels = []
            for q,(Z,pos) in zip(charges, atomlist):
                x,y,z = pos
                if q < 0.0:
                    color = (0.,0.,1.)
                elif q > 0.0:
                    color = (1.,0.,0.)
                else:
                    color = (1.,1.,1.)
                # color does not work so far
                txt = "%+2.3f" % q
                if self.show_flags["labels"] == True:
                    # maybe charges should not overlap with labels
                    txt = "%s" % txt
                label = mlab.text(x,y, txt, z=z, figure=self.scene.mayavi_scene)
                label.actor.set(text_scale_mode='none', width=0.05, height=0.1)
                label.property.set(justification='centered', vertical_justification='centered')
                self.charge_labels.append( label )
            self.shown_charges = True
        else:
            self.shown_charges = False
            
        if self.show_flags["frag. charges"] == True and type(charges) != type(None):
            self.frag_charge_labels = []
            for ifrag,(fragment_indeces, fragment_atomlist, fragment_box) in enumerate(fragments):
                # compute the charges on fragment 
                qfrag = np.sum(charges[fragment_indeces])
                print "Fragment charges: %s" % charges[fragment_indeces]
                # compute the center of the molecule,
                pos_frag = XYZ.atomlist2vector(fragment_atomlist)
                masses_frag = AtomicData.atomlist2masses(fragment_atomlist)
                com = MolCo.center_of_mass(masses_frag, pos_frag)
                #
                print "Fragment %d  charge = %s" % (ifrag, qfrag)
                txt = "%+2.3f" % qfrag
                label = mlab.text(com[0],com[1], txt, z=com[2], line_width=0.8, figure=self.scene.mayavi_scene)
                label.actor.set(text_scale_mode='none', width=0.05, height=0.1)
                label.property.set(justification='centered', vertical_justification='centered')
                self.frag_charge_labels.append( label )
            self.shown_frag_charges = True
        else:
            self.shown_frag_charges = False
            
        if self.show_flags["charge clouds"] == True and type(charges) != type(None):
            self.charge_clouds = []

            vec = XYZ.atomlist2vector(atomlist)
            x, y, z = vec[::3], vec[1::3], vec[2::3]
            s = abs(charges)

            # The charge clouds represent surfaces of equal charge density around each atoms.
            # In DFTB the charge fluctuations are modelled by a Gaussian:
            #  F(r) = 1/(2*pi*sA^2)^(3/2) * exp(-r^2/(2*sA^2))
            # The radii of the charge clouds are scaled by the charge on the atom:
            #  r = q * r0
            # The radius r0 belongs to a charge cloud containing exactly 1 electron, it depends
            # on the hubbard parameter through sA and on the isoValue:
            #  F(r0) = F(0) * isoValue
            #
            r0s = self.charge_cloud_radii.getRadii(atomlist)
            s *= r0s
            
            cloud = mlab.quiver3d(x,y,z,s,s,s, scalars=charges, 
                                           mode="sphere", scale_factor=1.0, resolution=20,
                                           opacity = 0.4, 
                                           figure=self.scene.mayavi_scene)

            # atoms are coloured by their atomic number
            cloud.glyph.color_mode = "color_by_scalar"
            cloud.glyph.glyph_source.glyph_source.center = [0,0,0]
            self.lut = cloud.module_manager.scalar_lut_manager.lut.table.to_array()
            red = np.array((255.0, 0.0, 0.0)).astype('uint8')
            blue = np.array((0.0, 0.0, 255.0)).astype('uint8')
            for i in range(0, 255):
                if i < 128:
                    color = blue
                else:
                    color = red
                self.lut[i,0:3] = color
            cloud.module_manager.scalar_lut_manager.lut.table = self.lut
            cloud.module_manager.scalar_lut_manager.data_range = (-1.0, 1.0)

            self.cloud = cloud
            self.charge_clouds.append( cloud )
            
            self.shown_charge_clouds = True
            
        else:
            self.shown_charge_clouds = False

        if self.show_flags["dipole moment"] == True and type(charges) != type(None):
            self.dipole_vectors = []
            # The dipole vector is placed at the center of mass
            pos = XYZ.atomlist2vector(atomlist)
            masses = AtomicData.atomlist2masses(atomlist)
            com = MolCo.center_of_mass(masses, pos)
            # compute dipole moment from charge distribution
            dipole = np.zeros(3)
            for i,q in enumerate(charges):
                dipole += q * pos[3*i:3*(i+1)]
            print "Dipole moment  D = %s a.u." % dipole
            # For plotting the dipole vector is converted to Debye
            dipole *= AtomicData.ebohr_to_debye
            print "Dipole moment  D = %s Debye" % dipole
            print "Length of dipole moment |D| = %s Debye" % la.norm(dipole)
            
            quiver3d = mlab.quiver3d(com[0],com[1],com[2],
                                     dipole[0], dipole[1], dipole[2], line_width=5.0,
                                     scale_mode='vector',
                                     color=(0,0,1), 
                                     scale_factor=1.0)
            
            self.dipole_vectors.append(quiver3d)
        else:
            self.shown_dipole_moment = False    
            
        if self.show_flags["enclosing box"] == True:
            self.frag_enclosing_boxes = []
            for ifrag,(fragment_indeces, fragment_atomlist, fragment_box) in enumerate(fragments):
                box = fragment_box
                # plot edges of the enclosing box
                for edge in box.edges:
                    l = mlab.plot3d(box.vertices[edge,0], box.vertices[edge,1], box.vertices[edge,2],
                                    color=(1,0,0),
                                    figure=self.scene.mayavi_scene)
                    self.frag_enclosing_boxes.append(l)
                # plot axes
                for axis, color in zip(box.axes, [(1.,0.,0.),(0.,1.,0.), (0.,0.,1.)]):
                    ax = mlab.quiver3d(float(box.center[0]), float(box.center[1]), float(box.center[2]),
                                       float(axis[0]), float(axis[1]), float(axis[2]),
                                       color=color, scale_factor=3.0, mode='arrow', resolution=20,
                                       figure=self.scene.mayavi_scene)
                    self.frag_enclosing_boxes.append(ax)
            self.shown_enclosing_boxes = True
        else:
            self.shown_enclosing_boxes = False    

        if self.show_flags["screening charges (COSMO)"] == True:
            # read parameter for solvent model from command line
            # or configuration file
            parser = OptionParserFuncWrapper([ImplicitSolvent.SolventCavity.__init__], "",
                                             unknown_options="ignore")
            # extract only arguments for constructor of solvent cavity
            (solvent_options, args) = parser.parse_args(ImplicitSolvent.SolventCavity.__init__)
            # 
            cavity = ImplicitSolvent.SolventCavity(**solvent_options)
            cavity.constructSAS(atomlist)
            area = cavity.getSurfaceArea()
            print "solvent accessible surface area: %8.6f bohr^2   %8.6f Ang^2" % (area, area*AtomicData.bohr_to_angs**2)
            points = cavity.getSurfacePoints()
            x,y,z = points[:,0], points[:,1], points[:,2]
            if type(charges) != type(None):
                # If there are Mulliken charges we can compute the
                # induced charges on the surface of the cavity
                # according to COSMO
                cavity.constructCOSMO()
                induced_charges = cavity.getInducedCharges(charges)
                screening_energy = cavity.getScreeningEnergy(charges)
                print "screening energy: %10.6f Hartree  %10.6f kcal/mol" \
                    % (screening_energy, screening_energy * AtomicData.hartree_to_kcalmol)
                # The surface points are colored and scaled
                # according to their charge
                #  negative -> blue   positive -> red
                points3d = mlab.points3d(x,y,z,induced_charges,
                                         colormap="blue-red",
                                         mode='2dcross',
                                         scale_factor=10)
#                                         mode='2dvertex')
                points3d.glyph.color_mode = "color_by_scalar"
            else:
                # 
                points3d = mlab.points3d(x,y,z,color=(0.0,0.0,0.8), mode='2dvertex')
            self.sas_points.append( points3d )
        else:
            self.shown_solvent_screening_charges = False

        if self.show_flags["non-adiab. coupling vectors"] == True:
            vec = XYZ.atomlist2vector(atomlist)
            x, y, z = vec[::3], vec[1::3], vec[2::3]
            
            # assume that charges are transition charges
            transition_charges = charges
            
            # Because we don't know the energy here, the energy difference
            # in the denominator is set to 1, so only the direction and relative
            # lengths are correct.
            nac = NACsApprox.coupling_vector(atomlist, transition_charges, 1.0)
            # directions of NAC vectors
            u,v,w = nac[0,:], nac[1,:], nac[2,:]
            # Add a NAC vector at each atom
            quiver3d = mlab.quiver3d(x,y,z, u,v,w, line_width=5.0)
            
            self.nac_vectors.append(quiver3d)
            
        else:
            self.shown_nacs = False
            
    def _update_visualization(self, atomlist, charges, fragments):
        """
        only change the underlying data but do not recreate visualization.
        This is faster and avoids jerky animations.
        """
        mlab = self.scene.mlab
        if self.show_flags["charges"] == True:
            i = 0
            for q,(Z,pos) in zip(charges, atomlist):
                x,y,z = pos
                if q < 0.0:
                    color = (0.,0.,1.)
                elif q > 0.0:
                    color = (1.,0.,0.)
                else:
                    color = (1.,1.,1.)
                # color does not work so far
                txt = "%+2.3f" % q
                if self.show_flags["labels"] == True:
                    # maybe charges should not overlap with labels
                    txt = "%s" % txt
                # update label position and text
                label = self.charge_labels[i]
                label.remove()
                label = mlab.text(x,y, txt, z=z, figure=self.scene.mayavi_scene)
                self.charge_labels[i] = label
                label.actor.set(text_scale_mode='none', width=0.05, height=0.1)
                label.property.set(justification='centered', vertical_justification='centered')

                i += 1
        if self.show_flags["frag. charges"] == True:
            i = 0
            for ifrag,(fragment_indeces, fragment_atomlist) in enumerate(fragments):
                # compute the charges on fragment 
                qfrag = np.sum(charges[fragment_indeces])
                print "Fragment charges: %s" % charges[fragment_indeces]
                # compute the center of the molecule,
                pos_frag = XYZ.atomlist2vector(fragment_atomlist)
                masses_frag = AtomicData.atomlist2masses(fragment_atomlist)
                com = MolCo.center_of_mass(masses_frag, pos_frag)
                #
                print "Fragment %d  charge = %s" % (ifrag, qfrag)
                txt = "%+2.3f" % qfrag
                label = self.frag_charge_labels[i]
                label.remove()
                label = mlab.text(com[0],com[1], txt, z=com[2], line_width=0.8, figure=self.scene.mayavi_scene)
                self.frag_charge_labels[i] = label
                label.actor.set(text_scale_mode='none', width=0.05, height=0.1)
                label.property.set(justification='centered', vertical_justification='centered')

                i += 1
        if self.show_flags["charge clouds"] == True:

            vec = XYZ.atomlist2vector(atomlist)
            x, y, z = vec[::3], vec[1::3], vec[2::3]
            s = abs(charges)

            # The charge clouds represent surfaces of equal charge density around each atoms.
            # In DFTB the charge fluctuations are modelled by a Gaussian:
            #  F(r) = 1/(2*pi*sA^2)^(3/2) * exp(-r^2/(2*sA^2))
            # The radii of the charge clouds are scaled by the charge on the atom:
            #  r = q * r0
            # The radius r0 belongs to a charge cloud containing exactly 1 electron, it depends
            # on the hubbard parameter through sA and on the isoValue:
            #  F(r0) = F(0) * isoValue
            #
            r0s = self.charge_cloud_radii.getRadii(atomlist)
            s *= r0s

            self.cloud.mlab_source.set(x=x,y=y,z=z,u=s,v=s,w=s, scalars=charges, scale_factor=1.0)
            # atoms are coloured by their atomic number
            self.cloud.glyph.color_mode = "color_by_scalar"
            self.cloud.glyph.glyph_source.glyph_source.center = [0,0,0]
            self.cloud.module_manager.scalar_lut_manager.lut.table = self.lut
            self.cloud.module_manager.scalar_lut_manager.data_range = (-1.0, 1.0)

            
#The actual visualization
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())
    
    # list of geometries that are displayed as ball and stick models
    molecules = []
    # list of volumetric data sets
    cubes = []
    # list of lists of charges
    charge_lists = []
    # isovalue
    isoValue = 0.0004
    # list of drawables
    cube_objects = []
    molecule_objects = []
    charge_objects = []
    label_objects = []
    #
    have_fragments = False
    fragment_lists = []
    #
    animation = False
    #
    show_flags = {"atoms": True, "labels": False, "Lewis structure": False, "isosurface": True,
                  "X cut plane": False, "Y cut plane": False, "Z cut plane": False,
                  "charges": False, "frag. charges": False, "charge clouds": False, "dipole moment": False, 
                  "enclosing box": False, "screening charges (COSMO)": False, "non-adiab. coupling vectors": False}
    # charge clouds
    charge_cloud_radii = ChargeCloudRadii()
    def setFlag(self, flagname):
        self.show_flags[flagname] = True
    def unsetFlag(self, flagname):
        self.show_flags[flagname] = False
    def setAnimation(self, flag):
        self.animation = flag
    def setGeometries(self, molecules):
        self.molecules = molecules
        self.have_fragments = False
        self.fragment_lists = [[] for atomlist in self.molecules]
        self.charge_lists = [None for atomlist in self.molecules]
    def setVolumetricData(self, cubes):
        self.cubes = cubes
    def setIsoValue(self, isoValue):
        self.isoValue = isoValue
        self.charge_cloud_radii.setIsoValue(isoValue)
    def setCharges(self, charge_lists):
        self.charge_lists = charge_lists
    def calcFragments(self):
        self.fragment_lists = []
        for atomlist in self.molecules:
            fragment_graphs = MolecularGraph.atomlist2graph(atomlist)
            fragments = []
            for g in fragment_graphs:
                fragment_indeces = MolecularGraph.graph2indeces(g)
                fragment_atomlist = MolecularGraph.graph2atomlist(g, atomlist)
                fragment_box = MoleculeBox(fragment_atomlist)
                fragments.append( (fragment_indeces, fragment_atomlist, fragment_box) )
            self.fragment_lists.append(fragments)
        self.have_fragments = True
    @disable_rendering        
    def plotCube(self, cube, iso=0.04):
        #print "visualize cubes"
        #
        mlab = self.scene.mlab
        x,y,z = cube.grid
        amplitude = mlab.pipeline.scalar_field(x,y,z, cube.data, figure=self.scene.mayavi_scene)
        self.cube_objects += [amplitude]
        
        if self.show_flags["isosurface"] == True:
            # red position lobe
            positive_lobe = mlab.contour3d(x,y,z,cube.data, contours=[+iso],
                                           color=(1.0, 0.0, 0.0),
                                           opacity=0.5,
                                           figure=self.scene.mayavi_scene)
            # blue negative lobe
            negative_lobe = mlab.contour3d(x,y,z,cube.data, contours=[-iso],
                                           color=(0.0, 0.0, 1.0),
                                           opacity=0.5,
                                           figure=self.scene.mayavi_scene)
            #
            lobes = [positive_lobe, negative_lobe]
            for lobe in lobes:
                lobe.actor.property.representation = "surface" #"wireframe"
            self.cube_objects += [positive_lobe, negative_lobe] 
            
        nx,ny,nz = cube.data.shape
        cuts = []
        if self.show_flags["X cut plane"] == True:
            xcut = mlab.pipeline.image_plane_widget(amplitude, plane_orientation='x_axes', slice_index=nx/2, opacity=0.5, figure=self.scene.mayavi_scene)
            cuts.append(xcut)
        if self.show_flags["Y cut plane"] == True:
            ycut = mlab.pipeline.image_plane_widget(amplitude, plane_orientation='y_axes', slice_index=ny/2, opacity=0.5, figure=self.scene.mayavi_scene)
            cuts.append(ycut)
        if self.show_flags["Z cut plane"] == True:
            zcut = mlab.pipeline.image_plane_widget(amplitude, plane_orientation='z_axes', slice_index=nz/2, opacity=0.5, figure=self.scene.mayavi_scene)
            cuts.append(zcut)

        for cut in cuts:
            cut.module_manager.scalar_lut_manager.data_range = [-iso, iso]

        self.cube_objects += cuts
        #mlab.outline()
    @disable_rendering        
    def plotLabels(self):
        #print "visualize labels"
        # draw atom labels
        if self.show_flags["labels"] == True:
            mlab = self.scene.mlab
            labels = []
            molecules = []
            for cube in self.cubes:
                molecules.append(cube.atomlist)
            for atomlist in self.molecules:
                molecules.append(atomlist)    
            for atomlist in molecules:
                nat = len(atomlist)
                vec = XYZ.atomlist2vector(atomlist)
                x, y, z = vec[::3], vec[1::3], vec[2::3]
                atom_names = [AtomicData.atom_names[Z-1] for (Z,pos) in atomlist]
                for i in range(0, nat):
                    name = atom_names[i].upper() + "-%d" % i
                    label = mlab.text(x[i],y[i],name, z=z[i], figure=self.scene.mayavi_scene)
                    label.actor.set(text_scale_mode='none', width=0.05, height=0.1)
                    label.property.set(justification='centered', vertical_justification='centered')
                    labels.append(label)
            self.label_objects += labels
    def clear(self):
        for o in (self.molecule_objects + self.cube_objects + self.charge_objects + self.label_objects):
            o.remove()
        self.molecule_objects = []
        self.cube_objects = []
        self.charge_objects = []
        self.label_objects = []
    @on_trait_change('scene.activated')
    def update_all(self):
        # This function is called when the view is opened. We don't
        # populate the scene when the view is not yet open, as some
        # VTK features require a GLContext.
        self.scene.foreground = (0.0, 0.0, 0.0)  # black text
        self.scene.background = (1.0, 1.0, 1.0)  # white background
        self.update_molecules()
        self.update_cubes()
        self.update_charges()
        self.update_labels()
        
    def update_molecules(self):
        if self.animation == True and len(self.molecule_objects) > 0:
            for molob, atomlist in zip(self.molecule_objects, self.molecules):
                molob.draw(atomlist, animation=True)
        else:
            for o in self.molecule_objects:
                o.remove()
            self.molecule_objects = []
            for atomlist in self.molecules:
                molob = MoleculeObject(self.scene, self.show_flags)
                molob.draw(atomlist, animation=False)
                self.molecule_objects.append(molob)
    def update_cubes(self):
        for o in self.cube_objects:
            o.remove()
        self.cube_objects = []

        for cube in self.cubes:
            self.plotCube(cube, self.isoValue)
    def update_charges(self):
        if self.animation == True and len(self.charge_objects) > 0:
            # compute fragments if necessary
            if ((self.show_flags["frag. charges"] == True or self.show_flags["enclosing box"] == True) and self.have_fragments == False):
                self.calcFragments()
            for chrgob, atomlist, charges, fragments in zip(self.charge_objects, self.molecules, self.charge_lists, self.fragment_lists):
                chrgob.draw(atomlist, charges, fragments, animation=True)
        else:
            for o in self.charge_objects:
                o.remove()
            self.charge_objects = []
            # compute fragments if necessary
            if ((self.show_flags["frag. charges"]) == True or self.show_flags["enclosing box"] == True) and (self.have_fragments == False):
                self.calcFragments()
            for atomlist,charges,fragments in zip(self.molecules, self.charge_lists, self.fragment_lists):
                chrgob = ChargesObject(self.scene, self.show_flags, self.charge_cloud_radii)
                chrgob.draw(atomlist, charges, fragments, animation=False)
                self.charge_objects.append(chrgob)
    def update_labels(self):
        for o in self.label_objects:
            o.remove()
        self.label_objects = []

        self.plotLabels()

    # the layout of the dialog is created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True # We need this to resize with the parent widget
                )

################################################################################
# The QWidget containing the visualization, this is pure PyQt4 code.

class QCubeViewerWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)

        # field for setting iso value
        isoFrame = QtGui.QFrame()
        layout.addWidget(isoFrame)
        isoLayout = QtGui.QHBoxLayout(isoFrame)
        self.isoLabel = QtGui.QLabel("Iso Value:")
        isoLayout.addWidget(self.isoLabel)
        self.isovalue = QtGui.QLineEdit()
        self.isovalue.setText(str(self.visualization.isoValue))
        self.isovalue.returnPressed.connect(self.redraw)
        isoLayout.addWidget(self.isovalue)
        
        self.showButton = QtGui.QPushButton("Show...")
        self.menu = QtGui.QMenu()
        self.options = {}
        # general options
        self.options["atoms"] = QtGui.QAction("atoms", self.menu, checkable=True, checked=True)
        self.options["labels"] = QtGui.QAction("labels", self.menu, checkable=True)
        self.options["Lewis structure"] = QtGui.QAction("Lewis structure", self.menu, checkable=True)
        self.options["enclosing box"] = QtGui.QAction("enclosing box", self.menu, checkable=True, checked=False)
        
        self.menu.addAction(self.options["atoms"])
        self.menu.addAction(self.options["labels"])
        self.menu.addAction(self.options["Lewis structure"])
        self.menu.addAction(self.options["enclosing box"])

        self.menu.addSeparator()
        # options for volumetric data
        self.options["isosurface"] = QtGui.QAction("isosurface", self.menu, checkable=True, checked=True)
        self.options["X cut plane"] = QtGui.QAction("X cut plane", self.menu, checkable=True)
        self.options["Y cut plane"] = QtGui.QAction("Y cut plane", self.menu, checkable=True)
        self.options["Z cut plane"] = QtGui.QAction("Z cut plane", self.menu, checkable=True)

        self.menu.addAction(self.options["isosurface"])
        self.menu.addAction(self.options["X cut plane"])
        self.menu.addAction(self.options["Y cut plane"])
        self.menu.addAction(self.options["Z cut plane"])

        # options for charges
        self.options["charges"] = QtGui.QAction("charges", self.menu, checkable=True, checked=True)
        self.options["frag. charges"] = QtGui.QAction("frag. charges", self.menu, checkable=True, checked=False)
        self.options["charge clouds"] = QtGui.QAction("charge clouds", self.menu, checkable=True, checked=False)
        self.options["dipole moment"] = QtGui.QAction("dipole moment", self.menu, checkable=True, checked=False)
        self.options["screening charges (COSMO)"] = QtGui.QAction("screening charges (COSMO)", self.menu, checkable=True)
        self.options["non-adiab. coupling vectors"] = QtGui.QAction("non-adiab. coupling vectors", self.menu, checkable=True)
        
        self.menu.addSeparator()
        self.menu.addAction(self.options["charges"])
        self.menu.addAction(self.options["frag. charges"])
        self.menu.addAction(self.options["charge clouds"])
        self.menu.addAction(self.options["dipole moment"])
        self.menu.addAction(self.options["screening charges (COSMO)"])
        self.menu.addAction(self.options["non-adiab. coupling vectors"])
        
        for action in self.options.values():
            action.triggered.connect(self.redraw)
        self.showButton.setMenu(self.menu)
        # by default only the options for surfaces are shown
        self.selectShowOptions(options=["surfaces"])
        
        isoLayout.addWidget(self.showButton)
        # fullscreen
        icon = QtGui.QIcon.fromTheme("view-fullscreen")
        #print "icon = %s" % icon
        self.fullscreenButton = QtGui.QPushButton(icon, 'fullscreen on/off')
        isoLayout.addWidget(self.fullscreenButton)
        self.fullscreen = False
        self.fullscreenButton.clicked.connect(self.toggle_fullscreen)
                
    def toggle_fullscreen(self):
        self.fullscreen = not self.fullscreen
        if self.fullscreen == True:
            print "fullscreen: ON"
            self.setWindowFlags(QtCore.Qt.Window)
            self.setWindowState(QtCore.Qt.WindowFullScreen)
            self.show()
        else:
            print "fullscreen: OFF"
            self.setWindowState(QtCore.Qt.WindowNoState)
            self.setWindowFlags(QtCore.Qt.Widget)
            self.show()

    def getShowOptions(self):
        for k,v in self.options.iteritems():
            if v.isChecked():
                self.visualization.setFlag(k)
            else:
                self.visualization.unsetFlag(k)
                
    def setCubes(self, cubes):
        molecules = [cube.atomlist for cube in cubes]
        self.visualization.setGeometries(molecules)
        self.visualization.setVolumetricData(cubes)
        self.getShowOptions()
        self.visualization.update_all()

    def setAnimation(self, flag):
        self.visualization.setAnimation(flag)

    def savefig(self, filename):
        self.visualization.scene.save(filename)
        
    def clear(self):
        self.visualization.clear()
        
    def setGeometries(self, molecules):
        self.visualization.setGeometries(molecules)
        self.getShowOptions()
        self.visualization.update_all()
        
    def setGeometriesAndCharges(self, molecules,charge_lists):
        self.visualization.setGeometries(molecules)
        self.visualization.setCharges(charge_lists)
        self.getShowOptions()
        self.visualization.update_all()

    def selectShowOptions(self, options=["surfaces", "charges"]):
        """select which options should be shown when the 'Show...' menu is opened"""
        # hide all options
        for o in self.options.values():
            o.setVisible(False)
        # activate selected option blocks
        self.options["atoms"].setVisible(True)
        self.options["labels"].setVisible(True)
        self.options["Lewis structure"].setVisible(True)
        if "surfaces" in options:
            self.options["isosurface"].setVisible(True)
            self.options["X cut plane"].setVisible(True)
            self.options["Y cut plane"].setVisible(True)
            self.options["Z cut plane"].setVisible(True)
        if "charges" in options:
            self.options["charges"].setVisible(True)
            self.options["frag. charges"].setVisible(True)
            self.options["charge clouds"].setVisible(True)
            self.options["dipole moment"].setVisible(True)        
            self.options["screening charges (COSMO)"].setVisible(True)
            self.options["enclosing box"].setVisible(True)
        if "transition charges" in options:
            self.options["non-adiab. coupling vectors"].setVisible(True)
            
    def hideIsoControls(self):
        self.isoLabel.hide()
        self.isovalue.hide()

    def redraw(self):
        self.getShowOptions()
        isoValue = abs(float(str(self.isovalue.text())))
        self.visualization.setIsoValue(isoValue)

        if self.options["charges"].isVisible() == True \
           or self.options["Lewis structure"].isVisible() == True:
            # this is a HACK, charges and Lewis structures are not drawn
            # unless the molecule is updated too
            self.visualization.update_molecules()
        self.visualization.update_cubes()
        self.visualization.update_charges()
        self.visualization.update_labels()

        #self.visualization.update_all()
        
    def setIsoValue(self, isoValue):
        self.visualization.setIsoValue(isoValue)
        self.isovalue.setText(str(self.visualization.isoValue))
        
if __name__ == "__main__":
    import sys
    import os.path
    if len(sys.argv) < 2:
        print " "
        print "Usage: %s <list of cube- or xyz-files>" % os.path.basename(sys.argv[0])
        print " "
        print "  displays molecular geometries and volumetric data using Mayavi"
        print " "
        exit(-1)

    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()

    window = QtGui.QMainWindow()
    viewer = QCubeViewerWidget(window)

    files = sys.argv[1:]
    cubes = []
    molecules = []
    charge_lists = []
    for f in files:
        if ".cub" in f:
            cube = CubeData()
            cube.loadFromFile(f)
            cubes.append(cube)
            molecules.append(cube.atomlist)
        elif ".xyz" in f:
            atomlist = XYZ.read_xyz(f)[-1]
            molecules.append(atomlist)
        elif ".chg" in f:
            atomlist, charges = XYZ.read_charges(f)
            molecules.append(atomlist)
            charge_lists.append(charges)
        else:
            print "Unknown file type: %s" % f
        window.setWindowTitle(window.windowTitle() + " %s " % f)
    viewer.setCubes(cubes)
    viewer.setGeometries(molecules)
    
    if len(charge_lists) > 0:
        viewer.setGeometriesAndCharges(molecules, charge_lists)
        viewer.selectShowOptions(options=["surfaces", "charges", "transition charges"])
        
    """
    # plot E-field
    mlab = viewer.visualization.scene.mlab
    mlab.quiver3d(0,0,0, 0,0,10.0, color=(0.0,1.0,0.0), scale_factor=1.0, mode='arrow', resolution=20)
    mlab.text(0.0, 0.0,"E-field",z=10.0) 
    """

    window.setCentralWidget(viewer)
    window.show()
    
    # Start the main event loop.
    app.exec_()

    print "FINISHED"
