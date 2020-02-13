"""
The Dyson orbitals of neutral water clusters can be approximated as linear combinations of 
MOs of water on different monomers. This script builds the cluster geometries from templates
and constructs the cluster orbitals. 

"""
import numpy as np

from DFTB import XYZ, AtomicData
from DFTB.Analyse import MolecularGraph, Cube
from DFTB.BasisSets import load_pseudo_atoms, AtomicBasisSet
from DFTB.Modeling import MolecularCoords as MolCo
from DFTB.Modeling import OrbitalRotations
from DFTB.Scattering.SlakoScattering import save_dyson_orbitals, load_dyson_orbitals
from os.path import join

def combine_oriented_fragments(frag_orientations):
    atomlist_combined = []
    for atomlist_std, (a,b,g), cm in frag_orientations:
        atomlist_combined += OrbitalRotations.transform_molecule(atomlist_std, (a,b,g), cm)
    return atomlist_combined


def build_water_cluster(atomlist_template, orbital_names, phases):
    """
    build a cluster orbital as a linear combination of monomer orbitals with the correct orientation.
    The position and orientation of the water molecules is taken from the template. 
    
    Parameters:
    ===========
    atomlist_template: molecular geometry for cluster with n water molecules
    orbital_names: list of orbital names ('1b1', '3a1', '1b2' or '2a1') for each of the n water molecules
    phases: list of +1 or -1's for each orbital

    Returns:
    ========
    water_cluster: molecular geometry of the cluster
    orb_cluster: cluster orbital, that is a linear combination of the monomer orbitals.
    """
    # water geometry in bohr
    water_std = [(1, (0.0,  0.8459947982381987,  -1.4473477675908173)),
                 (8, (0.0, -0.21392561795490195,  0.0               )),
                 (1, (0.0,  0.8459947982381987,   1.4473477675908173))]
    valorbs, radial_val = load_pseudo_atoms(water_std)

    # Orbitals for H2O monomer in standard orientation:
    # molecular plane = yz-plane, oxygen lies on the negative y-axis, H-H bond is parallel to z-axis
        #                  H1-1s   O-2s    O-2py   O-2pz   O-2px   H2-1s  
    monomer_orbitals = {
        "1b1": np.array([ 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000]),
        "3a1": np.array([ 0.7816,-0.6842, 0.6989, 0.0000, 0.0000, 0.7816]),
        "1b2": np.array([-0.7264, 0.0000, 0.0000, 0.9429, 0.0000, 0.7264]),
        "2a1": np.array([ 0.1730, 0.8662, 0.0393, 0.0000, 0.0000, 0.1730])
    }

    # First the individual water molecules in the template are identified and their positions
    # and orientations are extracted.
    fragments = MolecularGraph.disconnected_fragments(atomlist_template)
    assert len(fragments) == len(orbital_names) == len(phases), "For each water fragment in the cluster you need to specify one fragment orbital ('1b1', '3a1', '1b2' or '2a1') and its phase"
    orb_cluster = []    # list of MO coefficients
    water_cluster = []  # list of atoms
    for i,water in enumerate(fragments):
        # the Euler angles (a,b,g) specify the orientation of the i-th water molecule
        water_std_i, (a,b,g), cm = MolCo.molecular_frame_transformation(water)
        print "WATER STANDARD"
        for (Zi,posi) in water_std:
            print "  %s   %8.6f  %8.6f  %8.6f" % (Zi, posi[0], posi[1], posi[2])
        print "WATER STANDARD %d" % i
        for (Zi,posi) in water_std_i:
            print "  %s   %8.6f  %8.6f  %8.6f" % (Zi, posi[0], posi[1], posi[2])
        # The desired orbital is placed on the i-th water molecule and is
        # rotated to match the orientation of the molecule.
        orb_monomer = monomer_orbitals[orbital_names[i]]
        # rotate orbital
        orb_monomer_rot = OrbitalRotations.rotate_orbitals(water_std, valorbs, orb_monomer, (a,b,g))
        # add orbital with desired phase to cluster MOs
        orb_cluster += list(phases[i] * orb_monomer_rot)
        # rotate geometry
        water_rot = MolCo.transform_molecule(water_std, (a,b,g), cm)
        XYZ.write_xyz("/tmp/water_std.xyz", [water_std, water_rot])
        water_cluster += water_rot
    # Assuming that the overlap between orbitals on different water molecules is negligible,
    # we can normalize the cluster orbital by dividing through sqrt(number of water molecules)
    n = np.sum(abs(np.array(phases))) # a phase of 0 indicates no orbital
    orb_cluster /= np.sqrt(n)
    
    return water_cluster, orb_cluster

def build_monomer(output_dir):
    cluster_template = XYZ.read_xyz("water_monomer.xyz")[0]
    names = ["1b1", "3a1", "1b2", "2a1"]
    ionization_energies = [12.5, 14.8, 18.5, 32.6]
    orbitals_list = [["1b1"], ["3a1"], ["1b2"], ["2a1"]] 
    phases_list   = [[ +1  ], [ +1  ], [ +1  ], [ +1  ]]
                
    dyson_orbs = []
    for name,IE,orbitals,phases in zip(names, ionization_energies,orbitals_list,phases_list):
        cluster, orb = build_water_cluster(cluster_template, orbitals, phases)
        dyson_orbs.append(orb)

    # save geometry of cluster
    XYZ.write_xyz(join(output_dir, "water_cluster_1.xyz"), [cluster])
    # and dyson orbitals
    dyson_orbs = np.array(dyson_orbs).transpose()
    save_dyson_orbitals(join(output_dir, "water_cluster_1.dyson"), names, ionization_energies, dyson_orbs)
    

def build_dimer(output_dir):
    cluster_template = XYZ.read_xyz("water_dimer.xyz")[0]
    names = ["S0->D0", "S0->D1", "S0->D2", "S0->D3", "S0->D4", "S0->D5"]
    ionization_energies = [11.55, 11.86, 13.43, 13.98, 17.71, 18.03]
    orbitals_list = [["1b1", "1b1"], ["1b1", "1b1"], ["3a1", "3a1"], ["3a1", "3a1"], ["1b2", "1b2"], ["1b2", "1b2"]]
    phases_list   = [[ -1  ,   -1 ], [ +1,    -1  ], [ +1  ,  +1  ], [ -1  ,  +1  ], [  0,     +1 ], [  +1 ,  0   ]]
                
    dyson_orbs = []
    for name,IE,orbitals,phases in zip(names, ionization_energies,orbitals_list,phases_list):
        cluster, orb = build_water_cluster(cluster_template, orbitals, phases)
        dyson_orbs.append(orb)

    # save geometry of cluster
    XYZ.write_xyz(join(output_dir, "water_cluster_2.xyz"), [cluster])
    # and dyson orbitals
    dyson_orbs = np.array(dyson_orbs).transpose()
    save_dyson_orbitals(join(output_dir, "water_cluster_2.dyson"), names, ionization_energies, dyson_orbs)
    
import sys
if __name__ == "__main__":
    build_monomer("/local_scratch/humeniuka/PHOTOANGULAR_DISTRIBUTIONS/H2O_CLUSTERS/1mer")
    #build_dimer("/local_scratch/humeniuka/PHOTOANGULAR_DISTRIBUTIONS/H2O_CLUSTERS/2mer")
    
    
    """
    cluster_template = XYZ.read_xyz("water_dimer.xyz")[0]
    cluster, orb = build_water_cluster(cluster_template, ["1b1", "1b1"], [+1,+1])
    # save cluster and template to xyz-file
    XYZ.write_xyz("/tmp/water_cluster_2.xyz", [cluster_template, cluster])
    # save cluster MO to cube file
    bs = AtomicBasisSet(cluster)
    Cube.orbital2grid(cluster, bs.bfs, orb, filename="/tmp/water_cluster_mo.cube", ppb=3.0)
    # save cluster MO coefficients
    from DFTB.Scattering.SlakoScattering import save_dyson_orbitals, load_dyson_orbitals
    orbs = np.zeros((len(orb),1))
    orbs[:,0] = orb
    save_dyson_orbitals("/tmp/water_cluster_2.dyson", ["S0->D0"], [-10.0], orbs)
    names, IEs, orbs = load_dyson_orbitals("/tmp/water_cluster_2.dyson")
    print names
    print IEs
    print orbs
    """

    """
    if len(sys.argv) < 2:
        print "Usage: python %s <.xyz file with template>" % sys.argv[0]
        exit(-1)

    xyz_file = sys.argv[1]
    atomlist_orig = XYZ.read_xyz(xyz_file)[-1]

    fragments = MolecularGraph.disconnected_fragments(atomlist_orig)

    fragments_std = []
    for i,atomlist in enumerate(fragments):
        atomlist_std, (a,b,g), cm = MolCo.molecular_frame_transformation(atomlist)
        print "Fragment %d" % i
        print "center of mass: %s" % cm
        print "Euler angles: %s %s %s" % (a,b,g)
        print "Standard geometry:"
        for Zi,posi in atomlist_std:
            print "   %s   %8.6f %8.6f %8.6f" % (AtomicData.atom_names[Zi-1], posi[0], posi[1], posi[2])
        fragments_std.append(atomlist_std)
        XYZ.write_xyz("/tmp/fragment_std.xyz", [atomlist_std])
    XYZ.write_xyz("/tmp/fragments_water.xyz", fragments_std)
    
    atomlist_combined = combine_oriented_fragments(
        [MolCo.molecular_frame_transformation(atomlist) for atomlist in fragments])

    XYZ.write_xyz("/tmp/fragments_combined.xyz", [atomlist_orig, atomlist_combined])
    """
