# NOT FINISHED YET
"""
This script analyses the inter- and intra-molecular redistribution of vibrational energy. 
In a trajectory the molecular fragments are identified and the kinetic energies
of the centers of mass (COM) are put into relation with the total kinetic energy. The ratio of the
two types of energies shows how directed motion is converted into vibrational 'random' motion.
"""

from DFTB.Analyse import MolecularGraph
from DFTB.Modeling import MolecularCoords
from DFTB import XYZ, AtomicData

import numpy as np

def velocities_finite_diff(atomlists, dt):
    """
    compute the velocities from finite differences between two time-steps

    Parameters:
    ===========
    atomlists: list of geometries along trajectories for each time step
    dt: nuclear time step in a.u.

    Returns:
    ========
    positions: list of vectors with positions for each time step
    velocities: list of vectors with velocities for each time step
    """
    Nt = len(atomlists)
    positions = []
    velocities = []
    pos0 = XYZ.atomlist2vector(atomlists[0])
    positions.append(pos0)
    for i in range(1, Nt):
        pos1 = XYZ.atomlist2vector(atomlists[i])
        vel = (pos1-pos0)/dt
        positions.append(pos1)
        velocities.append(vel)
        pos0 = pos1
    velocities.append(vel)

    return positions, velocities
        
def partition_kinetic_energy(atomlists, dt):
    Nt = len(atomlists)
    # determine the indeces of the molecular fragments at the first
    # time-step. 
    fragments_graph = MolecularGraph.atomlist2graph(atomlists[0])    
    fragments_indeces = [MolecularGraph.graph2indeces(g) for g in fragments_graph]
    fragment_atomlists = [[ [atomlist[ii] for ii in I] for atomlist in atomlists] for I in fragments_indeces]
    # compute vibrational and center of mass kinetic energy for each fragment
    Nfrag = len(fragment_atomlists)
    ekin_com = np.zeros((Nt,Nfrag))
    ekin_tot = np.zeros((Nt,Nfrag))
    for f in range(0, Nfrag):
        masses = AtomicData.atomlist2masses(fragment_atomlists[f][0])
        print fragment_atomlists[f][0]
        # total mass
        M = np.sum(masses)/3.0
        positions, velocities = velocities_finite_diff(fragment_atomlists[f], dt)
        for i in range(0, Nt):
            vm = velocities[i]*masses
            vel_com = np.array([np.sum(vm[0::3]), np.sum(vm[1::3]), np.sum(vm[2::3])])/M
            ekin_com_f = 1.0/2.0 * M * np.sum(vel_com**2)
            ekin_tot_f = 1.0/2.0 * np.sum(masses * velocities[i]**2)

            ekin_com[i,f] += ekin_com_f
            ekin_tot[i,f] += ekin_tot_f
    return ekin_com, ekin_tot


if __name__ == "__main__":
    import sys
    traj_file = sys.argv[1]
    dt = float(sys.argv[2])
    atomlists = XYZ.read_xyz(traj_file)
    ekin_com, ekin_tot = partition_kinetic_energy(atomlists, 3.0)

    print ekin_com[:100,0]
    print ekin_com[:100,1]
    import matplotlib.pyplot as plt

    for f in range(0,2):
        plt.plot(ekin_com[:,f]/ekin_com[:,f].max(), ls="-.")
    plt.show()
#    fragtraj, fragdic = MolecularGraph.fragment_trajectory(atomlists[::100])
#    print fragtraj
#    print fragdic
