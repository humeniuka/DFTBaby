# DUMMY REPULSIVE POTENTIAL for H-Cl
import numpy as np
Z1 = 1
Z2 = 17
# grid for distance d between atomic centers in bohr
d = np.linspace(0.6, 5.0, 200)
# repulsive potential in hartree/bohr
Vrep = np.zeros(200)
