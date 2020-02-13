# DUMMY REPULSIVE POTENTIAL
import numpy as np
Z1 = 0
Z2 = 0
# grid for distance d between atomic centers in bohr
d = np.linspace(0.6, 5.0, 200)
# dummy repulsive potential in hartree/bohr, is zero everywhere
Vrep = np.zeros(200)
