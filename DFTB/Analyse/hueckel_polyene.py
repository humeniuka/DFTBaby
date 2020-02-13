"""
The shift of the lowest excited state in a wire of meso-meso,beta-beta,beta-beta fused
porphyrins can be explained with the Hueckel model for linear polyenes.

The energy of the first excited state roughly equals the HOMO-LUMO gap, which for
a linear conjugated polymer consisting of N monomer units is given by:

E(HOMO-LUMO) = -4 beta sin(pi/(2 (N+1) )) ~ -4 beta pi/2 1/(N+1)

where beta = <a|H|a+1> is the interaction between the relevant orbital on adjacent monomers.
"""
import numpy as np
import numpy.linalg as la
from matplotlib import pyplot as plt

from DFTB import AtomicData

# energy of linear wire with n units
en_max_eV = np.array([4.07, 3.23, 1.06, 0.86, 0.72, 0.62, 0.52, 0.48, 0.42, 0.38, 0.34, 0.32])
en_max=en_max_eV / AtomicData.hartree_to_eV
wavelength_max=AtomicData.hartree_to_nm/en_max
N = np.array(range(1, len(en_max)+1))

# fit a line m*x + c through the data
x = N
y = wavelength_max
A = np.vstack([x, np.ones(len(x))]).transpose()
m, c = la.lstsq(A, y)[0]

plt.xlabel("$N$", fontsize=15)
plt.ylabel("$\lambda_{max}$ / nm", fontsize=15)

plt.plot(N, wavelength_max, "o", label="")
plt.plot(x, m*x+c, ls="-", lw=2, label="%3.1f N + %3.1f nm" % (m,c))
plt.legend(loc="lower center")
plt.show() 
