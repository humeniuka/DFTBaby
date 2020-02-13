#!/usr/bin/env python
"""
Shows energies of exciton states as a function of the length
of a linear chain with Hamiltonian

 H = sum_(n=1)^N { E |n><n| - J (|n><n+1| + |n+1><n|) }

In the basis state |n>, the n-th monomer unit is excited while all others 
are in the ground state

The coupling J is a positive number and only acts between nearest neighbours.
"""
import numpy as np
import matplotlib.pyplot as plt

# excitation energy of a single monomer
en=1.0
# coupling
J=0.5*en

Nmax=20
N=np.linspace(1,Nmax,Nmax-2)

for j in range(1,4):
    Ej=en-2*J*np.cos((np.pi*j)/(N+1.0))
    plt.plot(N,Ej,lw=2,label="S%d" % j)

plt.legend()
plt.xlabel("N, length of chain")
plt.ylabel("Exc. energy")
plt.savefig("linear_excitonic_chain.png")
plt.show()
             
