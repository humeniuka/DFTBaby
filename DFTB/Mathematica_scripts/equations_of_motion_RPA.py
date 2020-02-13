import sympy
from sympy.physics.secondquant import *

p,q,r,s = sympy.symbols("p q r s")

Spq = Fd(p)*F(q)
Srs = Fd(r)*F(s)
com = Commutator(Spq, Srs)
print com.doit()

