"""
generate a set of different nanotubes
"""
from DFTB.Modeling.nanotube_builder import capped_uff_nanotube
from os.path import join

L = 200.0

#
d="~/HERTEL/NANOTUBES/KATAURA/geometries/"

n1 = 8
n2 = 2
capped_uff_nanotube(n1,n2,L,out_xyz=join(d,"cnt_%d-%d.xyz" % (n1,n2)))

#
"""
for n1 in range(1,10):
    for n2 in range(1,10):
        if 8 <= n1+n2 <= 16:
            capped_uff_nanotube(n1,n2,L,
                                out_xyz="~/HERTEL/NANOTUBES/KATAURA/geometries/cnt_%d-%d.xyz" % (n1,n2))
"""
