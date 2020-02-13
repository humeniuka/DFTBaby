"""
Slater-Koster tables with overlaps and hamiltonian matrix elements as functions of the interatomic distance. The files <atom 1>_<atom 2>.py contain the Slater-Koster tables for the atom pair 1-2.
"""

"""
import h_h, h_c, h_n, h_o, h_ag, \
            c_c, c_n, c_o, c_ag, \
                 n_n, n_o, n_ag, \
                      o_o, o_ag, \
                          ag_ag

atompairs = {(1,1): h_h, (1,6): h_c, (1,7): h_n, (1,8): h_o, (1,47): h_ag,\
                          (6,6): c_c, (6,7): c_n, (6,8): c_o, (6,47): c_ag,\
                                      (7,7): n_n, (7,8): n_o, (7,47): n_ag,\
                                                  (8,8): o_o, (8,47): o_ag,\
                                                             (47,47):ag_ag}
"""
atompairs = {}

def load_atompairs(atpairs):
    """import all atom pairs that are needed automatically"""
    from DFTB.AtomicData import atom_names
    for (Zi,Zj) in atpairs:
        atompairs[(Zi,Zj)] = __import__("DFTB.SlaterKoster.slako_tables.%s_%s" % (atom_names[Zi-1], atom_names[Zj-1]), fromlist=['Z1','Z2'])

__all__ = ["atompairs"]
