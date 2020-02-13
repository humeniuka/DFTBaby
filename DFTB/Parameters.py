
"""
Paramters of DFTB calculations
* for electronic part:
 - confinement radii. After changing these the confined pseudo-atoms and the Slater-Koster
   files have to be regenerated.
 - Hubbard U
"""
import AtomicData

############ XC Functional for atomic orbitals ##########
pseudo_orbital_x = "gga_x_pbe" 
pseudo_orbital_c = "gga_c_pbe"

############ ELECTRONIC PART ############################
# radius for confinement potential Vconf(r) = (r/r0)^2
# is chosen to be r0=r0_factor*r_cov where r_cov is the covalent radius
#r0_factor = 1.85 # recommended by hotbit
r0_factor = 3.0 
confinement_radii_byZ = dict([(AtomicData.atomic_number(atom), r0_factor*rcov) \
                                for (atom, rcov) in AtomicData.covalent_radii.iteritems()])

# Hubbard U from experimental ionization energies and electron affinities
# U = IE-EA (in hartree)
hubbard_U_byZ = {}
for atom in AtomicData.ionization_energies_eV.keys():
    IE = AtomicData.ionization_energies_eV[atom]
    if AtomicData.electron_affinities_eV.has_key(atom):
        EA = AtomicData.electron_affinities_eV[atom]
        hubbard_U_byZ[AtomicData.atomic_number(atom)] = (IE-EA)/AtomicData.hartree_to_eV 

hubbard_U_byZ[7] = 0.53 # from hotbit N.elm, could not find experimental EA for nitrogen
#0.42 # taken from figure 2 in "J Chem Theory Comput. Apr 10, 2012; 7(4): 931-948 by Gaus and Elstner

## for testing we use hotbit parameters
#confinement_radii_byZ[1] = 1.482 # H
#hubbard_U_byZ[1] = 0.420
##confinement_radii_byZ[6] = 2.011 # 1.234 # C
#hubbard_U_byZ[6] = 0.365
## confinement radius for transition metals Sc,Ti,Fe,Co and Ni taken from J. Chem. Theory Comput. 2007, 3, 1349-1367
#confinement_radii_byZ[26] = 4.36
#hubbard_U_byZ[26] = 0.2005 # reference publication has different Hubbard parameters for each orbital: Us (0.2005) Up (0.2005) and Ud (0.35522)
## electronic zinc parameters are taken from Elstner et al. J Comput Chem (2003), 24, 565-581
#confinement_radii_byZ[30] = 4.9 
hubbard_U_byZ[30] = 0.153852 # Zinc has no stable negative ion, so there is no experimental electron affinity. 
                             # For U_H I take the value for U_p from the article 
                             # "DFTB Parameters for the Periodic Table Part I - Electronic Structure" 
                             # (J. Chem. Theory Comput. 2013, 9, 4006-4017)

# According to the 'NIST Atomic Spectra Database - Ionization Energies Data' the ionization energy of Ru^(2+) is
#  IE(Ru^2+) = 16.76 eV
# and the electron affinity is
#  EA(Ru^2+) = IE(Ru^2+) - IE(Ru^+) = 9.40 eV
# This gives a Hubbard parameters of  U_H(Ru^2+) = 7.36 eV, whereas for neutral Ruthenium U_H(Ru) = 6.31 eV
hubbard_U_byZ[44] = 7.36 / AtomicData.hartree_to_eV

# UGLY: hubbard_U_byZ is a global variable
def get_hubbard_parameters(parameter_set):
    global hubbard_U_byZ
    if parameter_set == "mio":
        # Hubbard Parameters by Niehaus from Phys. Rev. B 63 085108
        # for H,C,N and O
        hubbard_U_byZ[1] = 11.425 / AtomicData.hartree_to_eV
        hubbard_U_byZ[6] = 9.921 / AtomicData.hartree_to_eV
        hubbard_U_byZ[7] = 11.725 / AtomicData.hartree_to_eV
        hubbard_U_byZ[8] = 13.481 / AtomicData.hartree_to_eV
        return hubbard_U_byZ
    else:
        # default Hubbard parameters from experimental IE and EA
        # as defined above
        return hubbard_U_byZ

############# REPULSIVE POTENTIAL #########################

if __name__ == "__main__":
    import string
    # print default parameters for all atoms
    print "       DFTB Parameters"
    print "       ==============="
    print "Atom | r0/bohr  | U/hartree (eV)"
    print "-----|----------|-----------------"
    for atom in AtomicData.atom_names:
        try:
            r0 = "%.3f" % confinement_radii_byZ[AtomicData.atomic_number(atom)]
        except KeyError:
            r0 = "---"
        try:
            Uhub = hubbard_U_byZ[AtomicData.atomic_number(atom)]
            U = "  %.3f (%2.3f)" % (Uhub, Uhub*AtomicData.hartree_to_eV)
        except KeyError:
            U = "---"
        print string.center(atom, 5) + "|" + string.center(r0, 10) + "|" + string.center(U,10)
