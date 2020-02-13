
The files "pyrene_crystal.xyz" and "pyrene_crystal.ff" contain the same geometry. The bonding topology
for constructing the periodic force field is deduced from the atom positions given in "pyrene_crystal.ff",
while positions in "pyrene_crystal.xyz" are used to calculate the energy and gradients.
This separation into two files allows you to evaluate the force field even for very distorted geometries
where a correct assignment of bonds based on the atomic positions would not be possible.

The file "pyrene_crystal.ff" is structured similarly to a normal XYZ-file, with the following differences:
 - The 5th column specifies the atom type (an integer).
 - At the end of the geometry section, the lattice vectors a,b,c are given on the last 3 lines (in Angstrom):
      Tv  ax ay az
      Tv  bx by bz
      Tv  cx cy cz

To find the geometry with minimum energy of the combined QM/MM system on the ground state, run

   optimize.py pyrene_crystal.xyz 0 --verbose=0

in the QMMM folder. The configuration for DFTB with the indeces of the QM atoms
will be read from dftbaby.cfg, and the force field definition from 'pyrene_crystal.ff'.

