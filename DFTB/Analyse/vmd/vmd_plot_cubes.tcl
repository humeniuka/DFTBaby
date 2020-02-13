#
# run this script in VMD as
#
#  vmd -e vmd_plot_cubes.tcl -args volumetric.cube  isovalue
#
# to plot the geometry and volumetric data in a cube file
#
# Arguments:
#   - `volumetric.cube` should in the Gaussian cube file format
#   - isovalue for plotting isosurface (floating point number)

axes location off
display projection orthographic
display depthcue off
color Display Background white

# input files
set cube_file [lindex $argv 0]
set isovalue  [lindex $argv 1]

# load geometry and volumetric data
mol new $cube_file

# representation for molecule
mol delrep 0 top
mol representation CPK 0.2 0.1 12 12
mol addrep top

# positive isosurface
# Representation 1: Positive isosurface
#mol material Transparent
mol material Opaque
mol color ColorID 0
mol representation Isosurface $isovalue 0.0 0.0 0.0
mol selection {all}
mol addrep top

# Representation 2: Negative isosurface
mol color ColorID 1
mol representation Isosurface -$isovalue 0.0 0.0 0.0
mol selection {all}
mol addrep top

## render image and save it to a file
#render TachyonInternal vectors.tga

#exit

