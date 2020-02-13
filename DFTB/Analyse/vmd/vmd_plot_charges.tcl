#
# run this script in VMD as
#
#  vmd -e vmd_plot_charges.tcl -args  geometry.xyz  charges.chg
#
# to show charges on each atom.
#
# Arguments:
#   - `geometry.xyz` should contain the molecular geometry in XYZ-format
#   - `charges.chg` should contain the same molecular geometry in XYZ-format
#     with an additional 5th column containing charges.
#

proc align_axes {sel} {
    #
    # aligns the principal axes of inertia of the molecule
    # with the global axes
    #
    # Returns
    # =======
    #   m4     :  the 4x4 transformation matrix
    
    # measure tensor of inertia (a 3x3 matrix)
    set m3 [lindex [measure inertia $sel] 1]
    # build 4x4 transformation matrix
    set v0 "[lindex $m3 0] 0"
    set v1 "[lindex $m3 1] 0"
    set v2 "[lindex $m3 2] 0"
    set v3 {0 0 0  1}
    set m4 [list $v0 $v1 $v2 $v3]
    puts "Transformation matrix"
    puts $m4
    # rotate molecule so that the principal axes of inertia
    # coincide with the x,y and z axes
    $sel move $m4

    return $m4
}


axes location off
display projection orthographic
display depthcue off
color Display Background white

# input file
set xyz_file [lindex $argv 0]
set chg_file [lindex $argv 1]

# load geometry
mol new $xyz_file

# choose representation
mol delrep 0 top
mol representation CPK 0.2 0.1 12 12
mol addrep top


set atoms [atomselect top "all"]

# rotate molecule, so that it is aligned with the global axes
set m4 [align_axes $atoms]
# top view on molecule
rotate y by -90

set charges_fh [open $chg_file]
# read two lines with number of atoms and comment
gets $charges_fh
gets $charges_fh

set i 0
foreach coord [$atoms get {x y z}] name [$atoms get {name}] {
    # read next line from file which is not a comment
    set line "#"
    while {[regexp "#" $line] == 1} {
	# skip comment
	set line [gets $charges_fh]
    }
    # remove multiple whitespaces
    set words [join $line " "]
    # split line into words and assign them to variables
    lassign [split $words] atname xx yy zz charge 

    puts "charge $atname-$i = $charge"
    # put label with charge
    if {$charge < 0} {
	# blue for negative charges
	graphics 0 color blue
    } else {
	# red for positive charges
	graphics 0 color red
    }
    set label [format "%+2.2f" $charge]
    lassign $coord x y z
    set x [expr { $x + 0.25} ]
    set z [expr { $z + 0.25} ]
    graphics 0 text "$x $y $z" $label size 2 thickness 4
    
    incr i
}

close $charges_fh

# render image and save it to a file
#render TachyonInternal vectors.tga
#
#exit

