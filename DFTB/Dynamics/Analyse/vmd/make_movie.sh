#!/bin/bash
# render frames
vmd -e pyrene_crystal_trajectories.tcl
# make animated gif
echo "converting frames to animated gif...PLEASE WAIT..."
convert -delay 0.01 snap.0*.rgb movie.gif
# remove temporary files
rm -f snap.*.rgb
