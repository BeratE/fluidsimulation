#!/bin/sh

# Usage ./reconstruct.sh (surface program) (directory with vtk files)
# Execute the surface programm for each file in the given directy

find . "$(pwd)/$2" -maxdepth 1 -regex ".*\.\(vtk\)" -exec $1 {} \;
