#!/bin/bash

module purge
module load openfoam-2112-gcc-11.2.0-lhrpyq4

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions


# refined mesh case
caseName="cavityFine"
cd $caseName
runApplication reconstructPar -latestTime

paraFoam