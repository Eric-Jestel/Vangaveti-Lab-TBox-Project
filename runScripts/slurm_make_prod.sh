#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Error in $0 - Invalid Argument Count"
    echo "Usage : $0 dir nproc gromacs-bindir"
    exit
fi

module load mpi/openmpi-x86_64


cd $1
nproc=$2
grodir=$3
maxtime=$4

$grodir/gmx grompp -f md.mdp -c nvpt.gro -t nvpt.trr -p topol.top -o md_0_1.tpr -n index.ndx -maxwarn 18

