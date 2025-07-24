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

$grodir/gmx grompp -f nvpt.mdp -c em.gro -r em.gro -p topol.top -o nvpt.tpr -n index.ndx -maxwarn 18

