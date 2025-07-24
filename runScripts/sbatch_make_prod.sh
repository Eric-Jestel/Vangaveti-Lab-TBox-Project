#!/bin/bash

# Note change simpart for each restart.
# Simpart = 0 is the starting point
# maxtime is in minutes

# sbatch uses "minutes" format 
# gromacs -maxh option (argument 4 to slurm) uses hours

#if [ `hostname` == "tesla1.rit.albany.edu" ]
#    then
#    maxtime=24 # 2 days  
#    maxhour=`echo $maxtime | awk '{print $1/60}'`
#fi

##############################################
# Time limits
maxhour=24 #1 day
maxtime=`echo $maxhour | awk '{print $1*60}'`
##############################################

######################
# Processor allocation
nproc=256

simpart=0
######################
dir=`pwd`
grodir=/network/rit/lab/ChenRNALab/bin/gromacs-2025.2/bin/bin/
source ${grodir}/GMXRC.bash
Input_Var=$1
for i in ${1} #RI_RG RI_RA RI_RC RI_RU
do
    curr=$dir/$i
    #echo $curr
    sbatch -p minerva -t 7-0 -n 1 -c 256 --mem=0 --exclude=minerva-[01-04] --job-name=$i-equib\
	  slurm_make_prod.sh $curr $nproc $grodir $maxhour
done
