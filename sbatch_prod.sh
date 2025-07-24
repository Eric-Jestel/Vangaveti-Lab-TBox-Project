#!/bin/bash

# Note change simpart for each restart.
# Simpart = 0 is the starting point
# maxtime is in minutes

# sbatch uses "minutes" format 
# gromacs -maxh option (argument 4 to slurm) uses hours

#if [ `hostname` == "tesla1.rit.albany.edu" ]
#    then
#    maxtime=2880 # 2 days  
#    maxhour=`echo $maxtime | awk '{print $1/60}'`
#fi

##############################################
# Time limits
maxhour=200 #2 days
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
#grodir=/network/rit/lab/ChenRNALab/bin/gromacs-2016.4.hn7/bin/bin
Input_Var=$1
for i in ${Input_Var} #CYT_A CYT_C CYT_G CYT_U 
do
    curr=$dir/$i
    cd $curr
	    
    sbatch -p minerva -t 7-0 -n 1 -c 256 --mem=0 --exclude=minerva-[01-04] --job-name=$i-prod $dir/slurm_prod.sh $curr $nproc $grodir $maxhour
   
done




