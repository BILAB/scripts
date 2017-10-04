#!/bin/sh

source /lustre/gi58/i58005/apps/gmx2016_3/bin/GMXRC.bash

gro_file=$1
cycle=$2
threads=$3

gmx grompp \
    -f runset.mdp \
    -c $gro_file \
    -p ../../grotop/LysYZW.top \
    -n ../../grotop/index.ndx \
    -o cyc${cycle}/${cycle}-${threads}.tpr \
    -po cyc${cycle}/${cycle}-${threads}.out.mdp \
   -maxwarn 1
