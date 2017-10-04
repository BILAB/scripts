#! /bin/sh

source /lustre/gi58/i58005/apps/gmx2016_3/bin/GMXRC.bash

cycle=$1
threads=$2

gmx mdrun -v \
    -s cyc${cycle}/${cycle}-${threads}.tpr \
    -deffnm cyc${cycle}/${cycle}-${threads}
