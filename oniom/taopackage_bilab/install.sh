#!/bin/bash
#
# script to intall Toolkit to Assist ONIOM (TAO) calculations
#
# Peng Tao 2009 tao.21@osu.edu
#


# ATTENTION USERS:
# Please replace path ~/bin with the folder you want the links to
# TAO scripts to be installed. Please make sure this path is in your
# search PATH

mkdir bin
USERPATH=~/apps/taopackage/bin


DIR=$(pwd)

# Add full path 
cd $DIR/ESPT


# ATTENTION USERS:
# Please do NOT change anything beyond this point

perl -i -pe 's:REPLACETHISPATHPLEASE:'"$DIR"':g' AMBERFF.pm AMBERPARM.pm

cd $DIR/scripts

perl -i -pe 's:REPLACETHISPATHPLEASE:'"$DIR"':g' *.pl

cd $DIR


# Add link to each script to assigned path
for file in $DIR/scripts/*.pl; do

 filename=${file##*/}
 base=${filename%%.*}

 ln -s $file $USERPATH/$base

done

