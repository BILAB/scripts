#!/bin/bash

for file in *.pl; do

 filename=${file##*/}
 base=${filename%%.*}

 echo $filename
 echo $base
 pod2man $filename | groff -man -Tps -t > $base.ps
 ps2pdf $base.ps


done

