#!/bin/bash

str=`ls -d 0*` # possible automatic finding of NEB images (if less than 9) #
#str='00 01 02 03 04' # optional hard definition #

mv energy_tmp.txt energy_tmp-old.txt

for i in $str;
do
    echo $i
    j=`grep "energy  without entropy" $i/OUTCAR | tail -1`
    echo $i " " $j >> energy_tmp.txt
done

echo "done (energy greped) " 