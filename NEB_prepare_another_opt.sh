#!/bin/bash

str='01 02 03' # adjust according to name of directories of inner images #
j=1

for i in $str;
do
    echo $i
    cd $i
    cp CONTCAR CONTCAR-$j
    mv OUTCAR OUTCAR-$j
    cd ..
done

# uncomment following line if running on super computers #
#module load ase 

# both python & python3 can be used here #
python NEB_contcar2poscar.py