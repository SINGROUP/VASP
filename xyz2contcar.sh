#!/bin/bash
#
# Convert an xyz file to VASP CONTCAR direct format.
#

# Usage
if [ ! -n "$2" ]; then
    echo "Usage: `basename $0` [XYZ file] [CONTCAR file]"
    exit
fi

# the input and output file names given as command-line parameters
xyz=$1
contcar=$2

rm -f $contcar

#
# Get the scaling factor and the lattice vectors. We'll be using the
# direct format.
#

scalfact=1.0
xsize=$(awk '{n++}(n==2){print $(NF-2)}' $xyz)
ysize=$(awk '{n++}(n==2){print $(NF-1)}' $xyz)
zsize=$(awk '{n++}(n==2){print $NF}' $xyz)

# get the number of atoms
natoms=$(awk '{print $1; exit}' $xyz)

#
# sort the atomic coordinates according to chemical species
#

tail -$natoms $xyz > $xyz.coords
sort $xyz.coords > $xyz.coords.sorted

#
# get the atomic species and the number of each species, print the
# results out to file types.dat
#

awk 'BEGIN{type=NA; ntypes=0}{

if ($1!=type) {type=$1; ntypes++; ntypearray[ntypes]=0; typearray[ntypes]=type};


}

($1==type){

ntypearray[ntypes] = ntypearray[ntypes] + 1;

}

END{

typestring="";
numberstring="";

for (i = 1; i <= ntypes; i++) {

typestring = typestring "  " typearray[i];
numberstring = numberstring "  " ntypearray[i];

}

print typestring;
print numberstring;

}' $xyz.coords.sorted > types.dat

#
# create the CONTCAR header
#

rm -f $contcar
echo "Created with xyz2contcar.sh" > $contcar
echo $scalfact >> $contcar
echo $xsize 0.0 0.0 >> $contcar
echo 0.0 $ysize 0.0 >> $contcar
echo 0.0 0.0 $zsize >> $contcar
head -1 types.dat >> $contcar
tail -1 types.dat >> $contcar
echo "Direct" >> $contcar

#
# then append the coordinates in the direct format
#

awk -v xsize=$xsize -v ysize=$ysize -v zsize=$zsize '{

print $2/xsize, $3/ysize, $4/zsize;

}' $xyz.coords.sorted >> $contcar

# clean up
rm -f $xyz.* types.dat

exit 0;
