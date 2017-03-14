#!/bin/bash
#
# Convert a VASP POSCAR file into the xyz format.
#

# Usage
if [ ! -n "$2" ]; then
    echo "Usage: `basename $0` [CONTCAR file] [xyz file]"
    exit
fi

# the input and output file names given as command-line parameters
contcar=$1
xyz=$2

rm -f $xyz

# get the scaling factor and the lattice vectors
scalfact=$(awk '{n=n+1; if (n==2){print $1; exit} }' $contcar)
xsize=$(awk '{n=n+1; if(n==3){print $1; exit} }' $contcar)
ysize=$(awk '{n=n+1; if(n==4){print $2; exit} }' $contcar)
zsize=$(awk '{n=n+1; if(n==5){print $3; exit} }' $contcar)

# multiply the lattice vectors by the scaling factor
xsize=$( echo "$xsize*$scalfact" | bc -l)
ysize=$( echo "$ysize*$scalfact" | bc -l)
zsize=$( echo "$zsize*$scalfact" | bc -l)

# get the number of atoms
natoms=$(awk 'BEGIN{sum=0}{n=n+1; if(n==7){ for(i=1; i<=NF; i++){sum=sum+$i}}}END{print sum}' $contcar)

# are we direct or cartesian
if [[ -n "$(grep -i selective $contcar)" ]];
    then
    format=$(awk '{n=n+1; if(n==9){print $1} }' $contcar)
    selectiveDynamics=1
else
    format=$(awk '{n=n+1; if(n==8){print $1} }' $contcar)
    selectiveDynamics=0
fi

if [ "$format" == "Direct" ]; then
    isdirect=1
else
    isdirect=0
fi

# create the xyz file

# number of atoms
echo $natoms > $xyz
# header with cell size
printf "Frame number 1 1 fs boxsize %8.3f %8.3f %8.3f\n" $xsize $ysize $zsize >> $xyz

#
# selective dynamics or not
# 
if [[ selectiveDynamics -eq 0 ]];
then
    
    #
    # coordinates, from either direct or cartesian format
    #
    
    if [ $isdirect -eq 1 ]; then
	
	awk -v selectiveDynamics=$selectiveDynamics -v xsize=$xsize -v ysize=$ysize -v zsize=$zsize 'BEGIN{j=1; k=1; n=0; natoms=0}{

n=n+1;

if(n==6){for(i=1; i<=NF; i++){types[i]=$i}};
if(n==7){for(i=1; i<=NF; i++){ntypes[i]=$i; natoms=natoms+$i}; numberoftypes=i};

# create a simple natoms long array typelist with the types in the right order

# loop over types
for(i=1; i<=numberoftypes; i++) {

# loop over number of atoms of this type
for(l=1; l<=ntypes[i]; l++) {
 typelist[j] = types[i];
 j++;
}

}

if(n>8 && n<=8+natoms){ print typelist[k], $1*xsize, $2*ysize, $3*zsize, "1"; k++};

}' $contcar >> $xyz
	
    else
	
	awk -v scalfact=$scalfact -v xsize=$xsize -v ysize=$ysize -v zsize=$zsize 'BEGIN{j=1; k=1; n=0; natoms=0}{

n=n+1;

if(n==6){for(i=1; i<=NF; i++){types[i]=$i}};
if(n==7){for(i=1; i<=NF; i++){ntypes[i]=$i; natoms=natoms+$i}; numberoftypes=i};

# create a simple natoms long array typelist with the types in the right order

# loop over types
for(i=1; i<=numberoftypes; i++) {

# loop over number of atoms of this type
for(l=1; l<=ntypes[i]; l++) {
 typelist[j] = types[i];
 j++;
}

}

if(n>8 && n<=8+natoms){ print typelist[k], $1*scalfact, $2*scalfact, $3*scalfact, "1"; k++};

}' $contcar >> $xyz
	
    fi
    
else
    
    #
    # coordinates, from either direct or cartesian format
    #
    
    if [ $isdirect -eq 1 ]; then
	
	awk -v selectiveDynamics=$selectiveDynamics -v xsize=$xsize -v ysize=$ysize -v zsize=$zsize 'BEGIN{j=1; k=1; n=0; natoms=0}{

n=n+1;

if(n==6){for(i=1; i<=NF; i++){types[i]=$i}};
if(n==7){for(i=1; i<=NF; i++){ntypes[i]=$i; natoms=natoms+$i}; numberoftypes=i};

# create a simple natoms long array typelist with the types in the right order

# loop over types
for(i=1; i<=numberoftypes; i++) {

# loop over number of atoms of this type
for(l=1; l<=ntypes[i]; l++) {
 typelist[j] = types[i];
 j++;
}

}

if(n>9 && n<=9+natoms){ print typelist[k], $1*xsize, $2*ysize, $3*zsize, "1"; k++};

}' $contcar >> $xyz
	
    else
	
	awk -v scalfact=$scalfact -v xsize=$xsize -v ysize=$ysize -v zsize=$zsize 'BEGIN{j=1; k=1; n=0; natoms=0}{

n=n+1;

if(n==6){for(i=1; i<=NF; i++){types[i]=$i}};
if(n==7){for(i=1; i<=NF; i++){ntypes[i]=$i; natoms=natoms+$i}; numberoftypes=i};

# create a simple natoms long array typelist with the types in the right order

# loop over types
for(i=1; i<=numberoftypes; i++) {

# loop over number of atoms of this type
for(l=1; l<=ntypes[i]; l++) {
 typelist[j] = types[i];
 j++;
}

}

if(n>9 && n<=9+natoms){ print typelist[k], $1*scalfact, $2*scalfact, $3*scalfact, "1"; k++};

}' $contcar >> $xyz
	
    fi
fi

exit 0;

