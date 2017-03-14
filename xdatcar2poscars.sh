#!/bin/bash
#
# xdatcar2poscars.sh
#
# Convert an XDATCAR file into a set of POSCAR files.
#
# NB: Assumes an orthogonal simulation cell!
#
# Eero Holmstrom 2013
#

# usage
if [ -z "$1" ]
then
    echo "Usage: $(basename $0) [XDATCAR file]"
    exit 1;
fi

#
# get number of atoms from XDATCAR header
#
natoms=$(awk 'BEGIN{n=0};{n++};(n==7){natoms=0; for(i=1; i<=NF; i++){natoms=natoms+$i}; print natoms; exit;}' $1)

# assign input variables
f=$1

# prepare directory for the 1.coords, 2.coords, etc. files
rm -rf coordsframes
mkdir coordsframes

# prepare directory for the POSCAR.1 POSCAR.2 etc. files
rm -rf poscars
mkdir poscars

echo ""
echo "Cutting file $f into separate coordinate files..."

#
# remove 8-line header to make creating frames easier
#

rm -f coords.temp
nlines=$(wc $f | awk '{print $1}')
tail -$(($nlines-8)) $f > coords.temp

#
# split coords.temp into separate frames 1.coords, 2.coords etc. into
# directory coordsframes/
#

awk -v natoms=$natoms 'BEGIN{iframe=1; n=0};

(n>=0){print $1, $2, $3 > "coordsframes/"iframe".coords"};

{n++;};

(n==natoms){n=-1; close("coordsframes/"iframe".coords"); iframe++; };' coords.temp

# clean up
rm -f coords.temp

echo "Done."

#
# Loop over the coords files. For each frame, turn it into a POSCAR
# format file.
#

# get total number of frames
nframes=$(ls coordsframes/ | wc | awk '{print $1}')

#
# loop over frames
#

echo ""
echo "Creating the POSCAR files..."

iframe=1
while [ $iframe -le $nframes ]
do

#
# create CONTCAR version of frame
#

# start with the header, use "direct" VASP format
head -7 $f > temp.contcar
echo "Direct" >> temp.contcar

# then append the coordinates
cat coordsframes/$iframe.coords >> temp.contcar

# save file as a POSCAR
mv temp.contcar poscars/POSCAR.$iframe

# move to next frame
iframe=$(($iframe+1))

done

echo "Done."
echo ""

rm -rf coordsframes

echo "Exiting."
echo ""

exit 0

