#!/bin/bash
#
# xdatcar2xyz_movie.sh
#
# Convert a VASP XDATCAR file to an xyz movie. Uses contcar2xyz.sh.
#
# NB: Assumes an orthogonal simulation cell!
#
# Eero Holmstrom 2013
#

# usage
if [ -z "$3" ]
then
    echo "Usage: $(basename $0) [input XDATCAR] [input OUTCAR] [output XYZ file]"
    exit 1;
fi

# assing input arguments
infile=$1
outcar=$2
outfile=$3

#
# get the time step POTIM from OUTCAR
#
potim=$(awk '/POTIM/{print $3; exit}' $outcar)

echo ""
echo "Found time step of $potim fs."

#
# get number of atoms from XDATCAR header
#
natoms=$(awk 'BEGIN{n=0};{n++};(n==7){natoms=0; for(i=1; i<=NF; i++){natoms=natoms+$i}; print natoms; exit;}' $infile)

echo "Found $natoms ions."

#
# create a separate file for each frame with just the raw coordinates
#

# prepare directory for the 1.coords, 2.coords, etc. files.
rm -rf coordsframes
mkdir coordsframes

#
# remove 8-line header to make creating frames easier
#

rm -f coords.temp
nlines=$(wc $infile | awk '{print $1}')
tail -$(($nlines-8)) $infile > coords.temp

#
# split coords.temp into separate frames 1.coords, 2.coords etc. into
# directory coordsframes/
#

echo ""
echo "Splitting $infile into coordinate files..."

awk -v natoms=$natoms 'BEGIN{iframe=1; n=0};

(n>=0){print $1, $2, $3 > "coordsframes/"iframe".coords"};

{n++;};

(n==natoms){n=-1; close("coordsframes/"iframe".coords"); iframe++; };' coords.temp

echo "Done."

# clean up
rm -f coords.temp

#
# Loop over the coords files. For each frame, turn it into a CONTCAR
# format file, convert it into XYZ, and append it to the movie file.
#

# get total number of frames
nframes=$(ls coordsframes/ | wc | awk '{print $1}')

#
# loop over frames
#

rm -f $outfile
iframe=1

echo ""
echo "Looping over the coordinate files, turning them into XYZ and combining them..."

while [ $iframe -le $nframes ]
do

#
# create CONTCAR version of frame
#

# start with the header, use "direct" VASP format
head -7 $infile > temp.contcar
echo "Direct" >> temp.contcar

# then append the coordinates
cat coordsframes/$iframe.coords >> temp.contcar

#
# convert CONTCAR frame into xyz
#
contcar2xyz.sh temp.contcar temp.xyz

#
# add physical time and frame counter to xyz header
#
awk -v potim=$potim -v iframe=$iframe '{

if($1=="Frame"){print "Frame number", iframe, iframe*potim, "fs boxsize", $(NF-2), $(NF-1), $(NF)} else {print $0};

}' temp.xyz > temp.xyz.2

# append xyz file to movie file
cat temp.xyz.2 >> $outfile

# clean up
rm -f temp.contcar temp.xyz temp.xyz.2

# move to next frame
iframe=$(($iframe+1))

done

echo "Done."
echo ""

# clean up
rm -rf coordsframes

echo "Exiting."
echo ""

exit 0;

