#
# A script for extracting the total and projected DOS for a
# non-spin-polarized VASP run with LORBIT = 11. Takes as input the
# OUTCAR and DOSCAR as well as the first (a) and last (b) ion number
# to do the projected DOS over. Outputs the total DOS as well as the
# (s-, p-, and d-) projected DOS, the latter averaged over ions from i
# to j. Energy is shifted so that the Fermi energy is at E = 0 eV.
#
# Eero Holmstrom (2012, 2014, 2015)
#

# usage
if [ -z "$4" ]
then
    echo "Usage: $(basename $0) [OUTCAR] [DOSCAR] [a] [b]"
    exit 1;
fi

# assign input arguments
outcar=$1
doscar=$2
a=$3
b=$4

# find number of ions
nions=$(awk '/NIONS/{print $NF}' $outcar)

# find the maximum energy in the DOS
emax=$(awk '{n++}(n==6){print $1; exit}' $doscar)

#
# get the Fermi energy
#

ISMEAR=$(awk '/ISMEAR/{print substr($3,1,2); exit;}' $outcar)
EFermi=$(grep 'E-fermi' $outcar | tail -1 | awk '{print $3}')

# get NEDOS
NEDOS=$(awk '/number of dos/{print $6; exit;}' $outcar)

# print out stats to the user
echo ""
echo "Found $nions ion(s)."
echo "Found ISMEAR = $ISMEAR."
echo "Found Fermi energy of $EFermi eV".
echo "Found maximum energy of $emax eV".
echo "Found NEDOS of $NEDOS."
echo ""

#
# 1. Extract total DOS
#

echo "Extracting total DOS..."
awk -v NEDOS=$NEDOS -v emax=$emax '{n++}(n>6 && NF==3 && $1<=emax && $3!=NEDOS){print $1, $2}' $doscar > gvd_total_dos
echo "Done."
echo ""

#
# 2. Get projected DOS averaged over ions from ion number a to b
#

echo "Getting projected DOS for ions $a to $b..."

awk -v a=$a -v b=$b -v NEDOS=$NEDOS 'BEGIN{i=0; e=0};

{n++};

(n>6 && NF==5 && $3==NEDOS){i++; if(i>=a && i<=b){e=0; print "Processing ion number", i"...";} else {next;}};

(n>6 && NF==10 && i>=a && i<=b){e++; energy[e]=$1;

#
# s-type
#
total_s[e] = total_s[e] + $2;

#
# p-type
#
total_p[e] = total_p[e] + $3 + $4 + $5;

#
# d-type
#
total_d[e] = total_d[e] + $6 + $7 + $8 + $9 + $10;

maxe=e;};

END{

print "All ions done."
print ""; print "Printing the data out to files...";

for(e=1; e <= maxe; e++){

#
# s-type
#
print energy[e], total_s[e] > "gvd_total_s_"a"_to_"b"_dos";

#
# p-type
#
print energy[e], total_p[e] > "gvd_total_p_"a"_to_"b"_dos";

#
# d-type
#
print energy[e], total_d[e] > "gvd_total_d_"a"_to_"b"_dos";

};

print "Done."

}' $doscar

#
# rename the unshifted dos files
#

for f in gvd*dos
do

cp $f $f".dat"

done

#
# shift fermi energy to zero and rename the shifted files
#

echo ""
echo "Shifting Fermi energy to zero for all data sets..."

for f in gvd*dos
do

awk -v EFermi=$EFermi '{print $1-EFermi, $2}' $f > $f"_shifted.dat"
rm $f

done

echo "Done."

# done
echo ""
echo "Exiting."
echo ""

exit 0;
