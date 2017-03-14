#
# A script for extracting the total and projected DOS for a
# spin-polarized VASP run with LORBIT = 11. Takes as input the OUTCAR
# and DOSCAR as well as the first (a) and last (b) ion number to do
# the projected DOS over. Outputs the total DOS as well as the (s-,
# p-, and d-) projected DOS, the latter averaged over ions from i to
# j. Energy is shifted so that the Fermi energy is at E = 0 eV.
#
# Eero Holmstrom (2012, 2014, 2015)
#

# usage
if [ -z "$4" ]
then
    echo "Usage: $(basename $0) [OUTCAR] [DOSCAR] [a] [b] [multiply spin down dos by -1]"
    exit 1;
fi

# assign input arguments
outcar=$1
doscar=$2
a=$3
b=$4

if [[ -n "$5" ]];
then
    echo ""
    echo "NB: Multiplying spin down DOS by -1."
    multiplier=-1
else
    multiplier=1
fi

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

echo "Extracting total DOS, spin up..."
awk -v NEDOS=$NEDOS -v emax=$emax '{n++}(n>6 && NF==5 && $1<=emax && $3!=NEDOS){print $1, $2}' $doscar > gvd_spin_up_total_dos
echo "Done."
echo ""

echo "Extracting total DOS, spin up..."
awk -v multiplier=$multiplier -v NEDOS=$NEDOS -v emax=$emax '{n++}(n>6 && NF==5 && $1<=emax && $3!=NEDOS){print $1, multiplier*$3}' $doscar > gvd_spin_down_total_dos
echo "Done."
echo ""

#
# 2. Get projected DOS averaged over ions from ion number a to b
#

echo "Getting projected DOS for ions $a to $b..."

awk -v multiplier=$multiplier -v a=$a -v b=$b -v NEDOS=$NEDOS 'BEGIN{i=0; e=0};

{n++};

(n>6 && NF==5 && $3==NEDOS){i++; if(i>=a && i<=b){e=0; print "Processing ion number", i"...";} else {next;}};

(n>6 && NF==19 && i>=a && i<=b){e++; energy[e]=$1;

#
# s-type
#
total_spinup_s[e] = total_spinup_s[e] + $2;
total_spindown_s[e] = total_spindown_s[e] + multiplier*$3;

#
# p-type
#
total_spinup_p[e] = total_spinup_p[e] + $4 + $6 + $8;
total_spindown_p[e] = total_spindown_p[e] + multiplier*($5 + $7 + $9);

#
# d-type
#
total_spinup_d[e] = total_spinup_d[e] + $10 + $12 + $14 + $16 + $18;
total_spindown_d[e] = total_spindown_d[e] + multiplier*($11 + $13 + $15 + $17 + $19);

maxe=e;};

END{

print "All ions done."
print ""; print "Printing the data out to files...";

for(e=1; e <= maxe; e++){

#
# s-type
#
print energy[e], total_spinup_s[e] > "gvd_spin_up_total_s_"a"_to_"b"_dos";
print energy[e], total_spindown_s[e] > "gvd_spin_down_total_s_"a"_to_"b"_dos";

#
# p-type
#
print energy[e], total_spinup_p[e] > "gvd_spin_up_total_p_"a"_to_"b"_dos";
print energy[e], total_spindown_p[e] > "gvd_spin_down_total_p_"a"_to_"b"_dos";

#
# d-type
#
print energy[e], total_spinup_d[e] > "gvd_spin_up_total_d_"a"_to_"b"_dos";
print energy[e], total_spindown_d[e] > "gvd_spin_down_total_d_"a"_to_"b"_dos";

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
