
NEB_*.py are python3 (also working with python2) and ASE employing scripts for working with geometries for Nudget Elastic Band calculations in VASP.

NEB_contcar2poscar.py takes calculated CONTCAR from previous CALCULATIONS, checks whether some atoms didn't cross cell borders and adjust their position
according to the positition of initial image (00/POSCAR)

NEB_make_linear.py makes a linear path between POSCAR-in and POSCAR-out for wanted number of images. !!! It does not check for cell border crossing !!! 

# Feel free to add description for older *.sh scripts.
