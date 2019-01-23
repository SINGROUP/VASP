# Scripts for preparement & processing of VASP inputs/outputs/
_NEB_*.py are python3 (also working with python2) and ASE employing scripts for working with geometries for Nudget Elastic Band calculations in VASP._

__NEB_contcar2poscar.py__ takes calculated CONTCAR from previous CALCULATIONS, checks whether some atoms didn't cross cell borders and adjust their position
according to the positition of initial image (00/POSCAR)

__NEB_make_linear.py__ makes a linear path between POSCAR-in and POSCAR-out for wanted number of images. !!! It does not check for cell border crossing - you have to prepare POSCAR-in and POSCAR-out by yourself !!! 

__NEB_prolonge_linear.py__ put images in between already calculated images to prolonge the NEB path, to be more precise. The new images are calculated as linear combination of neighbouring images.

__NEB_view_poscar.py__ checking script to visualize currently prepared NEB path. There are 1 second delays between each visualization -- meaning that the visualized images should be ordered along the path.

### example of preparation of inputs for NEB:

Input geometry is stored in __POSCAR-pre-in__ from previous calculations:
> from ase.visualize import *; from ase import *; from ase.io import *; neco = read('POSCAR-pre-in'); view(neco); # to import libraries and see the input #

> import numpy as np ; move = np.zeros((neco.get_number_of_atoms(),3)); # create a field where translation of position of atoms is stored #

> cell = neco.get_cell();
> move[XX1]=move[XX2]=-cell[0]; # moving atoms to the other side of the cell #

> move[XX3]=move[XX4]=cell[0]

> move[XX5]=move[XX6]=cell[1]

> neco.translate(move) # application of atoms translation #

> view(neco); # to check the changes # 

> write('POSCAR-in',neco) # writing geometry in POSCAR-in#

Similar procedure has to be done for POSCAR-out. Don't forget to visually check everything while running __NEB_contcar2poscar.py__.

### Feel free to add description for older *.sh scripts.
