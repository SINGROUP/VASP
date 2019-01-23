#!/bin/python3
from ase import *
from ase.io import *
import numpy as np
import os
from shutil import copyfile

from ase.visualize import view

num_im=3

reactant = read('00/POSCAR')  # BEWARE: sometimes it is necessary to write elements into the POSCAR-files (note done automatically in ASE)#
product  = read('%02d/POSCAR' %(num_im+1))

view(reactant)

for i in range(num_im):
	print ("image", i+1, "of", num_im)
	image = read('%02d/POSCAR' %(i+1))
	view(image)

view(product)

print ("All done, go wherever you want")


