#!/bin/python3
from ase import *
from ase.io import *
import numpy as np
import os
from shutil import copyfile
import time
from ase.visualize import view

num_im= 3 # adjust according to number of images in NEB calc #
wrpov = False # write povray files #

reactant = read('00/POSCAR')  # BEWARE: sometimes it is necessary to write elements into the POSCAR-files (note done automatically in ASE)#
product  = read('%02d/POSCAR' %(num_im+1))

view(reactant)
time.sleep(1)

for i in range(num_im):
	print ("image", i+1, "of", num_im)
	image = read('%02d/POSCAR' %(i+1))
	if wrpov:
		write('neb_tmp_%02d.pov',image)
	view(image)
	time.sleep(1)

view(product)

print ("All done, go wherever you want")


