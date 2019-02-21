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
rotation = '0x,0z' # from which side povray files should be written: '0x,0z' from above ; '90x,180z' from x side #
povray_name = 'neb_tmp' # how the povray files should be named #

reactant = read('00/POSCAR')  # BEWARE: sometimes it is necessary to write elements into the POSCAR-files (note done automatically in ASE)#
product  = read('%02d/POSCAR' %(num_im+1))

view(reactant)
if wrpov:
	i=-1
	write(povray_name+'_%02d.pov' %(i+1),reactant, rotation=rotation)

time.sleep(1)

for i in range(num_im):
	print ("image", i+1, "of", num_im)
	image = read('%02d/POSCAR' %(i+1))
	if wrpov:
		write(povray_name+'_%02d.pov' %(i+1),image, rotation=rotation)
	view(image)
	time.sleep(1)

view(product)
if wrpov:
	i=num_im
	write(povray_name+'_%02d.pov' %(i+1),product, rotation=rotation)

print ("All done, go wherever you want")


