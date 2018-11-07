#!/bin/python3
from ase import *
from ase.io import *
import numpy as np
import os

num_im = 3

see = True


if see:
    from ase.visualize import view

reactant = read('POSCAR-in')
product  = read('POSCAR-out')


pos_in  = reactant.get_positions()
pos_out = product.get_positions()

os.mkdir('00')
write('00/POSCAR', reactant)

if see:
	view(reactant)

os.mkdir('%02d' %(num_im+1))
write('%02d/POSCAR' %(num_im+1), product)

for i in range(num_im):
	print ("image", i+1, "of", num_im)
	os.mkdir('%02d' %(i+1))
	image = reactant.copy()
	pos = pos_in + (i+1)*(pos_out-pos_in)/(num_im+1)
	image.set_positions(pos)
	if see:
		view(image)
	write('%02d/POSCAR' %(i+1),image)
	print ("written: %02d/POSCAR" %(i+1))
	del image, pos;

if see:
	view(product)

print ("All done, now you can run NEB calc. by your-computer-self")


