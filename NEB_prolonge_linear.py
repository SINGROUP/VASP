#!/bin/python3
from ase import *
from ase.io import *
import numpy as np
import os
import time

num_im_orig = np.array([0 ,2 ,4 ,6 ,8]) # has to have initial and final image as well #
num_im_new  = np.array([  1 ,3 ,5 ,7  ]) # now inly one middle image in between        #

see = True


if see:
    from ase.visualize import view


for i,ii in enumerate(num_im_new):
	print ("image", i+1, "of", len(num_im_new));
	image_in = read('%02d/POSCAR' %num_im_orig[i] ); image_out = read('%02d/POSCAR' %num_im_orig[i+1] ); 
	if see:
		view(image_in); time.sleep(1)
	image = image_in.copy()
	pos_in = image_in.get_positions(); pos_out = image_out.get_positions();
	pos = (pos_in + pos_out)/(2.)
	image.set_positions(pos)
	if see:
		view(image); time.sleep(1)
	write('%02d/POSCAR' %(ii),image)
	print ("written: %02d/POSCAR" %(ii))
	del image_in;

if see:
	view(image_out)

print ("All done, now you can run NEB calc. by your-computer-self")


