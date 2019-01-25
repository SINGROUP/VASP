#!/bin/python3
from ase import *
from ase.io import *
import numpy as np
import os
from shutil import copyfile
from time import sleep

# === options === #

see = True # visualize images?#
aimsGeom = True # save files also in 'geometry.in' files of FHI-AIMS? They are easy to visualize via JMOL .#
delay = True # To delay visualization in between images so they would be organized in proper order? 

num_im = 3 # Number of images # 
save_dist = 4.5 # Save distance of atoms between images #

# === script body === #
if see:
	from ase.visualize import view

reactant = read('00/POSCAR')  # BEWARE: sometimes it is necessary to write elements into the POSCAR-files (note done automatically in ASE)#
product  = read('%02d/POSCAR' %(num_im+1))

pos_in  = reactant.get_positions()
no_at   = len(pos_in)
cell    = reactant.get_cell()
#pos_out = product.get_positions()

if see:
	view(reactant)
if delay:
	sleep(1)


for i in range(num_im):
	print ("image", i+1, "of", num_im)
	image = read('%02d/CONTCAR' %(i+1))
	copyfile('%02d/CONTCAR' %(i+1), '%02d/CONTCAR-old' %(i+1))
	pos   = image.get_positions()
	for ia in range(no_at):
		if (np.linalg.norm(pos_in[ia]-pos[ia]) >= save_dist):
			print ("atom", ia, "has distance from previous image longer than", save_dist)
			tmp = np.array([ pos[ia],            (-1*cell[0])+(-1*cell[1])+pos[ia], (0*cell[0])+(-1*cell[1])+pos[ia],
				(1*cell[0])+(-1*cell[1])+pos[ia], (1*cell[0])+(0*cell[1])+pos[ia], (1*cell[0])+(1*cell[1])+pos[ia],
				(0*cell[0])+(1*cell[1])+pos[ia], (-1*cell[0])+(1*cell[1])+pos[ia], (-1*cell[0])+(0*cell[1])+pos[ia] ])
			#print "DEBUG: tmp",tmp
			tmp_dist = np.zeros(9);
			for j in range(9):
				tmp_dist[j] = np.linalg.norm(pos_in[ia]-tmp[j])
			#print "DEBUG: tmp_dist",tmp_dist
			ja = np.argmin(tmp_dist)
			print ("the new distance will be:", tmp_dist[ja])
			#print "DEBUG: tmp[ja]", tmp[ja]
			pos[ia]=tmp[ja]
	image.set_positions(pos)
	if see:
		view(image)
	if delay:
		sleep(1)
	write('%02d/POSCAR' %(i+1),image)
	print ("written: %02d/POSCAR" %(i+1))
	if aimsGeom:
		write('%02d/geometry.in' %(i+1),image)
		print ("written: %02d/geometry.in" %(i+1))
	pos_in = pos

if see:
	view(product)

print ("All done, go wherever you want")


