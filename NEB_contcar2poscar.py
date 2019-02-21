#!/bin/python3
from ase import *
from ase.io import *
import numpy as np
import os
from shutil import copyfile
from time import sleep

# === options === #

see = True # visualize images? #
aimsGeom = True # save files also in 'geometry.in' files of FHI-AIMS? They are easy to visualize via JMOL .#
delay = True # To delay visualization in between images so they would be organized in proper order? #

num_im = 3 # Number of images # 
save_dist = 4.5 # Save distance of atoms between images #

safety_check = True # Necessary only, when there are constrained atoms possibly moving through the cell border  #
                    # slows down everything, but !!! sometimes constrained atoms are automatically folded into the original cell!!! This procedure should overcome it. #

# === script body === #
if see:
	from ase.visualize import view

reactant = read('00/POSCAR')  # BEWARE: sometimes it is necessary to write elements into the POSCAR-files (note done automatically in ASE) #
product  = read('%02d/POSCAR' %(num_im+1))

pos_in  = reactant.get_positions()
no_at   = len(pos_in)
cell    = reactant.get_cell()
posd_in = np.array([np.linalg.solve(cell, ia) for ia in pos_in]) # unfolded DIRECT coordinates #

if see:
	view(reactant)
if delay:
	sleep(1)

#if aimsGeom: # uncomment this, if you want to write down geometry.in files for reactant #
#	write('00/geometry.in', reactant)
#	print ("written: 00/geometry.in" )

for i in range(num_im):
	print ("image", i+1, "of", num_im)
	image = read('%02d/CONTCAR' %(i+1))
	copyfile('%02d/CONTCAR' %(i+1), '%02d/CONTCAR-old' %(i+1))
	pos   = image.get_positions()
	posd  = np.array([np.linalg.solve(cell, ia) for ia in pos]) # unfolded DIRECT coordinates #
	tmpc  = np.zeros(pos.shape)
	for ia, posd_ia in enumerate(posd[:]):
		if (np.linalg.norm(pos_in[ia]-pos[ia]) >= save_dist):
			print ("atom", ia, "has distance from previous image longer than", save_dist)
			print ("the input    DIRECT coordinates are:", posd_in[ia] )
			print ("the original DIRECT coordinates are:", posd_ia )
			for j in range(3):
				hold = [abs(posd_ia[j] + b - posd_in[ia][j]) for b in [-1, 0, 1]]
				tmpc[ia,j] = [-1, 0, 1][np.argmin(hold)]
				posd_ia[j] += tmpc[ia,j]
			print ("the new DIRECT coordinates will be :", posd_ia )
			pos[ia]= np.dot(cell,posd_ia); posd[ia]=posd_ia;
			print ("the new cartesian coord:s will be  :", pos[ia] )
	image.set_positions(pos)
	if safety_check: # procedure to check if ASE does not fold constrain atoms to the original cell automatically #
		move = np.zeros(pos.shape)
		tmp_pos = image.get_positions(); folded = False;
		for ia, tmp_ia in enumerate(tmp_pos): 
			if (np.linalg.norm(tmp_ia-pos[ia]) >= 0.1):
				folded = True; print ("Wrongly folded atom:", ia);
				move[ia] = np.dot(cell,tmpc[ia])
		if folded:
			print ("moving wrongly folded atoms")
			image.translate(move)
	del tmpc; 
	# follows visualizations & writing procedures #
	if see:
		view(image)
	if delay:
		sleep(1)
	write('%02d/POSCAR' %(i+1),image)
	print ("written: %02d/POSCAR" %(i+1))
	if aimsGeom:
		write('%02d/geometry.in' %(i+1),image)
		print ("written: %02d/geometry.in" %(i+1))
	pos_in = pos; posd_in = posd; 

if see:
	view(product)
#if aimsGeom: # uncomment this, if you want to write down geometry.in files for product #
#	write('%02d/geometry.in' %(num_im+1),product)
#	print ("written: %02d/geometry.in" %(num_im+1))

print ()
print ("All done, go wherever you want")

