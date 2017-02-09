#Centers the protein at the origin and calculates Calpha displacement from the center starts at the startInd residues and ending at 
#the endInd residue

#assumes protein oriented along z-axis

#Example run script:
#from Calphascript_all import *
#import MDAnalysis as md
#u = md.Universe("your_favorite_protein.gro","your_favorite_protein_trajectory.xtc")
#computeDistCenterAxis(u, "CA_dist.dat",1,100)



import numpy as np
import MDAnalysis as md
import math

def computeDistCenterAxis(u,filename, startInd,endInd):	
	all=u.select_atoms("protein")



	allInds=[]

	f=open(filename,'w')
	#header
	f.write("Residue_name "+ "Residue_number "+ "CA_dist_to_center\n")


	#atom selection based on index
	for i in range(startInd,endInd+1):
		choice = "protein and resid " + str(i) + " and name CA"
		tmpInd = u.select_atoms(choice)
		allInds.append(tmpInd)

	time=0

	#arrays for calculation across time and averages
	alpha_matrix = np.zeros((len(allInds), 2))
	time_matrix = np.zeros((len(u.trajectory),len(allInds)))


	#loop through all frames
	for ts in u.trajectory:
		#center protein about the origin
		all.positions=all.positions-all.center_of_geometry()
	
		
		#get cordinates for each atom and calculate distance to center
		#assumes protein is oriented along z-axis
		for i in range(0,len(allInds)):
			pos=allInds[i].positions
			x=pos[0,0]
			y=pos[0,1]
			time_matrix[time,i]=math.sqrt(x**2+y**2)
		time=time+1


	for i in range(0,len(allInds)):
		#get resid number
		alpha_matrix[i,0]=allInds[i].atoms.residues.resids[0]
		#get the time average of distances to the central axis
		alpha_matrix[i,1]=np.mean(time_matrix,axis=0)[i]
		f.write(allInds[i].atoms.residues.resnames[0] + ' ' +str(alpha_matrix[i,0]) + ' ' + str(alpha_matrix[i,1]) +'\n')

	f.close()
	return


