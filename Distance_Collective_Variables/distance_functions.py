
# USAGE:
# from distance_functions import *

# PREAMBLE:

import numpy as np
import MDAnalysis

sums = np.sum
square = np.square
zeros = np.zeros

# SUBROUTINES/FUNCTIONS:

def RMSD(x,y,n):
	""" Calculates the Root Mean Squared Distance between two arrays of the same size

	Usage: rmsd = RMSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape (n X 3)
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance between x and y arrays

	"""
	
	return (sums(square(x-y))/n)**0.5	# the MSD should never be negative, so using **0.5 rather than np.sqrt is safe

def MSD(x,y,n):
	""" Calculates the Mean Squared Distance between two arrays of the same size

	Usage: msd = MSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance squared between x and y arrays

	"""

	return sums(square(x-y))/n

def wrapping(x,dim):
	""" Calculates the translation matrix needed to wrap a particle back into the original periodic box
	
	Usage: t = wrapping(x,dim)

	Arguments:
	x: a numpy array of size (3) that corresponds to the xyz coordinates of an ATOM/COM/COG of a residue
	dim: a numpy array of size (3) that holds the xyz dimensions of the periodic box at that timestep

	"""
	
	t = zeros(3)
	dim2 = dim/2.
	for i in range(3):
		if (x[i]<-dim2[i]) or (x[i]>dim2[i]):
			t[i] = -dim[i]*round(x[i]/dim[i])
	return t

def euclid_dist(x,y):
	""" Calculates the Euclidian Distance between two arrays of the same size
	Usage: dist,dist2 = euclid_dist(x,y)
		
	Arguments:
	x, y: numpy arrays with the same size
	"""
	
	dist2 = sums(square(x-y))
	dist = dist2**0.5	# the MSD should never be negative, so using **0.5 rather than np.sqrt is safe
	return dist, dist2

def dist_matrix_calc(pdb,atom_selections,traj_loc,start,end,system_descriptor,ignore_n_nearest_neighbors=0,step=1):
        """
        """

        # ----------------------------------------
        # CREATE OUTPUT FILE NAMING VARIABLES
        node_output_filename = system_descriptor + '.nodes.txt'
        selection_output_filename = system_descriptor + '.selections.txt'
        data_output_filename = system_descriptor + '.col_var.dat'

        # ----------------------------------------
        # LOAD IN AND CREATE ATOM SELECTIONS IN THE ANALYSIS UNIVERSE OBJECT
        u = MDAnalysis.Universe(pdb)
        col_var_selections = u.select_atoms(atom_selections)  # MDAnalysis atom selection string formatting required.
        nNodes = col_var_selections.n_atoms     # assumes each col var is a distance between a pair of atoms, i and j. 
        nNodes_range = range(nNodes)
	i_max = nNodes-1-ignore_n_nearest_neighbors   # max value of atom looping index i; to be used multiple times so why calculate it multiple times
        boolean_matrix = np.full((nNodes,nNodes),False)  # 2D matrix of False values; elements of this matrix will be set to True as we loop over our i,j atom pairs. This will be used later for the plotting of 1D collective variable vectors onto the respective 2D atom pair matrix.

        count = 0
        with open(selection_output_filename,'w') as W, open(node_output_filename,'w') as Y:
                Y.write('# Node description: atom name, index, residresname\n')
                for i in nNodes_range:
                        Y.write('%s %d %s%d\n'%(col_var_selections[i].name,col_var_selections[i].index+1,col_var_selections[i].resname,col_var_selections[i].resid))
                
                W.write('# Collective variable description: atom name, index, residresname to atom name, index, residresname\n')
		for i in nNodes_range[:-1-ignore_n_nearest_neighbors]:
                        for j in nNodes_range[i+1+ignore_n_nearest_neighbors:]:
                            boolean_matrix[i,j] = True
                            W.write('%s %d %s%d   to    %s %d %s%d\n'%(col_var_selections[i].name,col_var_selections[i].index+1,col_var_selections[i].resname,col_var_selections[i].resid,col_var_selections[j].name,col_var_selections[j].index+1,col_var_selections[j].resname,col_var_selections[j].resid)) 
                            count += 1
        
        print 'Number of nodes:', nNodes,', while skipping ', ignore_n_nearest_neighbors, ' nearest neighbors, creates ', count, 'number of collective variables to be analyzed.'

        # ----------------------------------------
        # TRAJECTORY ANALYSIS
        start = int(start)
        end = int(end)
        print 'Beginning trajectory analysis'
        with open(data_output_filename,'w') as W:
                while start <= end:
	                print 'Loading trajectory ', start, ', located at ', traj_loc%(start)
                        u.load_new(traj_loc%(start))
                        for ts in u.trajectory[::step]:
                                temp_positions = col_var_selections.positions
                                for i in nNodes_range[:-1-ignore_n_nearest_neighbors]
                                        for j in nNodes_range[i+1+ignore_n_nearest_neighbors:]:
                                                dist,dist2 = euclid_dist(temp_positions[i],temp_positions[j])
                                                W.write('%f '%(dist))
                                W.write('\n')
                        start += 1

        print 'Finished analyzing trajectory and outputting raw data. Onto calculating the covariance and correlation matrices for the collective variables analyzed.'
        return data_output_filename, boolean_matrix, count

def calc_collective_variable_boolean_matrix(nColVars,ignore_n_nearest_neighbors=0):
        """
        """
        nNodes = int(np.round(np.max(np.roots([1,-1,1-nColVars*2]))))+ignore_n_nearest_neighbors
        nNodes_range = range(nNodes)
        boolean_matrix = np.full((nNodes,nNodes),False)  # 2D matrix of False values; elements of this matrix will be set to True as we loop over our i,j atom pairs. This will be used later for the plotting of 1D collective variable vectors onto the respective 2D atom pair matrix.
	for i in nNodes_range[:-1-ignore_n_nearest_neighbors]:
                for j in nNodes_range[i+1+ignore_n_nearest_neighbors:]:
                        boolean_matrix[i,j] = True

        print boolean_matrix
        print 'Number of nodes (aka atoms) in the original atom selection is', nNodes, 'with', ignore_n_nearest_neighbors, 'nearest neighbors in this atom selection being ignored in the collective variable space. This creates', nColVars, 'collective variables. Do these numbers check out with the past analysis?'

        return boolean_matrix

