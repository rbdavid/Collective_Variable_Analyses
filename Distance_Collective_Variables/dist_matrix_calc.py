#!/home/rbdavid/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import os
import importlib
import numpy as np
import matplotlib.pyplot as plt

config_file = sys.argv[1]

flush = sys.stdout.flush
zeros = np.zeros

# ----------------------------------------
# FUNCTIONS:

necessary_parameters = ['output_directory','system_descriptor','distance_functions_file','pca_svd_clustering_functions_file','dist_calc_boolean','step_nFrames','ignore_n_nearest_neighbors','nCluster_list']
all_parameters = ['output_directory','system_descriptor','distance_functions_file','pca_svd_clustering_functions_file','dist_calc_boolean','step_nFrames','ignore_n_nearest_neighbors','nCluster_list','pdb','traj_loc','start','end','atom_selections','user_defined_data_file','corr_calc_boolean','user_defined_mean_vector_file','user_defined_variance_vector_file','user_defined_covariance_matrix_file','user_defined_correlation_matrix_file','equilib_time','delta_t','plotting_boolean','nProjections','figure_format','write_summary']

def config_parser(config_file):	
        """ Function to take config file and create/fill the parameter dictionary (created before function call). 
        
        Usage: 
            parameters = {}     # initialize the dictionary to be filled with keys and values
            config_parser(config_file)

        Arguments:
            config_file: string object that corresponds to the local or global position of the config file to be used for this analysis.

        """
        
        # NECESSARY PARAMETERS ARE INITIALIZED IN DICTIONARY WITH EMPTY STRINGS:
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['pdb'] = None
        parameters['traj_loc'] = None
        parameters['start'] = None
        parameters['end'] = None
        parameters['atom_selections'] = None
        parameters['user_defined_data_file'] = None
        parameters['corr_calc_boolean'] = True
        parameters['user_defined_mean_vector_file'] = None
        parameters['user_defined_variance_vector_file'] = None
        parameters['user_defined_covariance_matrix_file'] = None
        parameters['user_defined_correlation_matrix_file'] = None
        parameters['equilib_time'] = 200    # units of ns
        parameters['delta_t'] = 0.002       # units of ns
        parameters['plotting_boolean'] = False
        parameters['nProjections'] = 2
        parameters['figure_format'] = 'png'
        parameters['write_summary'] = True 

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)

        # TESTING IF ANY PARAMETER HAS BEEN LEFT EMPTY:
        for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

        # TESTING IF PARAMETERS HAVE BEEN DECLARED CORRECTLY AND WITH THEIR DEPENDENT PARAMETERS:
        if parameters['dist_calc_boolean'] == True and (parameters['pdb'] == None or parameters['traj_loc'] == None or parameters['start'] == None or parameters['end'] == None or parameters['atom_selections'] == None):
                print 'You have not read in the necessary parameters to perform the distance matrix analysis on trajectories. Please read in values for pdb, traj_loc, start, end, and atom_selections parameters.'
                sys.exit()

        if parameters['dist_calc_boolean'] == False and parameters['user_defined_data_file'] == None:
                print 'You have not read in the necessary parameters to read in a user defined data file for subsequent pca analysis. Please read in a value for user_defined_data_file.'
                sys.exit()

        if parameters['dist_calc_boolean'] == False and parameters['corr_calc_boolean'] == False and (parameters['user_defined_mean_vector_file'] == None or parameters['user_defined_variance_vector_file'] == None or parameters['user_defined_covariance_matrix_file'] == None or parameters['user_defined_correlation_matrix_file'] == None):
                print 'You have not read in file locations for the necessary parameters to skip the calculation of the mean, variance, covariance, and correlation of the raw collective variable data. Please check the user_defined_mean_vector_file, user_defined_variance_vector_file, user_defined_covariance_matrix_file, and user_defined_correlation_matrix_file.'
                sys.exit()

        if type(parameters['step_nFrames']) != int or parameters['step_nFrames'] <= 0:
                print 'You have not read in an acceptable value for the step_nFrames parameter, which should have an integer value that is equal to or greater than 1. Noninteger values or values of 0 or less do not make logical sense and will raise errors.'
                sys.exit()

        if type(parameters['ignore_n_nearest_neighbors']) != int or parameters['ignore_n_nearest_neighbors'] < 0:
                print 'You have not read in an acceptable value for the ignore_n_nearest_neighbors parameter, which should have an integer value that is equal to or greater than 0. Noninteger values or values less than zero do not make logical sense and will raise errors.'
                sys.exit()

def summary(summary_filename):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs the version of MDAnalysis, how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            summary(summary_filename)

        Arguments:
            summary_filename: string object of the file name to be written that holds the summary information.

        """
	with open(summary_filename,'w') as f:
                if parameters['dist_calc_boolean']:
		        f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
                else:
                        f.write('No distance calculations were performed on MD trajectories. Instead, the user read in a data file.\n')

        	f.write('To recreate this analysis, run this line:\n')
        	for i in range(len(sys.argv)):
        		f.write('%s ' %(sys.argv[i]))
        	f.write('\n\n')
		f.write('Parameters used:\n')
                
                for key, value in parameters.iteritems():
                        if key == '__builtins__':
                                continue
                        if type(value) == int or type(value) == float:
			        f.write("%s = %s\n" %(key,value))
                        else:
			        f.write("%s = '%s'\n" %(key,value))

		f.write('\n')

def main():
        
        system_descriptor = parameters['output_directory'] + parameters['system_descriptor']

        # ----------------------------------------
        # DECIDE IF DISTANCE COLLECTIVE VARIABLES NEED TO BE CALCULATED
        # ----------------------------------------
        if parameters['dist_calc_boolean']:
                print 'Distance calculations need to performed. Beginning trajectory analysis step.'
                data_file_to_be_analyzed, collective_variable_boolean_matrix, nColVars = dist_matrix_calc(parameters['pdb'],parameters['atom_selections'],parameters['traj_loc'],parameters['start'],parameters['end'],system_descriptor,ignore_n_nearest_neighbors=parameters['ignore_n_nearest_neighbors'],step=parameters['step_nFrames'])

        else:
                print 'User is reading in a data file containing collective variable data. Trajectory analysis step is being skipped.'
                data_file_to_be_analyzed = parameters['user_defined_data_file']
                nColVars = int(parameters['nColVars'])
        	# if the trajectory analysis was not performed, we need to prep the 2D boolean_matrix that will be used in the plotting functions
        	collective_variable_boolean_matrix = calc_collective_variable_boolean_matrix(nColVars,ignore_n_nearest_neighbors = parameters['ignore_n_nearest_neighbors'])


        equilib_frame = int(parameters['equilib_time']/(parameters['delta_t']*parameters['step_nFrames']))  # units of frames
        print 'Equilibration time frame number = ', equilib_frame, '. Everything before this frame number will not be analyzed in the subsequent PCA, projection, and clustering analysis.'

        # ----------------------------------------
        # Calculate SVD of raw data
        # ----------------------------------------
	if parameters['svd_boolean']:
        	eigenvector_output_filename = parameters['output_directory'] + '%0'+'%d'%(int(np.log10(nColVars))+1)+'d' + '.' + parameters['system_descriptor'] + '.pca_eigenvector.dat'
        	projected_data_figure_names = parameters['output_directory'] + '%0'+'%d'%(int(np.log10(parameters['nProjections']))+1)+'d' + '.' + parameters['system_descriptor'] + '.projected_data.1d_hist.' + parameters['figure_format']
		
		eigenvector_matrix, projected_data = svd_calc(data_file_to_be_analyzed,system_descriptor,eigenvector_output_filename,equilib_index=equilib_frame,eigenvec_projection_figure_names=projected_data_figure_names,nBins=250)
	
        # ----------------------------------------
        # Or calculate the covariance and correlation matrices of raw data, followed by PCA
        # ----------------------------------------
	else:
		if parameters['corr_calc_boolean']:
        	        data, mean_vector, variance_vector, covariance_matrix, correlation_matrix = covar_corr_matrix_calc(data_file_to_be_analyzed,system_descriptor,equilib_index=equilib_frame,step=1)
        	        print 'Finished calculating the mean, variance vectors, covariance matrix, and correlation matrix. Onto diagonalizing the covariance or correlation matrix and outputting the eigenvalues and eigenvectors.'

        	else:
                	print 'User is reading in data files containing the mean and varaiance vectors as well as the covariance and correlation matrices for the raw data. Calculation of these arrays is being skipped.'
                	data = np.loadtxt(parameters['user_defined_data_file'])
                	mean_vector = np.loadtxt(parameters['user_defined_mean_vector_file'])
                	variance_vector = np.loadtxt(parameters['user_defined_variance_vector_file'])
                	covariance_matrix = np.loadtxt(parameters['user_defined_covariance_matrix_file'])
                	correlation_matrix = np.loadtxt(parameters['user_defined_correlation_matrix_file'])

        	# ----------------------------------------
        	# Plotting Mean, Variance, Covariance, and Correlation results
        	if parameters['plotting_boolean']:
        	        # plotting mean vector as 2d heatmap of mean atom pair distances
        	        mean_vector_heatmap_figure_name = system_descriptor + '.mean_vector.heatmap.' + parameters['figure_format']
        	        plot_vector_as_2dheatmap(mean_vector,collective_variable_boolean_matrix,mean_vector_heatmap_figure_name,cbar_label='Average Distance ($\AA$)', plotting_cmap='Blues')

        	        # plotting variance vector as 2d heatmap of variance of atom pair distances
        	        variance_vector_heatmap_figure_name = system_descriptor + '.var_vector.heatmap.' + parameters['figure_format']
        	        plot_vector_as_2dheatmap(variance_vector,collective_variable_boolean_matrix,variance_vector_heatmap_figure_name,cbar_label='Variance of Distances ($\AA^{2}$)', plotting_cmap='Blues')
        	
        	        # plotting covariance matrix as 2d heatmap 
        	        covar_matrix_heatmap_figure_name = system_descriptor + '.covar_matrix.heatmap.' + parameters['figure_format']
        	        plot_2dmatrix(covariance_matrix,covar_matrix_heatmap_figure_name,cbar_label='Covariance of Collective Variables ($\AA^{2}$)',plotting_cmap='bwr',v_range=[-np.max(covariance_matrix),np.max(covariance_matrix)])
        	    
        	        # plotting covariance matrix as 2d heatmap 
        	        corr_matrix_heatmap_figure_name = system_descriptor + '.corr_matrix.heatmap.' + parameters['figure_format']
        	        plot_2dmatrix(correlation_matrix,corr_matrix_heatmap_figure_name,cbar_label='Correlation of Collective Variables',plotting_cmap='bwr',v_range=[-1.0,1.0])
            
        	# ----------------------------------------
        	# PCA analysis
        	eigenvector_output_filename = parameters['output_directory'] + '%0'+'%d'%(int(np.log10(len(mean_vector)))+1)+'d' + '.' + parameters['system_descriptor'] + '.pca_eigenvector.dat'
        	eigenvector_matrix = pca_calc(covariance_matrix,system_descriptor,eigenvector_output_filename)
        	print 'Finished the principle component analysis on the covariance matrix. Onto projecting the raw data onto the eigenvectors.'

        	# ----------------------------------------
        	# projection analysis
        	zero_padded_string_formatting = '%0'+'%d'%(int(np.log10(parameters['nProjections']))+1)+'d'
        	projected_data_figure_names = parameters['output_directory'] + zero_padded_string_formatting +'.' + parameters['system_descriptor'] + '.projected_data.1d_hist.' + parameters['figure_format']
        	
        	projected_data = data_projection(data,mean_vector,eigenvector_matrix,parameters['nProjections'],system_descriptor,plotting_bool = parameters['plotting_boolean'],eigenvec_projection_figure_names=projected_data_figure_names,nBins=250,test_eigenvec_projections=True)

        # ----------------------------------------
        # Plotting nProjections Eigenvectors onto 2D heatmaps
        if parameters['plotting_boolean']:
                eigenvector_heatmap_figure_name = parameters['output_directory'] + '%0'+'%d'%(int(np.log10(parameters['nProjections']))+1)+'d' + '.' + parameters['system_descriptor'] +  '.eigenvector.heatmap.' + parameters['figure_format']

                for i in range(parameters['nProjections']):
                        plot_vector_as_2dheatmap(np.square(eigenvector_matrix[:,i]),collective_variable_boolean_matrix,eigenvector_heatmap_figure_name%(i),cbar_label='Square of Component Magnitude', plotting_cmap='Blues')
                
                nProjections_eigenvector_heatmap_figure_name = system_descriptor + '.sum_squares_nProjections_eigenvectors.heatmap.' + parameters['figure_format']
                plot_vector_as_2dheatmap(np.sum(np.square(eigenvector_matrix[:,:parameters['nProjections']]),axis=1),collective_variable_boolean_matrix,nProjections_eigenvector_heatmap_figure_name,cbar_label='Sum of Squares', plotting_cmap='Blues')
        
        # ----------------------------------------
        # clustering analysis and plotting
        cluster_labels_output_string = parameters['output_directory'] + '%0'+'%d'%(int(np.log10(np.max(parameters['nCluster_list']))+1))+'d' + '.' + parameters['system_descriptor']
        cluster_figure_names = parameters['output_directory'] + zero_padded_string_formatting + '.' + parameters['system_descriptor'] + '.clustering.' + parameters['figure_format']
        
        kmeans_clustering(projected_data,equilib_frame,parameters['nCluster_list'],system_descriptor,cluster_labels_output_string,cluster_figure_names,step = parameters['step_nFrames'])
        print 'Finished clustering the data. Done with the analyses encoded by this script. How does the data look?'

        # ----------------------------------------
        # SUMMARY OUTPUT 
        if parameters['write_summary']:
                summary_filename = system_descriptor + '.dist_matrix_pca_clustering.summary'
        	summary(summary_filename)

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOADING IN NECESSARY FUNCTIONS FROM MODULE FILES
if parameters['dist_calc_boolean']:
        dist_matrix_calc = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).dist_matrix_calc
else:
        calc_collective_variable_boolean_matrix = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).calc_collective_variable_boolean_matrix

if parameters['svd_boolean']:
	svd_calc = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).svd_calc
elif parameters['corr_calc_boolean']:
        covar_corr_matrix_calc = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).covar_corr_matrix_calc
	pca_calc = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).pca_calc
	data_projection = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).data_projection
	plot_2dmatrix = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).plot_2dmatrix
else:
	pca_calc = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).pca_calc
	data_projection = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).data_projection
	plot_2dmatrix = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).plot_2dmatrix

plot_vector_as_2dheatmap = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).plot_vector_as_2dheatmap

kmeans_clustering = importlib.import_module(parameters['pca_svd_clustering_functions_file'].split('.')[0],package=None).kmeans_clustering

# ----------------------------------------
# CREATING OUTPUT DIRECTORY
if parameters['output_directory'][-1] != os.sep:
        parameters['output_directory'] += os.sep

if os.path.exists(parameters['output_directory']):
        print 'The output directory, ', parameters['output_directory'], ' already exists. Please delete this directory or select a different one for output before proceeding.'
        sys.exit()
else:
        os.mkdir(parameters['output_directory'])

# ----------------------------------------
# MAIN
if __name__ == '__main__':
	main()

