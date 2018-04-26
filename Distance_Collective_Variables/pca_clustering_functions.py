
# USAGE:
# from pca_clustering_functions import *

# PREAMBLE:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter, MultipleLocator
from sklearn.cluster import KMeans

sqrt = np.sqrt

# SUBROUTINES/FUNCTIONS:

def covar_corr_matrix_calc(data_file,system_descriptor,equilib_index=0,step=1):
        """ Calculates the mean and variance vectors as well as the covariance and correlation matrices for a data found in a data file. The data file is assumed to have a specific organization: rows correspond to nInstances of multivariate data; the nColVars columns correspond to the multivariate data over all nInstances. This analysis code calculates the mean, variance, covariance, and correlation arrays for the nColVars data (columns).

        Usage:
                data, mean_vector, var_vector, covariance_matrix, correlation_matrix = covar_corr_matrix_calc(data_file,system_descriptor, equilib_index = equilib_frame, step=step_nFrames)

        Arguments:
                data_file: string object that corresponds to the data file containing the dataset to be analyzed. Rows correspond to nInstances. Columns correspond to individual/specific/consistent collective variables.
                system_descriptor: string object that will be used in naming output files created during this analysis.
                equilib_index: OPTIONAL; DEFAULT VALUE = 0; int object with a minimum value of 0; rows with indices less than equilib_index will be ignored in the analysis. 
                step: OPTIONAL; DEFAULT VALUE = 1; POTENTIAL BUGS; int object used to step through the dataset loaded in from the data_file.

### COMMENTS
        Output:
                data:
                mean_vector:
                var_vector:
                covariance_matrix:
                correlation_matrix:

        """

        data = np.loadtxt(data_file)[equilib_index::step]    # assumes data file is formatted with timestep data in each row, a column corresponds to a individual/specific/consistent collective variable.
        nInstances = len(data)
        nColVars = len(data[0])
        print 'nInstances (aka nSteps) being analyzed: ', nInstances, ' each with ',nColVars, ' nColVars.'
        
        mean_vector = np.mean(data,axis=0)
        var_vector = np.mean(np.square(data),axis=0) - np.square(mean_vector)

        covariance_matrix = np.zeros((nColVars,nColVars),dtype=np.float64)
        for ts in range(nInstances):
                for i in range(nColVars):
                        for j in range(i,nColVars):
                                covariance_matrix[i,j] += data[ts,i]*data[ts,j]
        
        covariance_matrix /= nInstances

        for i in range(nColVars):
                for j in range(i,nColVars):
                        covariance_matrix[i,j] -= mean_vector[i]*mean_vector[j]
                        
        correlation_matrix = np.copy(covariance_matrix)
        for i in range(nColVars):
                for j in range(i,nColVars):
                        correlation_matrix[i,j] /= sqrt(var_vector[i]*var_vector[j]) 
                        covariance_matrix[j,i] = covariance_matrix[i,j]
                        correlation_matrix[j,i] = correlation_matrix[i,j]

        # SAVING RESULTS OUT TO FILE
        np.savetxt(system_descriptor+'.mean_vector.dat',mean_vector,delimiter='   ')
        np.savetxt(system_descriptor+'.variance_vector.dat',var_vector,delimiter='   ')
        np.savetxt(system_descriptor+'.covariance_vector.dat',covariance_matrix,delimiter='   ')
        np.savetxt(system_descriptor+'.correlation_vector.dat',correlation_matrix,delimiter='   ')

        # RETURNING RESULTS TO MAIN 
        return data, mean_vector, var_vector, covariance_matrix, correlation_matrix

def plot_vector_as_2dheatmap(vector,two_dimensional_boolean_matrix,figure_name,cbar_label='',plotting_cmap='bwr',v_range=None,minor_ticks=1,major_ticks=10):
        """
        """
        vector = np.copy(vector)
        two_dimensional_boolean_matrix = np.copy(two_dimensional_boolean_matrix)
        nNodes = len(two_dimensional_boolean_matrix)
        two_dimensional_matrix_from_vector = np.zeros((nNodes,nNodes),dtype=np.float32)
        counter = 0
        for i in range(nNodes):
                for j in range(nNodes):
                        if two_dimensional_boolean_matrix[i][j]:
                                two_dimensional_matrix_from_vector[i][j] = vector[counter]
                                two_dimensional_matrix_from_vector[j][i] = vector[counter]
                                counter += 1

        node_range = range(nNodes+1)
        fig, ax = plt.subplots()
        ax.tick_params(which='major',length=6,width=2)
        ax.tick_params(which='minor',length=3,width=1)
        ax.xaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.xaxis.set_major_locator(MultipleLocator(major_ticks))
        ax.yaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.yaxis.set_major_locator(MultipleLocator(major_ticks))

        if v_range != None:
                temp = plt.pcolormesh(node_range,node_range,two_dimensional_matrix_from_vector,cmap=plotting_cmap,vmin=v_range[0],vmax=v_range[1])
        else:
                temp = plt.pcolormesh(node_range,node_range,two_dimensional_matrix_from_vector,cmap=plotting_cmap)
        cb1 = plt.colorbar()
        cb1.set_label(r'%s'%(cbar_label))

        xlabels = [str(int(x)) for x in temp.axes.get_xticks()[:]]
        ylabels = [str(int(y)) for y in temp.axes.get_yticks()[:]]
        temp.axes.set_xticks(temp.axes.get_xticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
        temp.axes.set_yticks(temp.axes.get_yticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)
        temp.axes.set_xticklabels(xlabels)
        temp.axes.set_yticklabels(ylabels)

        plt.xlim((-0.5,nNodes+0.5))
        plt.ylim((-0.5,nNodes+0.5))
        plt.xlabel('Atom index in atom selection',size=14)
        plt.ylabel('Atom index in atom selection',size=14)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig(figure_name,dpi=600,transparent=True)
        plt.close()

def plot_2dmatrix(square_matrix,figure_name,cbar_label='',plotting_cmap='bwr',v_range=None,minor_ticks=10,major_ticks=100):
        """
        """
        nNodes = len(square_matrix)
        node_range = range(nNodes+1)
        fig, ax = plt.subplots()
        ax.tick_params(which='major',length=6,width=2)
        ax.tick_params(which='minor',length=3,width=1)
        ax.xaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.xaxis.set_major_locator(MultipleLocator(major_ticks))
        ax.yaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.yaxis.set_major_locator(MultipleLocator(major_ticks))
        
        if v_range != None:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap,vmin=v_range[0],vmax=v_range[1])
        else:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap)
        cb1 = plt.colorbar()
        cb1.set_label(r'%s'%(cbar_label))

        xlabels = [str(int(x)) for x in temp.axes.get_xticks()[:]]
        ylabels = [str(int(y)) for y in temp.axes.get_yticks()[:]]
        temp.axes.set_xticks(temp.axes.get_xticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
        temp.axes.set_yticks(temp.axes.get_yticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)

        plt.xlim((-0.5,nNodes+0.5))
        plt.ylim((-0.5,nNodes+0.5))
        plt.xlabel('Collective Variable Index',size=14)
        plt.ylabel('Collective Variable Index',size=14)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig(figure_name,dpi=600,transparent=True)
        plt.close()

def pca_calc(square_matrix,system_descriptor,eigenvector_output_filename):
        """
        """
        square_matrix = np.copy(square_matrix)
        # CALCULATE THE EIGENVALUES AND EIGENVECTORS OF THE SQUARE MATRIX
        eigval,eigvec = np.linalg.eig(square_matrix)    # NOTE: EIGENVEC IS ORGANIZED WHERE COMPONENTS OF ONE, INDIVIDUAL EIGENVECTOR ARE STORED IN THE COLUMN (SECOND INDEX)
        # SORT THE EIGENVALUES AND EIGENVECTORS FROM LARGEST VALUED TO SMALLEST VALUED
        idx = eigval.argsort()[::-1]    # get arrangement of eigval indices that sort from largest value to smallest value
        eigval = eigval[idx]    # rearrange eigval
        eigvec = eigvec[:,idx]  # rearrange the columns of eigvec
        
        # ANALYZE THE EIGENVALUES AND CUMULATIVE EIGENVALUE TO HELP THE USER DECIDE THE NUMBER OF EIGENVECTORS THAT ADEQUATELY DESCRIBE THE DATASET
        nVec = len(eigvec)
        cumulative_eigval = np.zeros(nVec,dtype=np.float64)
        total_eigval = 0
        for i in range(nVec):
                total_eigval += eigval[i]
                cumulative_eigval[i] = total_eigval

        # OUTPUT EIGENVALUES AND EIGENVECTORS TO FILE
        with open(system_descriptor+'.pca_eigenvalues.dat','w') as f:
                f.write('# Eigenvalue   Frac_Total   Cumulative   Frac_Cumulative\n')
                for i in range(nVec):
                        f.write('%f   %f   %f   %f\n' %(eigval[i],eigval[i]/total_eigval,cumulative_eigval[i],cumulative_eigval[i]/total_eigval))
                        np.savetxt(eigenvector_output_filename%(i),eigvec[:,i],fmt='%f')
        
        # RETURN THE EIGENVEC ARRAY
        return eigvec

def data_projection(data,mean_vector,var_vector,eigvec,nProjections,system_descriptor,standardize=False,plotting_bool=True,eigenvec_projection_figure_names='%d.projected_data.1d_hist.png',nBins=100,test_eigenvec_projections=True):
        """
        """

        data -= mean_vector
        if standardize:
                data /= np.sqrt(var_vector)
        projection_data = np.zeros((len(data),nProjections),dtype=np.float64)
        for i in range(nProjections):
                projection_data[:,i] = np.dot(data,eigvec[:,i])

                if plotting_bool:
                        events,edges,patches = plt.hist(projection_data[:,i],bins=nBins,histtype='bar',normed=True)
	                plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
                        plt.xlabel('Data projected onto Eigenvector %d'%(i))
                        plt.ylabel('Probability Density')
                        plt.savefig(eigenvec_projection_figure_names %(i),dpi=600,transparent=True)
                        plt.close()
        
        if test_eigenvec_projections:
                print 'The data has been mean centered :. the average of the projected data should be zero. The dot product of the eigenvectors should be zero. The slope of projected data should be close to zero as well.'
                for i in range(nProjections-1):
                        print 'Eigenvec', i, ' average = ', np.mean(projection_data[:,i])
                        for j in range(i+1,nProjections):
                                print 'Eigenvec', i, 'and eigenvec', j,': dot product = ', np.dot(eigvec[:,i],eigvec[:,j]), 'slope, intercept of linear least squares:', np.polyfit(projection_data[:,i],projection_data[:,j],deg=1)
                print 'Does everything look the way it should?'
        
        np.savetxt(system_descriptor+'.projected_data.dat',projection_data)

        return projection_data

def kmeans_clustering(projection_data,equilib_frame,nClusters_list,system_descriptor,cluster_labels_output_string,cluster_figure_names, step = 1):
        """
        """

### COMMENTS

        nFrames = len(projection_data)
        nColVars = len(projection_data[0])
        frame_numbers = [equilib_frame + i*step for i in range(nFrames)]
        elbow_w_list = []
        with open(system_descriptor+'.cluster_validation.dat','w') as W, open(system_descriptor+'.cluster_centers.dat','w') as X:
                X.write('# cluster number, minimum euclid dist to centroid, index of frame closest to centroid, number of frames in cluster\n')
                for nClusters in nClusters_list:
                        X.write('Number of clusters: %d\n'%(nClusters))
                	
                        #-------------------------------------
                	# Initialize the clusterer with desired kwargs
                        #-------------------------------------
                        clusterer = KMeans(n_clusters=nClusters,init='k-means++',n_init=100,max_iter=1000,tol=0.0000001,precompute_distances='auto',verbose=0,n_jobs=1)
                        cluster_labels = clusterer.fit_predict(projection_data)   # returns labels for each data point
                        np.savetxt(cluster_labels_output_string %(nClusters) + '.cluster_labels.dat',np.c_[frame_numbers,cluster_labels],fmt='%d   %d',header='Frame Number   Cluster Number')
                        
                        cluster_centers = clusterer.cluster_centers_    # returns n_clusters x n_features matrix, each row corresponds to the data of a single MD frame
                        print 'For nClusters =', nClusters,', finished labeling the full dataset into specific clusters. Outputting the frames closest to cluster centers, followed by calculating elbow data for the full data set.'
        
                        #-------------------------------------
                        # determining the trajectory frame corresponding to the data point closest to the centroids of each cluster
                        #-------------------------------------
                        for i in range(nClusters):
                                min_eu_dist2 = nColVars*10
                                min_index = 0.
                                for j in range(nFrames):
                                        if cluster_labels[j] == i:
                                                eu_dist2 = np.sum(np.square(projection_data[j] - cluster_centers[i]))
                                                if eu_dist2 < min_eu_dist2:
                                                        min_eu_dist2 = eu_dist2
                                                        min_index = j
                                X.write('%d   %f   %d   %d\n'%(i,np.sqrt(min_eu_dist2),int(min_index+equilib_frame),len(projection_data[cluster_labels == i])))
                       
                        #-------------------------------------
                        # Calculating the euclidian distance btw each centroid, using the projection space
                        #-------------------------------------
                        cluster_dist_matrix = np.zeros((nClusters,nClusters),dtype=np.float64)
                        for i in range(nClusters):
                                for j in range(i+1,nClusters):
                                        cluster_dist_matrix[i,j] = np.sqrt(np.sum(np.square(cluster_centers[i] - cluster_centers[j])))
                        np.savetxt(cluster_labels_output_string %(nClusters) + '.cluster_centers_dist_matrix.dat',cluster_dist_matrix)
        
                        #-------------------------------------
                        # CALCULATE THE TOTAL INTRA-CLUSTER VARIATION DATA FOR THE CLUSTERING; USED TO CREATE AN ELBOW PLOT TO DETERMINE THE APPROPRIATE NUMBER OF CLUSTERS.
                        #-------------------------------------
                        intra_cluster_variation = np.zeros(nClusters,dtype=np.float64)
                        for i in range(nClusters):
                                intra_cluster_variation[i] = np.sum(np.square(projection_data[cluster_labels == i] - cluster_centers[i]))    #   /len(projection_data[cluster_labels == i])
                        print intra_cluster_variation
                        elbow_w = np.sum(intra_cluster_variation)
                        elbow_w_list.append(elbow_w)
                        W.write('For n_clusters = %d, the total intra-cluster variation is: %f\n' %(nClusters,elbow_w))
                        print 'For nClusters =', nClusters,', finished calculating elbow data for the full data set. Moving onto plotting projected data (colored by cluster labeling).'
        
                        #-------------------------------------
                        # PLOTTING OF PROJECTED DATA
                        #-------------------------------------
                        colors = cm.spectral(cluster_labels.astype(float) / nClusters)
                        plt.scatter(projection_data[:,0],projection_data[:,1],marker='.',s=5,lw=0,alpha=0.7,c=colors)
                        plt.plot(cluster_centers[:,0],cluster_centers[:,1],marker='X',mec='k',mfc='r',ms=15,ls='None')
                        for i in range(nClusters):
                                plt.text(cluster_centers[i,0],cluster_centers[i,1],str(i),color='r',fontsize=15)
        	        plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
                        plt.xlabel("Mean-centered, Standardized Data Projected onto PC 0" )
                        plt.ylabel("Mean-centered, Standardized Data Projected onto PC 1" )
                        plt.savefig(cluster_figure_names %(nClusters),transparent=True,dpi=600)
                        plt.close()
        
        #-------------------------------------
        # PLOTTING OF TOTAL INTRA-CLUSTER VARIANCE RESULTS FOR THE RANGE OF N_CLUSTERS DATA
        #-------------------------------------
        plt.plot(nClusters_list, elbow_w_list,'ko')
        plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
        plt.xlabel('Number of Clusters',size=14)
        plt.ylabel('Total Intra-Cluster Variation',size=14)
        plt.savefig(system_descriptor + '.elbow_plot.png',transparent=True,dpi=600)
        plt.close()

