#!/home/rbdavid/bin/python

import sys 
import numpy as np
import matplotlib.pyplot as plt

nodes_string_txt_file = open(sys.argv[1],'r')
selection_string_txt_file = open(sys.argv[2], 'r')
nEigenvectors = int(sys.argv[3])
eigenvector_file_names = sys.argv[4]
nColVars = int(sys.argv[5])
ignore_n_nearest_neighbors = int(sys.argv[6])

nNodes = int(np.round(np.max(np.roots([1,-1,1-nColVars*2]))))+ignore_n_nearest_neighbors
print nNodes
nNodes_range = range(nNodes)

node_strings = []
for line in node_string_txt_file:
        if line.startswith('#'):
                continue
        else: 
                node_strings.append(line.rstrip('\n'))

if len(node_strings) != nNodes:
        print 'number of nodes does not match...'

selection_strings = []
for line in selection_string_txt_file:
        if line.startswith('#'):
                continue
        else: 
                selection_strings.append(line.rstrip('\n'))

legend_labels = []
for i in range(nEigenvectors):
        print 'Analyzing eigenvector', i

        eigenvector_data = np.loadtxt(eigenvector_file_names %(i))
        sq_eigenvector_data = np.square(eigenvector_data)
        
        mean = np.mean(sq_eigenvector_data)
        stdev = np.std(sq_eigenvector_data)
        print 'average square of eigenvector component = ', mean, ' standard deviation = ',stdev
        
        indices = np.argwhere(sq_eigenvector_data > mean+3*stdev)
        num_components = len(indices)
        indices.shape = (num_components,)
        
        important_eigenvector_components = sq_eigenvector_data[indices]
        idx = important_eigenvector_components.argsort()[::-1]
        important_eigenvector_components = important_eigenvector_components[idx]
        indices = indices[idx]
        for j in indices:
                print '%.4f' %(sq_eigenvector_data[j]), selection_strings[j]
        print '\n'
        
        histogram = np.zeros(nNodes)
        count = 0
        for a in nNodes_range[:-1-ignore_n_nearest_neighbors]:
                for b in nNodes_range[a+1+ignore_n_nearest_neighbors:]:
                        if sq_eigenvector_data[count] > mean+3*stdev:
                                histogram[a] += 1
                                histogram[b] += 1
                        count += 1

        plt.plot(histogram,ls='-',lw=1.0)
        legend_labels.append('Vec %d'%(i))

        mean = np.mean(histogram)
        stdev = np.std(histogram)

        indices = np.argwhere(histogram > mean+3*stdev)
        num_nodes = len(indices)
        indices.shape = (num_nodes,)
        
        important_nodes = histogram[indices]
        idx = important_nodes.argsort()[::-1]
        important_nodes = important_nodes[idx]
        indices = indices[idx]
        for j in indices:
                print '%d' %(histogram[j]), node_strings[j]
        print '\n'

plt.xlim((-0.25,nNodes+0.25))
plt.xlabel('Node Index',size=14)
plt.ylabel('Counts',size=14)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.legend(legend_labels)
plt.savefig('Eigenvector_component_analysis.png',dpi=600,transparent=True)
plt.close()

