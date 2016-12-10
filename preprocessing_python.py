import numpy as np, json, os, matplotlib.pyplot as plt
#========================================================================================#
def get_distance_matrix(M):
	'''
	##############################################
	Input
		M = Data matrix. Rows are datapoints (e.g. cell) and columns are features (e.g. genes)
		
	Output (D)
		D = All Pairwise euclidian distances between points in M
	##############################################		
	'''
	D = np.zeros((M.shape[0],M.shape[0]))
	for i in range(M.shape[0]):
		Mtiled = np.tile(M[i,:][None,:],(M.shape[0],1))
		D[i,:] = np.sqrt(np.sum((Mtiled - M)**2, axis=1))
	return D

#========================================================================================#
def filter_cells(E, min_reads):
	'''
	##############################################
	Filter out cells with total UMI count < min_reads
	
	Input
		E         = Expression matrix. Rows correspond to cells and columns to genes
		min_reads = Minimum number of reads required for a cell to survive filtering
		
	Output  (Efiltered, cell_filter)
		Efiltered   = Filtered expression matrix
		cell_filter = Boolean mask that reports filtering. True means that the cell is 
		              kept; False means the cells is removed
	##############################################		
	'''	
	total_counts = np.sum(E,axis=1)
	cell_filter = total_counts >= min_reads
	if np.sum(cell_filter) == 0:
		return None, cell_filter
	else: return E[cell_filter,:],cell_filter

#========================================================================================#	
def row_normalize(E):
	'''
	##############################################
	Normalize so that every cell has the same total read count. Only use genes that make
	up <5% of the total reads in every cell. 
	
	Input
		E = Expression matrix. Rows correspond to cells and columns to genes.
		
	Output (Enormalized)
		Enormalized = Normalized version of the original expression matrix. 
	##############################################
	'''
	total_counts = np.sum(E,axis=1)
	tc_tiled = np.tile(total_counts[:,None],(1,E.shape[1]))
	included = np.all(E < tc_tiled * 0.05, axis=0)	
	tc_include = np.sum(E[:,included],axis=1)
	tc_tiled = np.tile(tc_include[:,None],(1,E.shape[1]))+.000001
	return E / tc_tiled * np.mean(total_counts)

#========================================================================================#
def Zscore(E):
	'''
	##############################################
	Z-score standardize the expression of each gene. 
	
	Input
		E = Expression matrix. Rows correspond to cells and columns to genes.
		
	Output (EZ)
		EZ = Z-score standardized version of the original expression matrix.
	##############################################
	'''
	means = np.tile(np.mean(E,axis=0)[None,:],(E.shape[0],1))
	stds = np.tile(np.std(E,axis=0)[None,:],(E.shape[0],1))
	return (E - means) / (stds + .0001)

#========================================================================================#
def filter_genes(E, Ecutoff, Vcutoff):
	'''
	##############################################
	Filter genes based on expression level and variability. Only keep genes with mean
	expression > Ecutoff and Fano factor > Vcutoff
	
	Input
		E = Expression matrix. Rows correspond to cells and columns to genes.
		Ecutoff = minimum mean expression to keep a gene
		Vcutoff = minimum Fano factor to keep a genes
	
	Output (Efiltered, gene_filter)
		Efiltered   = Filtered expression matrix
		gene_filter = Boolean mask that reports filtering. True means that the gene is 
		              kept; False means the gene is removed		
	##############################################
	'''
	mean_filter = np.mean(E,axis=0)> Ecutoff
	var_filter = np.var(E,axis=0) / (np.mean(E,axis=0)+.0001) > Vcutoff
	gene_filter = np.nonzero(np.all([mean_filter,var_filter],axis=0))[0]
	return E[:,gene_filter], gene_filter

#========================================================================================#	
def get_knn_edges(dmat, k):
	'''
	##############################################
	Calculate knn-graph edges from a distance matrix.
	
	Input
		dmat = Square distance matrix. (dmat)_ij = the distance between i and k
		k    = Number of edges to assign each node (i.e. k in the knn-graph)
		
	Output (edge_list)
		edge_list = A list of unique undirected edges in the knn graph. Each edge comes in 
		            the form of a tuple (i,j) representing an edge between i and j.
	##############################################
	''' 
	edge_dict = {}
	for i in range(dmat.shape[0]):
		for j in np.nonzero(dmat[i,:] <= sorted(dmat[i,:])[k])[0]:
			if i != j:
				ii,jj = tuple(sorted([i,j]))
				edge_dict[(ii,jj)] = dmat[i,j]
	return edge_dict.keys()

#========================================================================================#
def get_PCA(E,numpc):
	'''
	##############################################
	Reduce to the dimension of a data matrix using principal components analysis (PCA)
	
	Input
		E     = Data matrix. Rows are data points (e.g. cells) and columns are features (e.g. genes)
		numpc = Final number of PCs to use
		
	Output (Epca)
		Epca  = PC-reduced data matrix. If E has N rows, then Epca has shape (N,numpc). 
	##############################################
	''' 
	from sklearn.decomposition import PCA
	pca = PCA(n_components=numpc)
	return pca.fit_transform(E)


#========================================================================================#
def save_spring_dir(E,D,k,gene_list,project_directory, custom_colors={},cell_groupings={}):
	pass
#========================================================================================#
def row_sum_normalize(A):
	print A.shape
	d = np.sum(A,axis=1)
	A = A / np.tile(d[:,None],(1,A.shape[1]))
	return A

#========================================================================================#
def write_graph(n_nodes, edges,path):
	nodes = [{'name':i,'number':i} for i in range(n_nodes)]
	edges = [{'source':i, 'target':j, 'distance':0} for i,j in edges]
	out = {'nodes':nodes,'links':edges}
	open(path+'/graph_data.json','w').write(json.dumps(out,indent=4, separators=(',', ': ')))
	
#========================================================================================#
def write_color_tracks(ctracks, fname):
	out = []
	for name,score in ctracks.items():
		line = ','.join([name]+[repr(round(x,1)) for x in score])
		out += [line]
	out = sorted(out,key=lambda x: x.split(',')[0])
	open(fname,'w').write('\n'.join(out))

#========================================================================================#
def frac_to_hex(frac):
	rgb = tuple(np.array(np.array(plt.cm.jet(frac)[:3])*255,dtype=int))
	return '#%02x%02x%02x' % rgb
