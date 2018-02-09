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
	gene_filter = np.all([mean_filter,var_filter],axis=0)
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
def save_spring_dir(E,D,k,gene_list,project_directory, custom_colors={},cell_groupings={}, use_genes=[], coarse_grain_X=1):
	''' 
	##############################################
	Builds a SPRING project directory and transforms data into SPRING-readable formats
	
	Input (Required)
		E                  = (numpy array) matrix of gene expression. Rows correspond to 
		                     celles and columns correspond to genes. 
		D                  = (numpy array) distance matrix for construction of knn graph.
		                     Any distance matrix can be used as long as higher values
		                     correspond to greater distances. 
		k                  = Number of edges assigned to each node in knn graph
		gene_list          = An ordered list of gene names with length length E.shape[1]
		project_directory  = Path to a directory where SPRING readable files will be 
							 written. The directory does not have to exist before running 
							 this function.
	
		Input (Optional)
		cell_groupings     = Dictionary with one key-value pair for each cell grouping. 
							 The key is the name of the grouping (e.g. "SampleID") and 
							 the value is a list of labels (e.g. ["sample1","sample2"...])
							 If there are N cells total (i.e. E.shape[0] == N), then the
							 list of labels should have N entries. 
		custom_colors      = Dictionary with one key-value pair for each custom color. 
							 The key is the name of the color track and the value is a 
							 list of scalar values (i.e. color intensities). If there are 
							 N cells total (i.e. E.shape[0] == N), then the list of labels 
							 should have N entries. 	
		coarse_grain_X	   = The fold-change reduction in nodes to achieve through graph
							 coarse graining. If 1, then no coarse graining will be 
							 performed
	##############################################
	'''
	# Convert arry data types
	E = np.array(E.tolist()) 
	D = np.array(D.tolist()) 
	
	os.system('mkdir '+project_directory)
	if not project_directory[-1] == '/': project_directory += '/'
	# Build graph
	#print 'Building graph'
	edges = get_knn_edges(D,k)
	

	# Coarse grain 
	if coarse_grain_X < 1:
		edges,E,custom_colors,cell_groupings = coarse_grain(edges,E,k,custom_colors,cell_groupings,coarse_grain_X)
	
	# save genesets
	#print 'Saving gene sets'
	custom_colors['Uniform'] = np.zeros(E.shape[0])
	write_color_tracks(custom_colors, project_directory+'color_data_gene_sets.csv')
	all = []
	
	# save gene colortracks
	#print 'Savng coloring tracks'
	os.system('mkdir '+project_directory+'gene_colors')
	II = len(gene_list) / 50 + 1
	for j in range(50):	
		fname = project_directory+'/gene_colors/color_data_all_genes-'+repr(j)+'.csv'
		if len(use_genes) > 0: all_gene_colors = {g : E[:,i+II*j] for i,g in enumerate(gene_list[II*j:II*(j+1)]) if g in use_genes}
		else: all_gene_colors = {g : E[:,i+II*j] for i,g in enumerate(gene_list[II*j:II*(j+1)]) if np.mean(E[:,i+II*j])>0.05}
		write_color_tracks(all_gene_colors, fname)
		all += all_gene_colors.keys()
	
	# Create and save a dictionary of color profiles to be used by the visualizer
	#print 'Color stats'
	color_stats = {}
	for i in range(E.shape[1]):
		mean = float(np.mean(E[:,i]))
		std = float(np.std(E[:,i]))
		max = float(np.max(E[:,i]))
		centile = float(np.percentile(E[:,i],99.6))
		color_stats[gene_list[i]] = (mean,std,0,max,centile)
	for k,v in custom_colors.items():
		color_stats[k] = (0,1,np.min(v),np.max(v)+.01,np.percentile(v,99))
	json.dump(color_stats,open(project_directory+'/color_stats.json','w'),indent=4, sort_keys=True)
		

	# save cell labels
	#print 'Saving categorical color data'
	categorical_coloring_data = {}
	for k,labels in cell_groupings.items():
		labels = [(l if type(l)==str else repr(l)) for l in labels]
		label_colors = {l:frac_to_hex(float(i)/len(set(labels))) for i,l in enumerate(list(set(labels)))}
		categorical_coloring_data[k] = {'label_colors':label_colors, 'label_list':labels}
	json.dump(categorical_coloring_data,open(project_directory+'/categorical_coloring_data.json','w'),indent=4)
	
	
	#print 'Writing graph'
	nodes = [{'name':i,'number':i} for i in range(E.shape[0])]
	edges = [{'source':i, 'target':j, 'distance':0} for i,j in edges]
	out = {'nodes':nodes,'links':edges}
	open(project_directory+'graph_data.json','w').write(json.dumps(out,indent=4, separators=(',', ': ')))

#========================================================================================#
def coarse_grain(edges,E,k,custom_colors,cell_groupings,coarse_grain_X):
	# make a list of base nodes, give each one a unique index
	# make -1 be a special index that means "unassigned"
	# loop through edges, each one member has a label but the other does not, propogate
	# repeat until all nodes are assigned
	N = len(E)
	node_index = np.ones(N) * -1
	base_filter = np.random.uniform(0,1,len(node_index)) < coarse_grain_X
	node_index[base_filter] = np.arange(np.sum(base_filter))
	
	for count,(i,j) in enumerate(edges):
		if node_index[i] != -1 and node_index[j] == -1:
			node_index[j] = node_index[i]
		if node_index[j] != -1 and node_index[i] == -1:
			node_index[i] = node_index[j]	

	# now make a new list of edges and color vector and coordinates
	new_E = []
	new_custom_colors = {k:[] for k in custom_colors}
	new_cell_groupings = {k:[] for k in cell_groupings}
	for i in range(np.sum(base_filter)):
		new_E.append(np.mean(E[node_index==i,:],axis=0))
		for k in new_custom_colors:
			new_custom_colors[k].append(np.mean(np.array(custom_colors[k])[node_index==i]))
		for k in new_cell_groupings:
			new_cell_groupings[k].append(most_common([cell_groupings[k][ii] for ii in np.nonzero(node_index==i)[0]]))
	new_E = np.array(new_E)
	
	new_edges = {}
	for i,j in edges:
		ne = tuple(sorted([node_index[i],node_index[j]]))
		if not ne in new_edges: new_edges[ne] = 0.
		new_edges[ne] += 1.
	edge_weights = [v for k,v in new_edges.items()]
	cutoff = np.percentile(edge_weights,100./len(edge_weights)*coarse_grain_X*len(edges)*1.5)
	new_edges = [e for e,v in new_edges.items() if v >= cutoff]
	
	return new_edges,new_E,new_custom_colors,new_cell_groupings

#========================================================================================#
def most_common(names):
	counts = {n:0 for n in names}
	for n in names: counts[n] += 1
	highest = np.max([v for k,v in counts.items()])
	for k,v in counts.items():
		if v==highest: return k

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

