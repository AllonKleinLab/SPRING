import numpy as np, json, os, matplotlib.pyplot as plt

def get_distance_matrix(M):
	D = np.zeros((M.shape[0],M.shape[0]))
	for i in range(M.shape[0]):
		Mtiled = np.tile(M[i,:][None,:],(M.shape[0],1))
		D[i,:] = np.sqrt(np.sum((Mtiled - M)**2, axis=1))
	return D

def filter_cells(E, min_reads):
	total_counts = np.sum(E,axis=1)
	cell_filter = total_counts >= min_reads
	if np.sum(cell_filter) == 0:
		return None, cell_filter
	else: return E[cell_filter,:],cell_filter
	
def row_normalize(E):
	total_counts = np.sum(E,axis=1)
	tc_tiled = np.tile(total_counts[:,None],(1,E.shape[1]))
	included = np.all(E < tc_tiled * 0.02, axis=0)	
	tc_include = np.sum(E[:,included],axis=1)
	tc_tiled = np.tile(tc_include[:,None],(1,E.shape[1]))+.000001
	return E / tc_tiled * np.mean(total_counts)

def Zscore(E):
	means = np.tile(np.mean(E,axis=0)[None,:],(E.shape[0],1))
	stds = np.tile(np.std(E,axis=0)[None,:],(E.shape[0],1))
	return (E - means) / (stds + .0001)


def filter_genes(E, Ecutoff, Vcutoff):
	mean_filter = np.mean(E,axis=0)> Ecutoff
	var_filter = np.var(E,axis=0) / (np.mean(E,axis=0)+.0001) > Vcutoff
	gene_filter = np.nonzero(np.all([mean_filter,var_filter],axis=0))[0]
	return E[:,gene_filter], gene_filter
	
def get_knn_edges(dmat, k):
	edge_dict = {}
	for i in range(dmat.shape[0]):
		for j in np.nonzero(dmat[i,:] <= sorted(dmat[i,:])[k])[0]:
			if i != j:
				ii,jj = tuple(sorted([i,j]))
				edge_dict[(ii,jj)] = dmat[i,j]
	return edge_dict.keys()

def row_sum_normalize(A):
	print A.shape
	d = np.sum(A,axis=1)
	A = A / np.tile(d[:,None],(1,A.shape[1]))
	return A
	


def write_graph(n_nodes, edges,path):
	nodes = [{'name':i,'number':i} for i in range(n_nodes)]
	edges = [{'source':i, 'target':j, 'distance':0} for i,j in edges]
	out = {'nodes':nodes,'links':edges}
	open(path+'/graph_data.json','w').write(json.dumps(out,indent=4, separators=(',', ': ')))
	
def get_PCA(A,numpc):
	from sklearn.decomposition import PCA
	pca = PCA(n_components=numpc)
	return pca.fit_transform(A)


def sequential_filters(filter_list):
	while len(filter_list) > 1:
		filter_list[0][filter_list[0] > 0] = filter_list[1]
		del filter_list[1]
	return filter_list[0]

def load_coords(DIRECTORY):
	xcoords = {}
	ycoords = {}
	for l in open(DIRECTORY + '/coordinates.txt').read().split('\n')[:-1]:
		l = l.split(',')
		cell_number = int(l[0])
		xcoords[cell_number] = float(l[1])
		ycoords[cell_number] = float(l[2])
	xx,yy = [],[]
	for k in sorted(xcoords.keys()):
		xx += [xcoords[k]]
		yy += [ycoords[k]]
	xx = np.array(xx)
	yy = np.array(yy)
	return xx,-yy


#========================================================================================#
def write_color_tracks(ctracks, fname):
	out = []
	for name,score in ctracks.items():
		line = ','.join([name]+[repr(round(x,1)) for x in score])
		line.replace('.0','.')
		out += [line]
	out = sorted(out,key=lambda x: x.split(',')[0])
	open(fname,'w').write('\n'.join(out))

def frac_to_hex(frac):
	rgb = tuple(np.array(np.array(plt.cm.jet(frac)[:3])*255,dtype=int))
	return '#%02x%02x%02x' % rgb

def save_spring_dir(E,D,k,gene_list,project_directory, geneset_colors={},cell_groupings={}):
	os.system('mkdir '+project_directory)
	if not project_directory[-1] == '/': project_directory += '/'
	# Build graph
	#print 'Building graph'
	edges = get_knn_edges(D,k)

	# save genesets
	#print 'Saving gene sets'
	geneset_colors['Uniform'] = np.zeros(E.shape[0])
	write_color_tracks(geneset_colors, project_directory+'color_data_gene_sets.csv')
	all = []
	
	# save gene colortracks
	#print 'Savng coloring tracks'
	os.system('mkdir '+project_directory+'gene_colors')
	II = len(gene_list) / 50 + 1
	for j in range(50):	
		fname = project_directory+'/gene_colors/color_data_all_genes-'+repr(j)+'.csv'
		all_gene_colors = {g : E[:,i+II*j] for i,g in enumerate(gene_list[II*j:II*(j+1)]) if np.mean(E[:,i+II*j]) > 0.5}
		write_color_tracks(all_gene_colors, fname)
		all += all_gene_colors.keys()
	
	# Create and save a dictionary of color profiles to be used by the visualizer
	#print 'Color stats'
	color_stats = {}
	for i in range(E.shape[1]):
		mean = np.mean(E[:,i])
		std = np.std(E[:,i])
		max = np.max(E[:,i])
		centile = np.percentile(E[:,i],99.6)
		color_stats[gene_list[i]] = (mean,std,0,max,centile)
	for k,v in geneset_colors.items():
		color_stats[k] = (0,1,np.min(v),np.max(v)+.01,np.percentile(v,99))
	json.dump(color_stats,open(project_directory+'/color_stats.json','w'),indent=4, sort_keys=True)
		

	# save cell labels
	#print 'Saving categorical color data'
	categorical_coloring_data = {}
	for k,labels in cell_groupings.items():
		label_colors = {l:frac_to_hex(float(i)/len(set(labels))) for i,l in enumerate(list(set(labels)))}
		categorical_coloring_data[k] = {'label_colors':label_colors, 'label_list':labels}
	json.dump(categorical_coloring_data,open(project_directory+'/categorical_coloring_data.json','w'),indent=4)
	
	
	#print 'Writing graph'
	nodes = [{'name':i,'number':i} for i in range(E.shape[0])]
	edges = [{'source':i, 'target':j, 'distance':0} for i,j in edges]
	out = {'nodes':nodes,'links':edges}
	open(project_directory+'graph_data.json','w').write(json.dumps(out,indent=4, separators=(',', ': ')))
