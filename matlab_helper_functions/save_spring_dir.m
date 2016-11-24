function save_spring_dir(E,D,k,gene_list,project_directory, varargin)
    cell_groupings = cell(0);
    custom_colors = cell(0);
    vargs = varargin;
    nargs = length(vargs);
    names = vargs(1:2:nargs);
    values = vargs(2:2:nargs);
  

    %{
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
    %}
end