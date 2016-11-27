function write_graph(n_nodes, edges, project_directory)
% Writes graph data to json file
%
% n_nodes              = Number of nodes in the graph
% edges                = Nx2 array where N is the total number of edges. 
%                        The first column gives each edge-source and the 
%                        second column each edge-target.
% project_directory    = Directory to save graph data in
%
%%
    out_nodes = [struct('name',1,'number',1)];
    out_edges = [struct('source',edges(1,1),'target',edges(1,2))];
    for i = 1:n_nodes
        out_nodes(i) = struct('name',i-1,'number',i-1);
    end
    for i = 1:size(edges,1)
        out_edges(i) = struct('source',edges(i,1)-1,'target',edges(i,2)-1);
    end
	out = struct('nodes',out_nodes,'links',out_edges);
    project_directory = strcat(project_directory,'/graph_data.json');
	json.write(out,project_directory);
end