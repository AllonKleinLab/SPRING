function write_graph(n_nodes, edges, out_path)
    out_nodes = [struct('name',1,'number',1)];
    out_edges = [struct('source',edges(1,1),'target',edges(1,2))];
    for i = 1:n_nodes
        out_nodes(i) = struct('name',i-1,'number',i-1);
    end
    for i = 1:size(edges,1)
        out_edges(i) = struct('source',edges(i,1)-1,'target',edges(i,2)-1);
    end
	out = struct('nodes',out_nodes,'links',out_edges);
    out_path = strcat(out_path,'/graph_data.json');
	json.write(out,out_path);
end